package BioFuse::BioInfo::FASTA::GetMutFA;

use strict;
use warnings;
use Getopt::Long;
use List::Util qw/ min max any first /;
use Data::Dumper;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::Interval qw/ Get_Two_Seg_Olen intersect /;
use BioFuse::LoadOn;
use BioFuse::BioInfo::FASTA qw/ read_fasta_file /;
use BioFuse::BioInfo::BED qw/ read_bed_file /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              GetMutFasta
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::FASTA::GetMutFA';
#----- version --------
$VERSION = "0.05";
$DATE = '2019-10-22';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        GetMutFasta
                        load_mutations
                        merge_mutations
                        alertDiffCN
                        get_mut_refseq
                        add_mut_on_chr
                        add_single_mut
                     /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} mut_fa <[Options]>

     Options:
         -f  [s]  reference genome fasta file. <required>
         -m  [s]  mutation list. <required>
         -o  [s]  output mutated fasta file. <required>
         -s  [i]  flanking size around mutation, minimum:200. [1000]
         -d  [i]  minimum distance allowed between neighbor mutations. [0, disabled]
                  check snv/ins/del
         -da [i]  action when found distance less than '-d'. [x]
                  'x': alert and eXit; 's': Skip latter mutation and continue.
         -b  [s]  region BED file to filter mutations.
         -one     given region is one-based coordinate, not BED format. [disabled]
         -h       show this help

     Version:
        $VERSION at $DATE

     Author:
        $AUTHOR ($EMAIL)
 \n";
}

#--- load variant of this module to public variant (V_Href in LoadOn.pm) ---
sub Load_moduleVar_to_pubVarPool{
    $V_Href->{ $_->[0] } = $_->[1] for
        map {
            if( !exists $V_Href->{$_->[0]} ){
                ( $_ );
            }
            else{
                warn_and_exit "<ERROR>\tkey $_->[0] is already in V_Href!\n";
            }
        }
        (
            # input/output
            ## use 'whole_genome' in BioFuse::LoadOn
            [ mut_fa => undef ],
            [ mut_list => undef ],
            [ reg_bed => undef ],
            # option
            [ flankSize => 1000 ],
            [ oneBase => 0 ], # one-start region, not BED
            [ minMutDist => 0 ], # neighbor mutations' min distance
            [ minMutDistAct => 'x' ],
            # container
            [ mut => {} ],
            [ reg => {} ],
            # pre-sets
            [ mtcd => {snv=>1,ins=>1,del=>1,mnp=>1} ], # mut type to check distance
            # list to abs-path
            [ ToAbsPath_Aref => [ ['mut_fa'],
                                  ['mut_list'],
                                  ['reg_bed'],
                                  ['whole_genome']  ] ]
        );
}

#--- get options from command line ---
sub Get_Cmd_Options{
    # get options
    GetOptions(
        # input/output
        "-f:s"  => \$V_Href->{whole_genome},
        "-o:s"  => \$V_Href->{mut_fa},
        "-m:s"  => \$V_Href->{mut_list},
        "-b:s"  => \$V_Href->{reg_bed},
        # option
        "-s:i"  => \$V_Href->{flankSize},
        "-one"  => \$V_Href->{oneBase},
        "-d:i"  => \$V_Href->{minMutDist},
        "-da:s" => \$V_Href->{minMutDistAct},
        # help
        "-h|help"   => \$V_Href->{HELP},
        # for debug
        "-debug"    => \$V_Href->{in_debug} # hidden option
    );
}

#--- test para and alert ---
sub para_alert{
    return  (   $V_Href->{HELP}
             || !file_exist(filePath=>$V_Href->{whole_genome})
             || !file_exist(filePath=>$V_Href->{mut_list})
             || $V_Href->{flankSize} < 200
             || !defined $V_Href->{mut_fa}
             || $V_Href->{minMutDist} < 0
            );
}

#--- get fasta file from reference genome based on given mutations ---
sub GetMutFasta{
    # load mutations
    &load_mutations;

    # filter by region if provided
    &filer_mutations;

    # merge mutations with flank size
    &merge_mutations;

    # modify seq from ref-fasta based on mutations and output
    &get_mut_refseq;
}

#--- load mutation list ---
#chr stP edP hapInfo gene
## hapInfo
## ref+ref  h1:ref,A,1;h2:ref,A,1
## snv+ref  h1:ref,A,1;h2:snv,G,1
## ins+ref  h1:ref,A,1;h2:ins,+CC,1
## del+ref  h1:ref,A,1;h2:ins,-3,1  or  h1:ref,A,1;h2:ins,-AGA,1
## cnv      single locus: change the last integer of each hap
## cnv      region locus: h1:cnv,-,2;h2:cnv:-:3
sub load_mutations{
    open (ML, Try_GZ_Read($V_Href->{mut_list})) || die "fail read mutation list: $!\n";
    # theme
    (my $theme_line = lc(<ML>)) =~ s/^#//;
    my @theme_tag = split /\s+/, $theme_line;
    while(<ML>){
        next if /^\#/;
        my @info = split;
        my %mut = map{ ($theme_tag[$_],$info[$_]) } (0 .. $#theme_tag);
        my $mut_Hf = { stp => $mut{stp}, edp => $mut{edp}, gene => $mut{gene}||'N/A' };
        # hapInfo
        my %hapInfo = map { my ($hapNO, $type, $allele, $cn) = split /[:,]/;
                            ($hapNO => {type => $type, allele => $allele, cn => $cn});
                          } split /;/, $mut{hapinfo};
        $mut_Hf->{hapinfo} = \%hapInfo;
        # record mut
        push @{$V_Href->{mut}->{$mut{chr}}}, $mut_Hf;
    }
    close ML;
    # inform
    stout_and_sterr "[INFO]\tload mutation list.\n";
}

#--- only keep mutations in given region ---
sub filer_mutations{
    # not provided
    return unless file_exist(filePath=>$V_Href->{reg_bed});

    # load bed region, do merge
    $V_Href->{reg} = read_bed_file(bedFile => $V_Href->{reg_bed}, nonName => 1, oneBase => $V_Href->{oneBase});
    # inform
    stout_and_sterr "[INFO]\tload required region.\n";

    # keep mut/region in required regions
    for my $chr (sort keys %{$V_Href->{mut}}){
        if(!exists $V_Href->{reg}->{$chr}){
            delete $V_Href->{mut}->{$chr};
        }
        else{
            my $chrMutAf = $V_Href->{mut}->{$chr};
            @$chrMutAf = sort {$a->{stp}<=>$b->{stp}} @$chrMutAf;
            my $chrRegAf = $V_Href->{reg}->{$chr};
            # check mut one by one
            my $shift = 0;
            for my $i (0 .. scalar(@$chrMutAf)-1){
                my $mut = $chrMutAf->[$i+$shift];
                # till find the first lagging region
                shift @$chrRegAf while scalar(@$chrRegAf) != 0 && $chrRegAf->[0]->[1] < $mut->{stp};
                # no more required regions
                if(scalar(@$chrRegAf) == 0){
                    splice @$chrMutAf, $i+$shift; # sweep remain mutations
                    last;
                }
                # try to get overlaps region (or)
                my $orAf = intersect(itvAfListAf => [[[$mut->{stp}, $mut->{edp}]], $chrRegAf], skipMerge => 1);
                # deal with mutations
                my @orMut = map {{stp=>$_->[0], edp=>$_->[1], gene=>$mut->{gene}, hapinfo=>$mut->{hapinfo}}} @$orAf;
                splice @$chrMutAf, $i+$shift, 1, @orMut; # replace or discard
                # update shift
                $shift += scalar(@$orAf) - 1;
            }
        }
    }
    # inform
    stout_and_sterr "[INFO]\tfilter mutation with required region.\n";
}

#--- merge mutations with flank size ---
sub merge_mutations{
    for my $chr (sort keys %{$V_Href->{mut}}){
        @{$V_Href->{mut}->{$chr}} = sort {$a->{stp}<=>$b->{stp}} @{$V_Href->{mut}->{$chr}};
        my $Hf_1st = shift @{$V_Href->{mut}->{$chr}};
        my @merged = ({ stp => max($Hf_1st->{stp} - $V_Href->{flankSize}, 1),
                        edp => $Hf_1st->{edp} + $V_Href->{flankSize},
                        mut => [$Hf_1st]
                       });
        for my $mut_Hf (@{$V_Href->{mut}->{$chr}}){
            my $mfstp = max($mut_Hf->{stp} - $V_Href->{flankSize}, 1);
            my $mfedp = $mut_Hf->{edp} + $V_Href->{flankSize};
            my $ovlen = Get_Two_Seg_Olen($merged[-1]->{stp}, $merged[-1]->{edp}, $mfstp, $mfedp);
            # alert neighbor distance
            if($V_Href->{minMutDist} > 0){
                my $formerMut = first { any {exists $V_Href->{mtcd}->{$_->{type}}} values %{$_->{hapinfo}} } reverse @{$merged[-1]->{mut}};
                if(    defined $formerMut
                    && (any{exists $V_Href->{mtcd}->{$_->{type}}} values %{$mut_Hf->{hapinfo}})
                    && &test_mut_dist(chr => $chr, a_mut => $formerMut, b_mut => $mut_Hf) # return 1, means too near
                ){
                    next;
                }
            }
            # merge or seperate
            if($ovlen){ # merge!
                # $merged[-1]->{stp} = min($merged[-1]->{stp}, $mfstp); # useless
                $merged[-1]->{edp} = max($merged[-1]->{edp}, $mfedp);
                push @{$merged[-1]->{mut}}, $mut_Hf;
            }
            else{ # seperate
                &alertDiffCN(mut_Af => $merged[-1]->{mut});
                push @merged, { stp => $mfstp, edp => $mfedp, mut => [$mut_Hf] };
            }
        }
        &alertDiffCN(mut_Af => $merged[-1]->{mut});
        # update
        $V_Href->{mut}->{$chr} = \@merged;
        # inform
        stout_and_sterr "[INFO]\tgroup mutations on $chr.\n";
    }
}

#--- test distance of given mutations, obey the mut type ---
sub test_mut_dist{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $chr = $parm{chr};
    my $a_mut = $parm{a_mut};
    my $b_mut = $parm{b_mut};
    if(max($a_mut->{stp}, $b_mut->{stp}) - min($a_mut->{edp}, $b_mut->{edp}) < $V_Href->{minMutDist}){
        my $warn =  "$chr neighbor mutations distance less than required ($V_Href->{minMutDist}).\n"
                   .Data::Dumper->Dump([$a_mut,$b_mut],[qw/a_mut b_mut/]);
        if($V_Href->{minMutDistAct} eq 'x'){
            warn_and_exit $warn;
        }
        else{
            warn $warn."skip the latter one\n";
            return 1;
        }
    }
    else{
        return 0;
    }
}

#--- test regional merged mut have different hap cn ---
sub alertDiffCN{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $mut_Af = $parm{mut_Af};
    for my $i (1 .. scalar(@$mut_Af)-1){
        my $mut_i_Hf = $mut_Af->[$i];
        my $mut_j_Hf = $mut_Af->[$i-1];
        if(   scalar keys %{$mut_i_Hf->{hapinfo}}
           != scalar keys %{$mut_j_Hf->{hapinfo}}
        ){
            warn_and_exit  "<ERROR>\tmerged mutations have different amount of haplotypes.\n"
                          .Data::Dumper->Dump([$mut_i_Hf, $mut_j_Hf], [qw(mut_i mut_j)]);
        }
        for my $hapNO (sort keys %{$mut_i_Hf->{hapinfo}}){
            if(   !exists $mut_j_Hf->{hapinfo}->{$hapNO}
               || $mut_i_Hf->{hapinfo}->{$hapNO}->{cn} != $mut_j_Hf->{hapinfo}->{$hapNO}->{cn}
            ){
                warn_and_exit  "<ERROR>\tmerged mutations have different haplotype $hapNO.\n"
                              .Data::Dumper->Dump([$mut_i_Hf, $mut_j_Hf], [qw(mut_i mut_j)]);
            }
        }
    }
}

#--- modify seq from ref-fasta based on mutations and output ---
sub get_mut_refseq{
    # read fasta file
    open (my $mfafh, Try_GZ_Write($V_Href->{mut_fa})) || die "fail write output mut_fa file: $!\n";
    read_fasta_file( FaFile => $V_Href->{whole_genome},
                     subrtRef => \&add_mut_on_chr,
                     subrtParmAref => [mfafh => $mfafh, mut_Hf => $V_Href->{mut}]
                   );
    close $mfafh;
}

#--- modify seq of given chr based on mutations ---
sub add_mut_on_chr{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $segName = $parm{segName};
    my $segSeq_Sref = $parm{segSeq_Sref};
    my $mfafh = $parm{mfafh};
    my $mut_Hf = $parm{mut_Hf};

    return unless exists $mut_Hf->{$segName};

    my $linebase = 50;
    for my $grp_Hf (@{$mut_Hf->{$segName}}){
        my $refseq = substr($$segSeq_Sref, $grp_Hf->{stp}-1, $grp_Hf->{edp}-$grp_Hf->{stp}+1);
        my %seq;
        my $mut_Af = $grp_Hf->{mut};
        my %hap_cn;
        # firstly SNV, then INDEL
        for my $t ('snv', '(ins)|(del)'){
            for my $mut (sort {$b->{stp}<=>$a->{stp}} @$mut_Af){
                my $pidx = $mut->{stp} - $grp_Hf->{stp};
                for my $hapNO (sort keys %{$mut->{hapinfo}}){
                    my $hap_Hf = $mut->{hapinfo}->{$hapNO};
                    $seq{$hapNO} = $refseq unless exists $seq{$hapNO};
                    $hap_cn{$hapNO} = $hap_Hf->{cn}; # update
                    next if $hap_Hf->{cn} == 0 || $hap_Hf->{type} !~ /^$t$/i;
                    &add_single_mut(type=>$hap_Hf->{type}, alt=>$hap_Hf->{allele}, pidx=>$pidx, seqSf=>\$seq{$hapNO});
                }
            }
        }
        # output
        for my $hapNO (sort keys %seq){
            for my $cn (1 .. $hap_cn{$hapNO}){
                print {$mfafh} ">$segName:$grp_Hf->{stp}-$grp_Hf->{edp}:$hapNO:cn$cn\n";
                print {$mfafh} substr($seq{$hapNO}, $_ * $linebase, $linebase)."\n" for ( 0 .. int( (length($seq{$hapNO})-1) / $linebase ) );
            }
        }
    }
    # inform
    stout_and_sterr "[INFO]\toutput mutated sequences of $segName.\n";
}

#--- introduce single mutation to given seq ---
sub add_single_mut{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $type = $parm{type};
    my $alt  = $parm{alt};
    my $pidx = $parm{pidx};
    my $seqSf = $parm{seqSf};

    if($type =~ /^snv$/i){
        return if $alt !~ /^[ACGT]$/i;
        $$seqSf =  substr($$seqSf, 0, $pidx)
                  .$alt
                  .substr($$seqSf, $pidx+1);
    }
    elsif($type =~ /^ins$/i){
        return if $alt !~ /^\+[ACGT]+$/i;
        $alt =~ s/^\+//;
        $$seqSf =  substr($$seqSf, 0, $pidx+1)
                  .$alt
                  .substr($$seqSf, $pidx+1);
    }
    elsif($type =~ /^del$/i){
        return if $alt !~ /^\-[ACGT]+$/i && $alt !~ /^\-\d+$/;
        $alt =~ s/^\-//;
        my $del_len = ($alt =~ /^\d+$/ ? $alt : length($alt));
        $$seqSf =  substr($$seqSf, 0, $pidx)
                  .substr($$seqSf, $pidx+$del_len);
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
