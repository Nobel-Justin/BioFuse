package BioFuse::BioInfo::FASTA::GetMutFA;

use strict;
use warnings;
use Getopt::Long;
use List::Util qw/ min max /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::Interval qw/ Get_Two_Seg_Olen /;
use BioFuse::LoadOn;
use BioFuse::BioInfo::FASTA qw/ read_fasta_file /;

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
$VERSION = "0.01";
$DATE = '2019-09-26';

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
            # option
            [ flankSize => 1000 ],
            # container
            [ mut => {} ],
            # list to abs-path
            [ ToAbsPath_Aref => [ ['mut_fa'],
                                  ['mut_list'],
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
        # option
        "-s:i"  => \$V_Href->{flankSize},
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
            );
}

#--- get fasta file from reference genome based on given mutations ---
sub GetMutFasta{
    # load mutations
    &load_mutations;

    # merge mutations with flank size
    &merge_mutations;

    # modify seq from ref-fasta based on mutations and output
    &get_mut_refseq;
}

#--- load mutation list ---
sub load_mutations{
    open (ML, Try_GZ_Read($V_Href->{mut_list})) || die "fail read mutation list: $!\n";
    # theme
    #chr stP edP type ref alt ref_cn alt_cn gene
    (my $theme_line = lc(<ML>)) =~ s/^#//;
    my @theme_tag = split /\s+/, $theme_line;
    while(<ML>){
        next if(/^\#/);
        my @info = split;
        my %mut = map{ ($theme_tag[$_],$info[$_]) } (0 .. $#theme_tag);
        my $mut_href = { stp => $mut{stp}, edp => $mut{edp},
                         ref => $mut{ref}, ref_cn => $mut{ref_cn},
                         alt => $mut{alt}, alt_cn => $mut{alt_cn},
                         type => $mut{type}, gene => $mut{gene}
                       };
        push @{$V_Href->{mut}->{$mut{chr}}}, $mut_href;
    }
    close ML;
    # inform
    stout_and_sterr "[INFO]\tload mutation list.\n";
}

#--- merge mutations with flank size ---
sub merge_mutations{
    for my $chr (sort keys %{$V_Href->{mut}}){
        @{$V_Href->{mut}->{$chr}} = sort {$a->{stp}<=>$b->{stp}} @{$V_Href->{mut}->{$chr}};
        my $href_1st = shift @{$V_Href->{mut}->{$chr}};
        my @merged = ({ stp => $href_1st->{stp} - $V_Href->{flankSize},
                        edp => $href_1st->{edp} + $V_Href->{flankSize},
                        mut => [$href_1st]
                       });
        for my $mut_href (@{$V_Href->{mut}->{$chr}}){
            my $mfstp = $mut_href->{stp} - $V_Href->{flankSize};
            my $mfedp = $mut_href->{edp} + $V_Href->{flankSize};
            my $ovlen = Get_Two_Seg_Olen($merged[-1]->{stp}, $merged[-1]->{edp}, $mfstp, $mfedp);
            if($ovlen){ # merge!
                # $merged[-1]->{stp} = min($merged[-1]->{stp}, $mfstp); # useless
                $merged[-1]->{edp} = max($merged[-1]->{edp}, $mfedp);
                push @{$merged[-1]->{mut}}, $mut_href;
            }
            else{ # seperate
                push @merged, { stp => $mfstp, edp => $mfedp, mut => [$mut_href] };
            }
        }
        # update
        $V_Href->{mut}->{$chr} = \@merged;
        # inform
        stout_and_sterr "[INFO]\tgroup mutations on $chr.\n";
    }
}

#--- modify seq from ref-fasta based on mutations and output ---
sub get_mut_refseq{
    # read fasta file
    open (my $mfafh, Try_GZ_Write($V_Href->{mut_fa})) || die "fail write output mut_fa file: $!\n";
    read_fasta_file( FaFile => $V_Href->{whole_genome},
                     subrtRef => \&add_mut_on_chr,
                     subrtParmAref => [mfafh => $mfafh, mut_href => $V_Href->{mut}]
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
    my $mut_href = $parm{mut_href};

    return unless exists $mut_href->{$segName};

    my $linebase = 50;
    for my $grp_Hf (@{$mut_href->{$segName}}){
        my $h1_seq = substr($$segSeq_Sref, $grp_Hf->{stp}-1, $grp_Hf->{edp}-$grp_Hf->{stp}+1);
        my $h2_seq = $h1_seq; # first prior to mut
        my $mut_Af = $grp_Hf->{mut};
        # SNV
        my @SNV = grep $_->{type} =~ /^snv$/i, @$mut_Af;
        for my $snv (@SNV){
            my $pidx = $snv->{stp} - $grp_Hf->{stp};
            &add_single_mut(type=>$snv->{type}, alt=>$snv->{alt}, pidx=>$pidx, seqSf=>\$h1_seq) if $snv->{ref_cn} == 0;
            &add_single_mut(type=>$snv->{type}, alt=>$snv->{alt}, pidx=>$pidx, seqSf=>\$h2_seq) if $snv->{alt_cn} != 0;
        }
        # INDEL, reversed position
        my @INDEL = sort {$b->{stp}<=>$a->{stp}} grep $_->{type} =~ /^(ins)|(del)$/i, @$mut_Af;
        for my $indel (@INDEL){
            my $pidx = $indel->{stp} - $grp_Hf->{stp};
            &add_single_mut(type=>$indel->{type}, alt=>$indel->{alt}, pidx=>$pidx, seqSf=>\$h1_seq) if $indel->{ref_cn} == 0;
            &add_single_mut(type=>$indel->{type}, alt=>$indel->{alt}, pidx=>$pidx, seqSf=>\$h2_seq) if $indel->{alt_cn} != 0;
        }
        # output
        print {$mfafh} ">$segName:$grp_Hf->{stp}-$grp_Hf->{edp}:hap1\n";
        print {$mfafh} substr($h1_seq, $_ * $linebase, $linebase)."\n" for ( 0 .. int( (length($h1_seq)-1) / $linebase ) );
        print {$mfafh} ">$segName:$grp_Hf->{stp}-$grp_Hf->{edp}:hap2\n";
        print {$mfafh} substr($h2_seq, $_ * $linebase, $linebase)."\n" for ( 0 .. int( (length($h2_seq)-1) / $linebase ) );
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
        $$seqSf =  substr($$seqSf, 0, $pidx)
                  .$alt
                  .substr($$seqSf, $pidx+1);
    }
    elsif($type =~ /^ins$/i){
        $$seqSf =  substr($$seqSf, 0, $pidx+1)
                  .$alt
                  .substr($$seqSf, $pidx+1);
    }
    elsif($type =~ /^del$/i){
        $$seqSf =  substr($$seqSf, 0, $pidx)
                  .substr($$seqSf, $pidx+length($alt));
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
