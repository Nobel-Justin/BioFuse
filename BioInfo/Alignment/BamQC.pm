package BioFuse::BioInfo::Alignment::BamQC;

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ file_exist /;
use BioFuse::Util::GZfile qw/ Try_GZ_Write /;
use BioFuse::LoadOn;
use BioFuse::BioInfo::Objects::SeqData::Bam_OB;
use BioFuse::BioInfo::Objects::SeqData::Reads_OB;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              BamQC
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Alignment::BamQC';
#----- version --------
$VERSION = "0.01";
$DATE = '2019-11-09';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        return_HELP_INFO
                        Load_moduleVar_to_pubVarPool
                        Get_Cmd_Options
                        para_alert
                        BamQC
                        prepare
                        calcIndex
                        calcIndexFromBam
                        report
                      /;

#--- return HELP_INFO ---
sub return_HELP_INFO{
 return "
     Usage:   perl $V_Href->{MainName} bam_qc <[Options]>

     Options:
         -b  [s]  bam file. <required>
         -o  [s]  output folder. <required>
         -i  [s]  sample ID. <required>
         -s  [s]  samtools path. <required>
         -t  [s]  target region BED file, can be used multiple times. <required>
                  1) Format: -t 'key1:xxx.bed' -t 'key2:yyy.bed'
                  2) Please use the non-N region BED file for WGS, see func 'get_nonN'.
         -d  [s]  accumulated depth to stat for region listed in BED files. ['4,10,20']
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
            [ bam_path => undef ],
            [ output_dir => undef ],
            [ sampleID => undef ],
            [ target_bed => [] ],
            [ samtools => undef ],
            # option
            [ accuDepth => '4,10,20' ],
            # container
            [ bam => undef ],
            [ bed => {} ],
            [ QC => {} ], # result
            [ accuDepAf => [] ],
            # list to abs-path
            [ ToAbsPath_Aref => [ ['bam_path'],
                                  ['output_dir'],
                                  ['samtools']  ] ]
        );
}

#--- get options from command line ---
sub Get_Cmd_Options{
    # get options
    GetOptions(
        # input/output
        "-b:s"  => \$V_Href->{bam_path},
        "-o:s"  => \$V_Href->{output_dir},
        "-i:s"  => \$V_Href->{sampleID},
        # "-n:s"  => \$V_Href->{nonN_bed},
        "-t:s"  => \@{$V_Href->{target_bed}},
        # tools
        "-s:s"  => \$V_Href->{samtools},
        # option
        "-d:s"  => \$V_Href->{accuDepth},
        # help
        "-h|help"   => \$V_Href->{HELP},
        # for debug
        "-debug"    => \$V_Href->{in_debug} # hidden option
    );
}

#--- test para and alert ---
sub para_alert{
    return  (   $V_Href->{HELP}
             || !file_exist(filePath=>$V_Href->{bam_path})
             || !file_exist(filePath=>$V_Href->{output_dir})
             || !defined $V_Href->{sampleID}
             || !defined $V_Href->{samtools}
             || scalar @{$V_Href->{target_bed}} == 0
            );
}

#--- summarize QC result of given bam file ---
sub BamQC{
    # prepare
    &prepare;

    # QC Index
    &calcIndex;

    # report
    &report;
}

#--- prepare work ---
sub prepare{
    # bam
    $V_Href->{bam} = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath => $V_Href->{bam_path}, tag => 'toQC');
    $V_Href->{bam}->addTool(samtools => $V_Href->{samtools});
    # bed files
    for my $bedOpt (@{$V_Href->{target_bed}}){
        my ($key, $bed) = (split /,/, $bedOpt);
        unless(defined $key && defined $bed && -e $bed){
            warn_and_exit "<ERROR>\tbad -t option: $bedOpt\n";
        }
        if(exists $V_Href->{bed}->{$key}){
            warn_and_exit "<ERROR>\t$key bed already exists.\n";
        }
        $V_Href->{bed}->{$key} = $bed;
    }
    # accumulate depth
    $V_Href->{accuDepAf} = [split /,/, $V_Href->{accuDepth}];
    # inform
    stout_and_sterr "[INFO]\tprepare well.\n";
}

#--- QC Index on total level ---
sub calcIndex{
    # total level. common stats
    &calcIndexFromBam(key => 'total', viewOpt => '-F 0x100'); # 0x100 secondary alignment
    # targets
    &calcIndexFromBam(key => $_, viewOpt => "-F 0x100 -L $V_Href->{bed}->{$_}") for keys %{$V_Href->{bed}};
}

#--- calculate QC Indexes from bam ---
sub calcIndexFromBam{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $key = $parm{key} || 'total';
    my $viewOpt = $parm{viewOpt} || '';
    # load bam
    my $fh = $V_Href->{bam}->start_read(viewOpt => $viewOpt);
    while(<$fh>){
        my $reads = BioFuse::BioInfo::Objects::SeqData::Reads_OB->new(ReadsLineText => $_, _rc_optfd => 1, _rc_rdstr => 1);
        # clean reads, except of supplementary alignments
        unless($reads->is_suppmap){
            $V_Href->{QC}->{$key}->{CleanReads}++;
            $V_Href->{QC}->{$key}->{CleanBases} += $reads->rlen;
        }
        # aligned bases
        next if $reads->is_unmap;
        my $mapLen = $reads->mReadLen;
        $V_Href->{QC}->{$key}->{AlignBases} += $mapLen;
        $V_Href->{QC}->{$key}->{UniqBases}  += $mapLen unless $reads->is_mltmap;
        # supplementary alignment only contributes to aligned bases
        next if $reads->is_suppmap;
        # intial alignment
        $V_Href->{QC}->{$key}->{AlignReads}++;
        $V_Href->{QC}->{$key}->{DupReads}++ if $reads->is_dup;
        unless($reads->is_mltmap){
            $V_Href->{QC}->{$key}->{UniqReads}++;
            $V_Href->{QC}->{$key}->{MismatchBases} += $reads->mmCount;
            $V_Href->{QC}->{$key}->{EditDistance} += $reads->editDist;
        }
    }
    close $fh;
    # ratios
    $V_Href->{QC}->{$key}->{AlignRate} = sprintf "%.2f%%", 100 * $V_Href->{QC}->{$key}->{AlignReads} / $V_Href->{QC}->{total}->{CleanReads};
    $V_Href->{QC}->{$key}->{UniqRate}  = sprintf "%.2f%%", 100 * $V_Href->{QC}->{$key}->{UniqReads}  / $V_Href->{QC}->{$key}->{AlignReads};
    $V_Href->{QC}->{$key}->{DupRate}   = sprintf "%.2f%%", 100 * $V_Href->{QC}->{$key}->{DupReads}   / $V_Href->{QC}->{$key}->{AlignReads};
    $V_Href->{QC}->{$key}->{MismatchRate} = sprintf "%.2f%%", 100 * $V_Href->{QC}->{$key}->{MismatchBases} / $V_Href->{QC}->{$key}->{UniqBases};
    # depth and coverage
    if($key ne 'total'){
        my $statHf = $V_Href->{bam}->get_regionCovStat(bed => $V_Href->{bed}->{$key}, accuDepAf => $V_Href->{accuDepAf});
        $V_Href->{QC}->{$key}->{AverageDepth} = $statHf->{meanDepth};
        $V_Href->{QC}->{$key}->{Coverage}->{1}  = sprintf "%.2f%%", 100 * $statHf->{coverage};
        $V_Href->{QC}->{$key}->{Coverage}->{$_} = sprintf "%.2f%%", 100 * $statHf->{accuDepPt}->{$_} for keys %{$statHf->{accuDepPt}};
    }
    # inform
    stout_and_sterr "[INFO]\tcalculate QC index for $key.\n";
}

#--- write QC report ---
sub report{
    my @total_idx = qw/ CleanReads
                        CleanBases /;
    my @share_idx = qw/ AlignReads
                        AlignRate
                        AlignBases
                        UniqReads
                        UniqRate
                        UniqBases
                        DupReads
                        DupRate   /;
    my @target_idx = qw/ MismatchBases
                         MismatchRate
                         EditDistance
                         AverageDepth
                         Coverage     /;

    my $report = File::Spec->catfile($V_Href->{output_dir}, "$V_Href->{sampleID}.BamQC.report");
    open (REPORT, Try_GZ_Write($report)) || die "fail write output report file: $!\n";
    print REPORT '##'.`date`;
    print REPORT "##AlignRate(xxx)=AlignReads(xxx)/CleanReads(total)\n";
    print REPORT "##UniqRate(xxx)=UniqReads(xxx)/AlignReads(xxx)\n";
    print REPORT "##DupRate(xxx)=DupReads(xxx)/AlignReads(xxx)\n";
    print REPORT "##MismatchRate(xxx)=MismatchBases(xxx)/UniqBases(xxx)\n";
    print REPORT "#index\tvalue\n";
    my @key = ('total', (grep $_ ne 'total', sort keys %{$V_Href->{QC}}));
    for my $key (@key){
        my @index = @share_idx;
        unshift @index, @total_idx if $key eq 'total';
        push @index, @target_idx if $key ne 'total';
        my $keyQC_hf = $V_Href->{QC}->{$key};
        for my $idx (@index){
            if($idx ne 'Coverage'){
                print REPORT "$idx($key)\t$keyQC_hf->{$idx}\n";
            }
            else{
                print REPORT "$idx(>=${_}X,$key)\t$keyQC_hf->{$idx}->{$_}\n" for sort {$a<=>$b} keys %{$keyQC_hf->{$idx}};
            }
        }
    }
    close REPORT;
    # inform
    stout_and_sterr "[INFO]\twrite QC report.\n";
}

#--- 
1; ## tell the perl script the successful access of this module.
