package BioFuse::BioInfo::Objects::SeqData::Bam_OB;

use strict;
use warnings;
use Cwd qw/ abs_path /;
use Data::Dumper;
use List::Util qw/ first sum min max uniq /;
use BioFuse::Util::Sys qw/ trible_run_for_success file_exist reset_folder remove_folder /;
use BioFuse::Util::Log qw/ cluck_and_exit stout_and_sterr /;
use BioFuse::Util::GZfile qw/ Try_GZ_Write /;
use BioFuse::BioInfo::Alignment::ReAlign qw/ GATK_ReAlign /;
use BioFuse::BioInfo::Alignment::SAMtools qw/ SAMtool_calmd SAMtools_mpileup /;
use BioFuse::BioInfo::Objects::SeqData::Reads_OB;
use BioFuse::BioInfo::Objects::SeqData::PairEnd_OB;
use BioFuse::BioInfo::Objects::SeqData::HicReads_OB;
use BioFuse::BioInfo::Objects::SeqData::HicPairEnd_OB;
use BioFuse::BioInfo::Objects::SeqData::ReadsGroup_OB;
# use BioFuse::BioInfo::Objects::Allele::RefPos_OB

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Objects::SeqData::Bam_OB';
#----- version --------
$VERSION = "0.22";
$DATE = '2022-01-09';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        verify_bam
                        verify_index
                        delete_file
                        addTool
                        filepath
                        tag
                        tissue
                        rgOB_Hf
                        header_Af
                        reheader
                        isNsort
                        toNsort
                        isCsort
                        toCsort
                        toIndex
                        toFixmate
                        toGATKrealn
                        toCalMD
                        toMpileup
                        merge
                        merge_and_sort
                        toMarkdup
                        start_read
                        start_write
                        stop_write
                        write
                        load_reads_for_ReadsGroup
                        rg_count_need_reads_ForIns
                        pick_rgOB
                        add_rgOB
                        reads_count
                        get_regionCovStat
                        get_region_depth
                        delete_regionDepthFile
                        get_region_alt_vcf_gz
                        get_bam_stats
                        get_pos_marker_stat
                        get_allele_marker_stat
                        smartBam_PEread
                        get_cigarIDSC
                        reads2fasta
                     /;

#--- structure of object
# bam -> filepath = $filepath
# bam -> tag = $tag
# bam -> tissue = $tissue
# bam -> tools = {samtools => $samtools, bcftools => $bcftools, gatk => $gatk_jar, java => $java}
# bam -> write_fh = $write_fh
# bam -> rgOB = { rgID->$rgOB }, check BioFuse::BioInfo::Objects::SeqData::ReadsGroup_OB
# bam -> regionDepthFile = $regionDepthFile

#--- construction of object ---
sub new{
    my $type = shift;
    my %parm = @_;

    my $bam = {};
    $bam->{filepath} = $parm{filepath} || undef;
    $bam->{tag} = $parm{tag} || undef;
    $bam->{tissue} = $parm{tissue} || undef;
    $bam->{rgOB} = {};

    # $bam->{filepath} = abs_path $bam->{filepath} if defined $bam->{filepath};

    bless($bam);
    return $bam;
}

#--- verify bam file ---
sub verify_bam{
    my $bam = shift;
    $bam->{filepath} = abs_path $bam->{filepath};
    if( ! file_exist(filePath => $bam->{filepath}) ){
        cluck_and_exit "<ERROR>\tCannot find bam:\n".Dumper($bam);
    }
}

#--- verify index (bai file) ---
sub verify_index{
    my $bam = shift;
     my $bai_filepath1 = $bam->{filepath} . '.bai';
    (my $bai_filepath2 = $bam->{filepath}) =~ s/bam$/bai/;
    if(    ! file_exist(filePath => $bai_filepath1)
        && ! file_exist(filePath => $bai_filepath2)
    ){
        cluck_and_exit "<ERROR>\tCannot find index bai file of bam:\n"
                            ."\t$bam->{filepath}\n";
    }
}

#--- delete files ---
sub delete_file{
    my $bam = shift;
    `rm -f $bam->{filepath}` if -e $bam->{filepath};
     my $bai_filepath1 = $bam->{filepath} . '.bai';
    (my $bai_filepath2 = $bam->{filepath}) =~ s/bam$/bai/;
    `rm -f $bai_filepath1` if -e $bai_filepath1;
    `rm -f $bai_filepath2` if -e $bai_filepath2;
}

#--- add tool ---
sub addTool{
    my $bam = shift;
    my %parm = @_;
    $bam->{tools}->{lc($_)} = $parm{$_} for keys %parm;
}

#--- return file path ---
sub filepath{
    my $bam = shift;
    return $bam->{filepath};
}

#--- return bam tag ---
sub tag{
    my $bam = shift;
    return $bam->{tag};
}

#--- return tissue ---
sub tissue{
    my $bam = shift;
    return $bam->{tissue};
}

#--- return Hash-ref bam's rgOB ---
sub rgOB_Hf{
    my $bam = shift;
    return $bam->{rgOB};
}

#--- return refer of SAM header content array ---
sub header_Af{
    my $bam = shift;
    my %parm = @_;
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};

    my @header;
    open (HEAD, "$samtools view -H $bam->{filepath} |") || die "fail reading SAM header: $!\n";
    while(<HEAD>){ push @header, $_; }
    close HEAD;

    return \@header;
}

#--- output reheader bam ---
sub reheader{
    my $bam = shift;
    my %parm = @_;
    my $reheaderBam = $parm{reheaderBam};
    my $headerSAM = $parm{headerSAM}; # header SAM file
    my $headerAf = $parm{headerAf};   # or see bam_OB->header_Af func
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};

    if(defined $headerAf){
        $headerSAM = $reheaderBam->filepath.'.reheader.temp.header.sam';
        open (HS,Try_GZ_Write($headerSAM)) || die "fail write $headerSAM: $!\n";
        print HS for @$headerAf;
        close HS;
    }

    my $cmd = "$samtools reheader -P $headerSAM $bam->{filepath} > $reheaderBam->{filepath}";
    trible_run_for_success($cmd, 'reheader', {esdo_Nvb=>1});
    `rm $headerSAM` if defined $headerAf;
}

#--- test whether is N-sort ---
sub isNsort{
    my $bam = shift;
    my $HDline = first {/^\@HD/} @{$bam->header_Af};
    return (defined $HDline && $HDline =~ /SO:queryname/);
}

#--- output N-sort bam ---
sub toNsort{
    my $bam = shift;
    my %parm = @_;
    my $nSortBam = $parm{nSortBam};
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};

    # already Nsort?
    if($bam->isNsort){
        `cp $bam->{filepath} $nSortBam->{filepath}`;
        stout_and_sterr "[INFO]\torigBam is already Nsort bam, Just Copy the bam file.\n"
                             ."\torigBam: $bam->{filepath}\n"
                             ."\tsortBam: $nSortBam->{filepath}\n";
    }
    # do Nsort
    my $sortTmpFolder = $nSortBam->filepath . '-nSortTmpFolder';
    reset_folder(folder => $sortTmpFolder);
    my $cmd = "$samtools sort --threads 4 -n -m 1G -T $sortTmpFolder/nSort -o $nSortBam->{filepath} $bam->{filepath}";
    trible_run_for_success($cmd, 'nSortBam', {esdo_Nvb=>1});
    remove_folder(folder => $sortTmpFolder);
    stout_and_sterr "[INFO]\ttoNsort bam ok.\n"
                         ."\torigBam: $bam->{filepath}\n"
                         ."\tsortBam: $nSortBam->{filepath}\n";
}

#--- test whether is C-sort ---
sub isCsort{
    my $bam = shift;
    my $HDline = first {/^\@HD/} @{$bam->header_Af};
    return (defined $HDline && $HDline =~ /SO:coordinate/);
}

#--- output C-sort bam ---
sub toCsort{
    my $bam = shift;
    my %parm = @_;
    my $cSortBam = $parm{cSortBam};
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};

    # already Csort?
    if($bam->isCsort){
        `cp $bam->{filepath} $cSortBam->{filepath}`;
        stout_and_sterr "[INFO]\torigBam is already Csort bam, Just Copy the bam file.\n"
                             ."\torigBam: $bam->{filepath}\n"
                             ."\tsortBam: $cSortBam->{filepath}\n";
    }
    # do Csort
    my $sortTmpFolder = $cSortBam->filepath . '-cSortTmpFolder';
    reset_folder(folder => $sortTmpFolder);
    my $cmd = "$samtools sort --threads 4 -m 1G -T $sortTmpFolder/nSort -o $cSortBam->{filepath} $bam->{filepath}";
    trible_run_for_success($cmd, 'cSortBam', {esdo_Nvb=>1});
    remove_folder(folder => $sortTmpFolder);
    stout_and_sterr "[INFO]\ttoCsort bam ok.\n"
                         ."\torigBam: $bam->{filepath}\n"
                         ."\tsortBam: $cSortBam->{filepath}\n";
}

#--- index this bam if is C-sort ---
sub toIndex{
    my $bam = shift;
    my %parm = @_;
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};

    # check
    cluck_and_exit "<ERROR>\tbam is not C-sort.\n".Dumper($bam) unless $bam->isCsort;
    my $cmd = "$samtools index $bam->{filepath}";
    trible_run_for_success($cmd, 'indexBam', {esdo_Nvb=>1});
    stout_and_sterr "[INFO]\tindex bam ok.\n"
                         ."\tbam: $bam->{filepath}\n";
}

#--- output fixmate bam ---
# requires samtools v1.9+
# -m: Add mate score tag
sub toFixmate{
    my $bam = shift;
    my %parm = @_;
    my $fixmateBam = $parm{fixmateBam};
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};

    my $cmd = "$samtools fixmate -m $bam->{filepath} $fixmateBam->{filepath}";
    trible_run_for_success($cmd, 'fixmateBam', {esdo_Nvb=>1});
    stout_and_sterr "[INFO]\ttoFixmate bam ok.\n"
                         ."\torigBam: $bam->{filepath}\n"
                         ."\tsortBam: $fixmateBam->{filepath}\n";
}

#--- output gatk realn bam ---
sub toGATKrealn{
    my $bam = shift;
    my %parm = @_;
    my $relanBam = $parm{relanBam};
    my $ref = $parm{ref};
    my $java = $parm{java} || $bam->{tools}->{java};
    my $gatk = $parm{gatk} || $bam->{tools}->{gatk};
    my $jmem = $parm{jmem} || '1g';

    GATK_ReAlign( ref => $ref,
                  sBamPath => $bam->{filepath},
                  oBamPath => $relanBam->{filepath},
                  gatk => $gatk,
                  java => $java,
                  jmem => $jmem
                );
}

#--- output samtools calmd bam ---
sub toCalMD{
    my $bam = shift;
    my %parm = @_;
    my $ref = $parm{ref};
    my $calmdBam = $parm{calmdBam};
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};

    SAMtool_calmd( ref => $ref,
                   sBamPath => $bam->{filepath},
                   oBamPath => $calmdBam->{filepath},
                   samtools => $samtools
                 );
}

#--- output mpileup vcf ---
sub toMpileup{
    my $bam = shift;
    my %parm = @_;
    my $ref = $parm{ref};
    my $vcf = $parm{vcf}; # BioFuse::BioInfo::Objects::Variant::VCF_OB
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};

    SAMtools_mpileup( ref => $ref,
                      sBamPath => $bam->{filepath},
                      vcfPath => $vcf->filepath,
                      samtools => $samtools
                    );
}

#--- merge several bams ---
## this bam is the merged bam
sub merge{
    my $bam = shift;
    my %parm = @_;
    my $bamAf = $parm{bamAf};
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};

    # more than one bam
    cluck_and_exit "<ERROR>\tonly one bam to merge.\n".Dumper($bamAf) if @$bamAf <= 1;
    # should be sorted by the same index (N or C)
    if($bamAf->[0]->isNsort){
        for my $source_bam (@$bamAf){
            cluck_and_exit "<ERROR>\tnot N-sorted bam.\n".Dumper($source_bam) unless $source_bam->isNsort;
        }
    }
    elsif($bamAf->[0]->isCsort){
        for my $source_bam (@$bamAf){
            cluck_and_exit "<ERROR>\tnot C-sorted bam.\n".Dumper($source_bam) unless $source_bam->isCsort;
        }
    }
    else{
        cluck_and_exit "<ERROR>\tnot sorted bam.\n".Dumper($bamAf->[0]);
    }
    # do merge
    my @bamPath = map {$_->filepath} @$bamAf;
    my $nSortPara = $bamAf->[0]->isNsort ? '-n' : '';
    my $cmd = "$samtools merge $nSortPara -f -c -p --threads 4 $bam->{filepath} @bamPath";
    trible_run_for_success($cmd, 'mergeBam', {esdo_Nvb=>1});
    # inform
    stout_and_sterr "[INFO]\tSamTools merge BAMs ok.\n"
                         ."\tbam: $bam->{filepath}\n";
}

#--- merge several bams and N/C-sort---
## this bam is the final merged-sorted bam
sub merge_and_sort{
    my $bam = shift;
    my %parm = @_;
    my $bamAf = $parm{bamAf};
    my $sortBy = $parm{sortBy}; # 'N' or 'C'
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};

    # temp merged bam
    my $tmp_mergeBam = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath => $bam->filepath.'.tmp_merge.bam', tag => 'tmp_mergeBam');
    $tmp_mergeBam->addTool(samtools => $samtools);
    $tmp_mergeBam->merge(bamAf => $bamAf, samtools => $samtools);
    # sort
    if($sortBy =~ /^C$/i){
        $tmp_mergeBam->toCsort(cSortBam => $bam, samtools => $samtools);
    }
    elsif($sortBy =~ /^N$/i){
        $tmp_mergeBam->toNsort(nSortBam => $bam, samtools => $samtools);
    }
    else{
        cluck_and_exit "<ERROR>\twrong sortBy opt: $sortBy. Should be 'N' or 'C'.\n";
    }
    # sweep
    $tmp_mergeBam->delete_file;
    # inform
    stout_and_sterr "[INFO]\tSamTools merge BAMs and sort ok.\n"
                         ."\tbam: $bam->{filepath}\n";
}

#--- output markdup bam ---
# requires samtools v1.9+
sub toMarkdup{
    my $bam = shift;
    my %parm = @_;
    my $markdupBam = $parm{markdupBam};
    my $max_rlen = $parm{max_rlen};
    my $markDupSupp = $parm{markDupSupp} || 0;
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};

    my $opt =  (defined $max_rlen ? "-l $max_rlen " : '')
              .($markDupSupp      ? '-S '           : '');
    my $cmd = "$samtools markdup $opt $bam->{filepath} $markdupBam->{filepath}";
    trible_run_for_success($cmd, 'markdupBam', {esdo_Nvb=>1});
    # inform
    stout_and_sterr "[INFO]\tSamTools markdup bam ok.\n"
                         ."\torigBam:  $bam->{filepath}\n"
                         ."\tmkdupBam: $markdupBam->{filepath}\n";
}

#--- open file-handle to start reading ---
sub start_read{
    my $bam = shift;
    my %parm = @_;
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};
    my $viewOpt = $parm{viewOpt} || '';
    my $viewReg = $parm{viewReg} || '';

    # check existence
    $bam->verify_bam;
    # open fh
    open (my $readFH,"$samtools view $viewOpt $bam->{filepath} $viewReg |") || die "fail reading: $!\n".Dumper($bam);

    return $readFH;
}

#--- open file-handle to start writing ---
sub start_write{
    my $bam = shift;
    my %parm = @_;
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};

    # check
    if(exists $bam->{write_fh}){
        cluck_and_exit "<ERROR>\tthe write file-handle of bam_OB already exists.\n".Dumper($bam);
    }
    # open fh
    open ($bam->{write_fh}, "| $samtools view -b -o $bam->{filepath}") || die "fail writing: $!\n".Dumper($bam);
}

#--- close file-handle to stop writing ---
sub stop_write{
    my $bam = shift;
    close  $bam->{write_fh} if exists $bam->{write_fh};
    delete $bam->{write_fh};
}

#--- write contect to bam ---
sub write{
    my $bam = shift;
    my %parm = @_;
    my $content = $parm{content};
    print {$bam->{write_fh}} $content;
}

#--- load reads and INS-stat for all ReadsGroup ---
sub load_reads_for_ReadsGroup{
    my $bam = shift;
    my %parm = @_;
    my $rgid2rgOB_Href = $parm{rgid2rgOB_Href};
    my $only_SoftClip  = $parm{only_SoftClip};
    my $pr_DropProb    = $parm{pr_DropProb} || 0;
    my $pr_AimCount    = $parm{pr_AimCount} || 10000;
    my $Tool_Tag       = $parm{Tool_Tag} || 'BF';
    my $samtools       = $parm{samtools} || $bam->{tools}->{samtools};

    # still has reads group to load reads?
    my $rgNeedReads = $bam->rg_count_need_reads_ForIns;
    return if( $rgNeedReads == 0 );

    # FLAG: -f 0x2 (P), -F 0x100(sd) + 0x400(d) + 0x800(sp) = 0xD00
    open (BAM,"$samtools view -f 0x2 -F 0xD00 $bam->{filepath} |") || die"fail read bam: $!\n";
    while(<BAM>){
        # load reads object
        my $reads_OB = BioFuse::BioInfo::Objects::SeqData::Reads_OB->new( ReadsLineText => $_, _rc_optfd => 0, _rc_rdstr => 0 );
        # check RG_ID
        my $rg_OB = $reads_OB->find_rgOB(rgid2rgOB_Href => $rgid2rgOB_Href);
        ## enough reads for ins evalue
        next if( $rg_OB->{EnoughReadsBool} );
        # extract read-id FS prefix to rg_OB, and prepare for update
        $reads_OB->update_rid_RGprefix(Tool_Tag => $Tool_Tag) unless defined $rg_OB->{rID_prefix} ;
        # count in this paired reads for RG ins evalue
        $rg_OB->load_reads_for_ins_evalue( reads_OB => $reads_OB, only_SoftClip => $only_SoftClip, pr_AimCount => $pr_AimCount, pr_DropProb => $pr_DropProb );
        # stop when all RG(s) are full
        $rgNeedReads-- if( $rg_OB->{EnoughReadsBool} );
        last if( $rgNeedReads == 0 );
    }
    close BAM;

    # inform
    stout_and_sterr "[INFO]\tload paired-end info for reads group (only_SoftClip:$only_SoftClip) from bam file.\n"
                         ."\t$bam->{filepath}\n";
}

#--- whether all RG(s) of this bam get enough reads for ins evaluation ---
sub rg_count_need_reads_ForIns{

    my $bam = shift;

    my $count = 0;
    for my $rg_OB ( values %{$bam->{rgOB}} ){
        if( $rg_OB->{EnoughReadsBool} == 0 ){
            $count++;
        }
    }

    return $count;
}

#--- pick reads group and create related objects ---
sub pick_rgOB{
    my $bam = shift;
    my %parm = @_;
    my $rgid2rgOB_Href = $parm{rgid2rgOB_Href};
    my $samtools       = $parm{samtools} || $bam->{tools}->{samtools};
    # my $Allow_Lack_RG  = $parm{Allow_Lack_RG};

    # read bam header
    open (BAMHEADER,"$samtools view -H $bam->{filepath} |") || die"fail read bam: $!\n";
    while(<BAMHEADER>){
        if(/^\@RG/){ # the read group (RG)
            my ($RG_ID, $LB_ID); # add LB_ID for supporting-reads de-duplication
            if(/\sID\:(\S+).*\sLB\:(\S+)/){
                ($RG_ID, $LB_ID) = ($1, $2);
            }
            else{
                cluck_and_exit "<ERROR>\tCannot get RG_ID or LB_ID from header.\n".
                                      "\t$bam->{filepath}\n";
            }
            # check recurrence
            if( exists $rgid2rgOB_Href->{$RG_ID} ){
                cluck_and_exit "<ERROR>\tThe RG_ID $RG_ID of bam file is recurrent.\n".
                                      "\t$bam->{filepath}\n";
            }
            # create reads group object (rg_OB)
            $rgid2rgOB_Href->{$RG_ID} = BioFuse::BioInfo::Objects::SeqData::ReadsGroup_OB->new(bam => $bam, RG_ID => $RG_ID, LB_ID => $LB_ID);
            # link the rg_OB with the bam
            $bam->add_rgOB(rgOB => $rgid2rgOB_Href->{$RG_ID});
        }
    }
    close BAMHEADER;

    # once the bam lacks reads group
    if( scalar(keys %{$bam->{rgOB}}) == 0 ){
        if( 0 ){ # $Allow_Lack_RG
            # future updates 
        }
        else{
            cluck_and_exit "<ERROR>\tbam file lacks reads group info.\n".
                                  "\t$bam->{filepath}\n";
        }
    }

    # inform
    stout_and_sterr "[INFO]\tcreate reads group objects from bam file.\n"
                         ."\t$bam->{filepath}\n";
}

#--- add reads group object(s) to this bam ---
sub add_rgOB{
    my $bam = shift;
    my %parm = @_;
    $bam->{rgOB}->{$parm{rgOB}->RG_ID} = $parm{rgOB};
}

#--- count reads of given region ---
sub reads_count{
    my $bam = shift;
    my %parm = @_;
    my $region = $parm{region};
    my $bed = $parm{bed};
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};
    my $onlyPropMap = $parm{onlyPropMap} || 0; # only properly mapped
    my $skipSuppMap = $parm{skipSuppMap} || 0; # skip supplement/second alignment
    my $skipDupRead = $parm{skipDupRead} || 0; # skip PCR duplicated reads
    my $minMQ = defined $parm{minMQ} ? $parm{minMQ} : 30; # # mapping quality
    my $maxClipRatio = $parm{maxClipRatio} || 0; # Sofe/Hart-clip part ratio

    # check existence
    $bam->verify_bam;

    # read bam
    my $reads_count = 0;
    my $viewOpt  = "-q $minMQ ";
       $viewOpt .= "-f 0x2 " if $onlyPropMap;
       $viewOpt .= "-F 0x900 " if $skipSuppMap; # -F 0x800 + 0x100
       $viewOpt .= "-F 0x400 " if $skipDupRead;
    my $viewReg  = defined $region ? $region : (defined $bed ? "-L $bed" : "");
    my $fh = $bam->start_read(samtools=>$samtools, viewOpt=>$viewOpt, viewReg=>$viewReg);
    while(<$fh>){
        my $reads_OB = BioFuse::BioInfo::Objects::SeqData::Reads_OB->new(ReadsLineText => $_);
        next if $maxClipRatio == 0 && $reads_OB->is_clip;
        $reads_count++ if ($reads_OB->biClipLen / $reads_OB->rlen) <= $maxClipRatio;
    }
    close $fh;

    return $reads_count;
}

#--- get coverage/depth stat of given region ---
sub get_regionCovStat{
    my $bam = shift;
    my %parm = @_;
    my $region = $parm{region};
    my $bed = $parm{bed};
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};
    my $notAllPos = $parm{notAllPos} || 0; # donot use '-a'
    my $minMQ = $parm{minMQ} || 0;
    my $minBQ = $parm{minBQ} || 0;
    my $maxDP = defined $parm{maxDP} ? $parm{maxDP} : 8000;
    my $accuDepAf = $parm{accuDepAf} || [10,30,50,100];
    my $origLen = $parm{origLen};
    my $circExtl = $parm{circExtl} || 0; # for such as circular virus
    my $availLen = $parm{availLen} || 0; # for such as non-N region of virus
    my $minDepth = $parm{minDepth} || 1; # pos less than such depth will be skipped
    my $onlyDepStat = $parm{onlyDepStat} || 0; # return hash of each pos depth

    # check existence
    $bam->verify_bam;

    # return hash
    my $statHf = { coverage => 0,
                   meanDepth => 0,
                   accuDepPt => {map {($_,0)} @$accuDepAf}
                 };

    # samtools depth options
    my $depthOpt = "-d $maxDP -q $minBQ -Q $minMQ";
    unless($notAllPos){
        $depthOpt .= ' -a';
    }
    if($region){
        $depthOpt .= " -r $region";
    }
    elsif($bed){
        $depthOpt .= " -b $bed";
    }

    # pos depth
    my %depthCount;
    my $basicPosCount = 0;
    my %pos2depth; # for circExtl or just_return_depth_stats
    open (DEPTH, "$samtools depth $depthOpt $bam->{filepath} |") || die "fail samtools depth: $!\n";
    while (<DEPTH>){
        my ($pos, $depth) = (split)[1,2];
        if(    $onlyDepStat
            || (   $circExtl
                && ($pos > $origLen || $pos <= $circExtl)
               )
        ){
            $pos2depth{$pos} = $depth;
        }
        else{
            $depthCount{$depth}++;
        }
        # position count
        if(    !$circExtl
            || $pos <= $origLen
        ){
            $basicPosCount++;
        }
    }
    close DEPTH;

    # sometimes, no depth info obtained
    ## e.g., all alignments in the bam are duplicated
    return $statHf unless (keys %pos2depth || keys %depthCount);

    # do for circular extend part
    if($circExtl){
        for my $i (1 .. $circExtl){
            my $iExt = $origLen + $i;
            # skip if both pos are not in pos list
            next if !exists $pos2depth{$i} || !exists $pos2depth{$iExt};
            $pos2depth{$i} += ($pos2depth{$iExt} || 0);
            delete $pos2depth{$iExt};
        }
        # into depthCount
        $depthCount{$pos2depth{$_}}++ for keys %pos2depth;
    }

    # if just want depth stats
    return \%pos2depth if $onlyDepStat;
    undef %pos2depth;

    # filter pos and 
    my @accuDepth = sort {$b<=>$a} keys %{$statHf->{accuDepPt}};
    for my $depth (keys %depthCount){
        # minDepth filter, donot take outlier into account
        if($depth < $minDepth){
            delete $depthCount{$depth};
            next;
        }
        # accumulate depth
        my $i = first {$depth >= $accuDepth[$_]} (0 .. $#accuDepth);
        if(defined $i){
            $statHf->{accuDepPt}->{$accuDepth[$_]} += $depthCount{$depth} for ($i .. $#accuDepth);
        }
    }
    # reset basicPosCount
    $basicPosCount = $availLen || $basicPosCount || 1;
    # get accu-depth-perct
    for my $depth (@accuDepth){
        $statHf->{accuDepPt}->{$depth} = sprintf "%0.4f", min($statHf->{accuDepPt}->{$depth}, $basicPosCount) / $basicPosCount;
    }
    # coverage
    $statHf->{coverage}  = sprintf "%0.4f", (sum(values %depthCount) || 0) / $basicPosCount;
    # mean depth
    $statHf->{meanDepth} = sprintf "%0.2f", (sum(map {$_*$depthCount{$_}} keys %depthCount) || 0) / $basicPosCount;
    # ordinary return
    return $statHf;
}

#--- get depth of certain region ---
sub get_region_depth{
    my $bam = shift;
    my %parm = @_;
    my $region_refseg = $parm{refseg};
    my $region_st_pos = $parm{st_pos};
    my $region_ed_pos = $parm{ed_pos};
    my $min_baseQ = exists($parm{min_baseQ}) ? $parm{min_baseQ} : 5;
    my $min_mapQ  = exists($parm{min_mapQ})  ? $parm{min_mapQ}  : 10;
    my $no_SoftClip = $parm{no_SoftClip} || 0;
    my $outFilePrefix = $parm{out_prefix};
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};

    unless( $region_refseg && $region_st_pos && $region_ed_pos && $outFilePrefix){
        cluck_and_exit "bam->get_region_depth( refseg=>?, st_pos=>?, ed_pos=>?, out_prefix=>?)\n";
    }

    # check existence
    $bam->verify_bam;

    my $bamPathForDepth = $bam->{filepath};
    my $viewRegParm = "$region_refseg:$region_st_pos-$region_ed_pos";

    # filter soft-clip
    if( $no_SoftClip ){
        my $tmp_bam = "$outFilePrefix.noSoftClip.bam";
        my $tmpbam_cmd = "( $samtools view $bam->{filepath} $viewRegParm | awk '\$6 !~ /S/' | $samtools view -S -b -o $tmp_bam - ) && ( $samtools index $tmp_bam )";
        trible_run_for_success($tmpbam_cmd, 'tmpNoSoftClipBam', {esdo_Nvb=>1});
        $bamPathForDepth = $tmp_bam;
    }

    # get depth
    $bam->{regionDepthFile} = "$outFilePrefix.depth.gz";
    my $depth_cmd = "( $samtools depth -a -r $viewRegParm -q $min_baseQ -Q $min_mapQ $bamPathForDepth | gzip -c > $bam->{regionDepthFile} )";
    trible_run_for_success($depth_cmd, 'bamDepth', {esdo_Nvb=>1});

    # sweep
    if( $bamPathForDepth ne $bam->{filepath} ){
        `rm -rf $bamPathForDepth`;
    }
}

#--- sweep region depth stat file ---
sub delete_regionDepthFile{
    my $bam = shift;
    `rm -rf $bam->{regionDepthFile}` if(-e $bam->{regionDepthFile});
}

#--- get mutation vcf.gz in certain region ---
sub get_region_alt_vcf_gz{
    my $bam = shift;
    my %parm = @_;
    my $vcfgz = $parm{vcfgz};
    my $cmd_name = $parm{cmd_name} || 'get_alt_vcf_gz';
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};
    my $bcftools = $parm{bcftools} || $bam->{tools}->{bcftools};
    my $pos_list = $parm{pos_list};
    my $idx_ref = $parm{idx_ref};
    my $ploidy_set = $parm{ploidy_set} || 'GRCh37';
    my $no_indel = $parm{no_indel} || 0;
    my $snpGap = $parm{snpGap} || 5;
    my $IndelGap = $parm{IndelGap} || 5;
    my $min_altQual = $parm{min_altQual} || 20;
    my $min_mapQual = $parm{min_mapQual} || 30;
    my $min_posDepth = $parm{min_posDepth} || 10;
    my $min_alleleDepth = $parm{min_alleleDepth} || 5;
    my $min_strdDepth = $parm{min_strdDepth} || 2;
    my $only_hetAlt = $parm{only_hetAlt} || 0;

    # check existence
    $bam->verify_bam;

    # commands
    my $hetSNP_cmd = "($samtools mpileup -t 'DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR'".
                     (   defined $pos_list
                       ? " -l $pos_list "
                       : ""
                     ).
                     "-uvf $idx_ref $bam->{filepath} | ".
                     "$bcftools call --ploidy $ploidy_set -cv".
                     (   $no_indel
                       ? " -V indels | "
                       : " | "
                     ).
                     "$bcftools filter -g $snpGap -G $IndelGap ".
                     "-e 'QUAL<$min_altQual ".
                     "|| DP<$min_posDepth ".
                     (   $only_hetAlt
                       ? "|| AD[0]<$min_alleleDepth ".
                         "|| DP4[0]<$min_strdDepth ".
                         "|| DP4[1]<$min_strdDepth "
                       : ""
                     ).
                     "|| AD[1]<$min_alleleDepth ".
                     "|| DP4[2]<$min_strdDepth ".
                     "|| DP4[3]<$min_strdDepth ".
                     "|| MQ<$min_mapQual' ".
                     "-O z -o $vcfgz) && ".
                     "($bcftools index $vcfgz)";
    # run
    trible_run_for_success( $hetSNP_cmd, $cmd_name );
}

#--- get bam stats from samtools ---
sub get_bam_stats{
    my $bam = shift;
    my %parm = @_;
    my $statsFile = $parm{statsFile};
    my $opt = $parm{opt} || undef;
    my $bed = $parm{bed} || undef;
    my $cmd_name = $parm{cmd_name} || 'get_bam_stats';
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};
    # check existence
    $bam->verify_bam;
    # commands
    my $bamStats_cmd =  "$samtools stats "
                       .(defined $opt ? "$opt " : '')
                       .(defined $bed ? "-t $bed " : '')
                       ."$bam->{filepath} > $statsFile";
    # run
    trible_run_for_success( $bamStats_cmd, $cmd_name );
}

#--- do stat on marker at given position ---
## maker: barcode (10x) or PE-id (Hi-C, MinIon and PacBio)
sub get_pos_marker_stat{

    my $bam = shift;
    my %parm = @_;
    my $chr = $parm{chr};
    my $pos = $parm{pos};
    my $marker = $parm{marker}; # '10x', 'Hi-C', 'MinIon', 'PacBio'
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};

    my $PmakerStat_Href = {};
    # FLAG: -F 0x100(sd) + 0x400(d) + 0x800(sp) = 0xD00
    open (BAM,"$samtools view -F 0xD00 $bam->{filepath} $chr:$pos-$pos |") || die"fail read bam: $!\n";
    while(<BAM>){
        # load reads object
        my $reads_OB = BioFuse::BioInfo::Objects::SeqData::Reads_OB->new( ReadsLineText => $_, _rc_optfd => 1, _rc_rdstr => 1 );
        next if( $reads_OB->is_unmap );
        # get marker of pos
        my $marker_str = ( $marker =~ /10x/i ? $reads_OB->barc_10x : $reads_OB->pid );
        $PmakerStat_Href->{$marker_str} ++ if( $marker_str ne 'NA' );
    }
    close BAM;

    return $PmakerStat_Href;
}

#--- do stat on marker of alleles at given position ---
## maker: barcode (10x) or PE-id (Hi-C, MinIon and PacBio)
sub get_allele_marker_stat{

    my $bam = shift;
    my %parm = @_;
    my $chr = $parm{chr};
    my $pos = $parm{pos};
    my $marker = $parm{marker}; # '10x', 'Hi-C', 'MinIon', 'PacBio'
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};

    my $AmakerStat_Href = {};
    # FLAG: -F 0x100(sd) + 0x400(d) + 0x800(sp) = 0xD00
    open (BAM,"$samtools view -F 0xD00 $bam->{filepath} $chr:$pos-$pos |") || die"fail read bam: $!\n";
    while(<BAM>){
        # load reads object
        my $reads_OB = BioFuse::BioInfo::Objects::SeqData::Reads_OB->new( ReadsLineText => $_, _rc_optfd => 1, _rc_rdstr => 1 );
        next if( $reads_OB->is_unmap );
        # get marker of pos allele
        my $alleleInfo_Href = $reads_OB->get_pos_allele( chr => $chr, pos => $pos );
        # record with marker
        my $alleleType = $alleleInfo_Href->{type};
        # if( $allele =~ /^[ACGT]$/ ){ # former version, only for refORsnv
        if( $alleleType ne 'miss' ){
            my $allele;
            if( $alleleType eq 'refORsnv' ){
                $allele = $alleleInfo_Href->{allele};
            }
            elsif( $alleleType eq 'ins' ){
                $allele = 'I,'.$alleleInfo_Href->{allele};
            }
            elsif( $alleleType eq 'del' ){
                $allele = 'D,'.$alleleInfo_Href->{delSize}.','.$alleleInfo_Href->{offset};
            }
            my $marker_str = ( $marker =~ /10x/i ? $reads_OB->barc_10x : $reads_OB->pid );
            $AmakerStat_Href->{$allele}->{$marker_str} ++ if( $marker_str ne 'NA' );
        }
    }
    close BAM;

    return $AmakerStat_Href;
}

#--- read smartPE BAM and do something on each pe_OB/pe_OB_pool ---
sub smartBam_PEread{
    my $bam = shift;
    my %parm = @_;
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};
    my $mark = $parm{mark} || $bam->tag || '';
    my $viewOpt = $parm{viewOpt} || '';
    my $readsType = $parm{readsType} || 'norm'; # norm/HiC
    my $subrtRef = $parm{subrtRef};
    my $subrtParmAref = $parm{subrtParmAref} || [];
    my $peC_ReportUnit = $parm{peC_ReportUnit} || 1E6;
    my $peOB_AbufferSize = $parm{peOB_AbufferSize} || 2E4; # 1E5
    my $deal_peOB_pool = $parm{deal_peOB_pool} || 0;
    my $quiet = $parm{quiet} || 0; # not inform PE Count
    my $simpleLoad = $parm{simpleLoad} || 0; # just record reads' basic information

    # objects modules
    my $rdObjSource = {norm => 'Reads_OB',   HiC => 'HicReads_OB'};
    my $rdObjModule = 'BioFuse::BioInfo::Objects::SeqData::' . $rdObjSource->{$readsType};
    my $peObjSource = {norm => 'PairEnd_OB', HiC => 'HicPairEnd_OB'};
    my $peObjModule = 'BioFuse::BioInfo::Objects::SeqData::' . $peObjSource->{$readsType};

    # read BAM
    my $pe_Count = 0;
    my $pe_CountInform = $peC_ReportUnit;
    my $last_pid = '';
    my $last_peOB = undef;
    my @peOB_pool = ();
    my $subrt_signal = 1;
    my $fh = $bam->start_read(samtools => $samtools, viewOpt => $viewOpt);
    while(<$fh>){
        my $reads_OB = $rdObjModule->new( ReadsLineText => $_, _rc_optfd => 1, _rc_rdstr => 1, simpleLoad => $simpleLoad );
        if( $reads_OB->pid ne $last_pid ){
            # deal last pe_OB
            $pe_Count++;
            push @peOB_pool, $last_peOB if( defined $last_peOB );
            if( $pe_Count % $peOB_AbufferSize == 0 ){
                if($deal_peOB_pool){
                    $subrt_signal = &{$subrtRef}(pe_OB_poolAf => \@peOB_pool, @$subrtParmAref);
                }
                else{
                    $subrt_signal = &{$subrtRef}(pe_OB => $_, @$subrtParmAref) for @peOB_pool;
                }
                # sweep
                @peOB_pool = ();
                # inform
                if($pe_Count >= $pe_CountInform){
                    stout_and_sterr "[INFO]\t".`date`
                                         ."\tload $pe_CountInform PE-reads from $mark bam.\n" unless $quiet;
                    $pe_CountInform += $peC_ReportUnit;
                }
                # stop?
                last if defined $subrt_signal && $subrt_signal == -1;
            }
            # create new pe_OB
            $last_peOB = $peObjModule->new;
            # update
            $last_pid = $reads_OB->pid;
        }
        # load this reads_OB to pe_OB
        $last_peOB->load_reads_OB(reads_OB => $reads_OB);
    }
    close $fh;
    # deal the last one
    if(    defined $last_peOB
        && (   !defined $subrt_signal
            || $subrt_signal != -1   # not stop initiatively
           )
    ){
        push @peOB_pool, $last_peOB;
        if($deal_peOB_pool){
            &{$subrtRef}(pe_OB_poolAf => \@peOB_pool, @$subrtParmAref);
        }
        else{
            &{$subrtRef}(pe_OB => $_, @$subrtParmAref) for @peOB_pool;
        }
        # sweep
        @peOB_pool = ();
    }
    # inform
    stout_and_sterr "[INFO]\t".`date`
                         ."\ttotally, load $pe_Count PE-reads from $mark bam.\n" unless $quiet;
}

#--- collect the InDel (ID) and SoftClip (SC) info from CIGAR ---
sub get_cigarIDSC{
    my $bam = shift;
    my %parm = @_;
    my $viewOpt = $parm{viewOpt} || '-F 0x404'; # 0x400+0x004
    my $viewReg = $parm{viewReg} || '';
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};
    my $refposHf = $parm{refposHf};
    my $cigarIDSC_Hf = $parm{cigarIDSC_Hf}; # stores InDel SoftClip info
    my $cigarSC_Hf = $parm{cigarSC_Hf}; # only store SoftClip info, currently only 'softclip3'
    my $min_mTD = $parm{min_mTD}; # min reads Tail Dist of mutation
    my $readsType = $parm{readsType} || 'norm'; # norm/HiC

    # objects modules
    my $rdObjSource = {norm => 'Reads_OB',   HiC => 'HicReads_OB'};
    my $rdObjModule = 'BioFuse::BioInfo::Objects::SeqData::' . $rdObjSource->{$readsType};

    # read BAM to seek InDel and SoftClip sites
    my $fh = $bam->start_read(samtools=>$samtools, viewOpt=>$viewOpt, viewReg=>$viewReg);
    while(<$fh>){
        my $rOB = $rdObjModule->new( ReadsLineText => $_, _rc_optfd => 1, _rc_rdstr => 1 );
        # filter, must have InDel or SoftClip
        my $CIGAR = $rOB->cigar;
        next if $CIGAR !~ /[IDS]/ || $CIGAR =~ /[^IDSM\d]/; # no others
        # deal mutations
        my $ReadSeq = $rOB->rseq;
        my $mseg = $rOB->mseg;
        my $mPos = $rOB->mpos;
        my $i = 0;
        while(length($CIGAR)!=0){
            my ($mut_pos, $mut_type, $mut_seq);
            my ($len, $type) = ($CIGAR =~ /^(\d+)(\D)/);
            $CIGAR =~ s/^\d+\D//;
            if($type eq 'M'){
                $i += $len;
                $ReadSeq = substr($ReadSeq, $len);
                next;
            }
            elsif($type eq 'I'){
                $mut_pos = $mPos + $i - 1; # assign insertion to the forefront-pos
                $mut_type = 'ins';
                $mut_seq = uc(substr($ReadSeq, 0, $len));
                $ReadSeq = substr($ReadSeq, $len);
                next if( $i < $min_mTD || length($ReadSeq) < $min_mTD); # located at tail of this read
            }
            elsif($type eq 'D'){
                $mut_pos = $mPos + $i; # assign deletion to the first pos of the deleted part
                $i += $len;
                $mut_type = 'del';
                $mut_seq = 'D' x $len; # assign the seq later
                next if( $i < $min_mTD || length($ReadSeq) < $min_mTD); # located at tail of this read
            }
            elsif($type eq 'S'){
                if( $i == 0 ){ # left prime of read
                    $mut_pos = $mPos;
                    $mut_type = 'softclip5';
                }
                else{ # right prime of read
                    $mut_pos = $mPos + $i - 1;
                    $mut_type = 'softclip3';
                    $i += $len;
                }
                $ReadSeq = substr($ReadSeq, $len);
            }
            # filter this site if it already has mutation
            ## BioFuse::BioInfo::Objects::Allele::RefPos_OB
            next if exists $refposHf->{$mut_pos} && $refposHf->{$mut_pos}->has_mutation;
            # record the mutated pos with supporting
            if(!exists $cigarIDSC_Hf->{$mseg}->{$mut_pos}){
                $cigarIDSC_Hf->{$mseg}->{$mut_pos} = BioFuse::BioInfo::Objects::Allele::RefPos_OB->new(pos=>$mut_pos, refAllele=>$refposHf->{$mut_pos}->refAllele);
            }
            my $cigarIDSC_refpos_OB = $cigarIDSC_Hf->{$mseg}->{$mut_pos};
            $cigarIDSC_refpos_OB->addDepth(add=>1); # pretend to be mut_whole_depth
            # store supporting of soft-clip for later subtraction from ref_allele_depth
            if($type eq 'S'){
                # only primer3 impact depth
                if($mut_type eq 'softclip3'){
                    $cigarSC_Hf->{$mseg}->{$mut_pos}->{$mut_type}->[0] ++;
                    $cigarSC_Hf->{$mseg}->{$mut_pos}->{$mut_type}->[1] += ($rOB->is_fw_map ? 1 : 0); # forward support
                    $cigarSC_Hf->{$mseg}->{$mut_pos}->{$mut_type}->[2] += ($rOB->is_rv_map ? 1 : 0); # reverse complemented support
                }
                next; # do not record softclip as mutation
            }
            # store insertion and deletion info
            my $mut_id = "$mut_type,$mut_seq";
            $cigarIDSC_refpos_OB->addMut(mut_id=>$mut_id, depthAf=>[0,0,0,0]) if !$cigarIDSC_refpos_OB->has_mutation(mut_id=>$mut_id);
            $cigarIDSC_refpos_OB->addMutDepth(mut_id=>$mut_id, add_fw=>1) if $rOB->is_fw_map; # forward support
            $cigarIDSC_refpos_OB->addMutDepth(mut_id=>$mut_id, add_rv=>1) if $rOB->is_rv_map; # reverse complemented support
        }
    }
    close $fh;
}

#--- extract reads to fasta ---
sub reads2fasta{
    my $bam = shift;
    my %parm = @_;
    my $fasta = $parm{fasta};
    my $viewOpt = $parm{viewOpt} || '';
    my $viewReg = $parm{viewReg} || '';
    my $samtools = $parm{samtools} || $bam->{tools}->{samtools};
    my $covPosHf = $parm{covPosHf} || {}; # key: position must cover
    my $covPosAf = $parm{covPosAf} || []; # list:position must cover
    my $readsType = $parm{readsType} || 'norm'; # norm/HiC

    # cover pos if given
    my @covPos = sort {$a<=>$b} uniq (keys %$covPosHf,@$covPosAf);
    my $onlyCP = scalar @covPos; # only extract reads cover these pos
    # check C-sort if needs
    if(    !$bam->isCsort
        && ($onlyCP || length($viewReg))
    ){
        cluck_and_exit "<ERROR>\tcovPos or viewReg provided, but bam is not C-sort.\n".Dumper($bam);
    }

    # objects modules
    my $rdObjSource = {norm => 'Reads_OB',   HiC => 'HicReads_OB'};
    my $rdObjModule = 'BioFuse::BioInfo::Objects::SeqData::' . $rdObjSource->{$readsType};

    # read BAM to find reads
    open (FA,Try_GZ_Write($fasta)) || die "fail write $fasta: $!\n";
    my $fh = $bam->start_read(samtools=>$samtools, viewOpt=>$viewOpt, viewReg=>$viewReg);
    while(<$fh>){
        my $rOB = $rdObjModule->new(ReadsLineText => $_, _rc_rdstr => 1);
        my $mpos = $rOB->mpos;
        my $mlen = $rOB->mRefLen;
        # covPos if given
        if($onlyCP){
            shift @covPos while (@covPos && $covPos[0] < $mpos); # rOB surpasses
            last unless @covPos; # no pos need
            next if $covPos[0] > $mpos+$mlen-1; # rOB falls behind
        }
        # output
        print FA $rOB->printFA."\n";
    }
    close $fh;
    close FA;
}

1; ## tell the perl script the successful access of this module.
