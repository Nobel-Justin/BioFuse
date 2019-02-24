package BioFuse::BioInfo::Objects::Bam_OB;

use strict;
use warnings;
use Cwd qw/ abs_path /;
use Data::Dumper;
use BioFuse::Util::Sys qw/ trible_run_for_success file_exist /;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;
use BioFuse::BioInfo::Objects::Reads_OB;
use BioFuse::BioInfo::Objects::PairEnd_OB;
use BioFuse::BioInfo::Objects::HicReads_OB;
use BioFuse::BioInfo::Objects::HicPairEnd_OB;
use BioFuse::BioInfo::Objects::ReadsGroup_OB;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Objects::Bam_OB';
#----- version --------
$VERSION = "0.10";
$DATE = '2018-12-30';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        verify_bam
                        verify_index
                        get_filepath
                        get_tag
                        get_SAMheader
                        start_read
                        start_write
                        stop_write
                        write
                        load_reads_for_ReadsGroup
                        rg_count_need_reads_ForIns
                        extract_ReadsGroup_OB
                        add_ReadsGroup_OBs
                        get_region_depth
                        delete_regionDepthFile
                        get_region_alt_vcf_gz
                        get_pos_marker_stat
                        get_allele_marker_stat
                        smartBam_PEread
                     /;

#--- structure of object
# bam -> filepath = $filepath
# bam -> tag = $tag
# bam -> tissue = $tissue
# bam -> write_fh = $write_fh
# bam -> rgOB = { rgID->$rgOB }, check BioFuse::BioInfo::Objects::ReadsGroup_OB
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

    $bam->{filepath} = abs_path $bam->{filepath} if defined $bam->{filepath};

    bless($bam);
    return $bam;
}

#--- verify bam file ---
sub verify_bam{
    my $bam = shift;

    $bam->{filepath} = abs_path $bam->{filepath};

    if( ! file_exist(filePath => $bam->{filepath}) ){
        warn_and_exit "<ERROR>\tCannot find bam:\n"
                            ."\t$bam->{filepath}\n";
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
        warn_and_exit "<ERROR>\tCannot find index bai file of bam:\n"
                            ."\t$bam->{filepath}\n";
    }
}

#--- return file path ---
sub get_filepath{
    my $bam = shift;
    return $bam->{filepath};
}

#--- return bam tag ---
sub get_tag{
    my $bam = shift;
    return $bam->{tag};
}

#--- return refer of SAM header content array ---
sub get_SAMheader{
    my $bam = shift;
    my %parm = @_;
    my $samtools = $parm{samtools};

    my @header;
    open (HEAD, "$samtools view -H $bam->{filepath} |") || die "fail reading SAM header: $!\n";
    while(<HEAD>){ push @header, $_; }
    close HEAD;

    return \@header;
}

#--- open file-handle to start reading ---
sub start_read{
    my $bam = shift;
    my %parm = @_;
    my $samtools = $parm{samtools};
    my $viewOpt = $parm{viewOpt} || '';

    # check existence
    unless(file_exist(filePath => $bam->get_filepath)){
        warn_and_exit "<ERROR>\tCannot find bam\n".Dumper($bam);
    }
    # open fh
    open (my $readFH,"$samtools view $viewOpt $bam->{filepath} |") || die "fail reading: $!\n".Dumper($bam);

    return $readFH;
}

#--- open file-handle to start writing ---
sub start_write{
    my $bam = shift;
    my %parm = @_;
    my $samtools = $parm{samtools};

    # check
    if(exists $bam->{write_fh}){
        warn_and_exit "<ERROR>\tthe write file-handle of bam_OB already exists.\n".Dumper($bam);
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
    my $Tool_Tag       = $parm{Tool_Tag} || '__NULL__';
    my $samtools       = $parm{samtools};

    # still has reads group to load reads?
    my $rgNeedReads = $bam->rg_count_need_reads_ForIns;
    return if( $rgNeedReads == 0 );

    # to store the pid
    my $tissue = $bam->{tissue};
    # FLAG: -f 0x2 (P), -F 0x100(sd) + 0x400(d) + 0x800(sp) = 0xD00
    open (BAM,"$samtools view -f 0x2 -F 0xD00 $bam->{filepath} |") || die"fail read bam: $!\n";
    while(<BAM>){
        # load reads object
        my $reads_OB = BioFuse::BioInfo::Objects::Reads_OB->new( ReadsLineText => $_, _rc_optfd => 0, _rc_rdstr => 0 );
        # check RG_ID
        my $rg_OB = $reads_OB->get_available_rgOB( bam => $bam, rgid2rgOB_Href => $rgid2rgOB_Href );
        ## enough reads for ins evalue
        next if( $rg_OB->{EnoughReadsBool} );
        # extract read-id FS prefix to rg_OB, and prepare for update
        $reads_OB->extract_FS_readid_prefix( rg_OB => $rg_OB, Tool_Tag => $Tool_Tag ) unless( defined $rg_OB->{rID_prefix} );
        # count in this paired reads for RG ins evalue
        $rg_OB->load_reads_for_ins_evalue( reads_OB => $reads_OB, only_SoftClip => $only_SoftClip, pr_AimCount => $pr_AimCount, pr_DropProb => $pr_DropProb );
        # stop when all RG(s) are full
        $rgNeedReads-- if( $rg_OB->{EnoughReadsBool} );
        last if( $rgNeedReads == 0 );
    }
    close BAM;

    # inform
    stout_and_sterr "[INFO]\tLoad paired-end info for reads group (only_SoftClip:$only_SoftClip) from bam file.\n"
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

#--- extract reads group and create related objects ---
sub extract_ReadsGroup_OB{

    my $bam = shift;
    my %parm = @_;
    my $rgid2rgOB_Href = $parm{rgid2rgOB_Href};
    my $samtools       = $parm{samtools};
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
                warn_and_exit "<ERROR>:\tCannot get RG_ID or LB_ID from header line $. of  bam file:\n".
                                      "\t$bam->{filepath}\n";
            }
            # check recurrence
            if( exists $rgid2rgOB_Href->{$RG_ID} ){
                warn_and_exit "<ERROR>:\tThe RG_ID $RG_ID from header line $. of $bam->{tissue} bam file is recurrent.\n".
                                      "\t$bam->{filepath}\n";
            }
            # create reads group object (rg_OB)
            $rgid2rgOB_Href->{$RG_ID} = BioFuse::BioInfo::Objects::ReadsGroup_OB->new( bam => $bam, RG_ID => $RG_ID, LB_ID => $LB_ID );
            # link the rg_OB with the bam
            $bam->add_ReadsGroup_OBs( rgOB_Aref => [ $rgid2rgOB_Href->{$RG_ID} ] );
        }
    }
    close BAMHEADER;

    # once the bam lacks reads group
    if( scalar(keys %{$bam->{rgOB}}) == 0 ){
        if( 0 ){ # $Allow_Lack_RG
            # future updates 
        }
        else{
            warn_and_exit "<ERROR>:\tThe $bam->{tissue} bam file lacks reads group info.\n".
                                  "\t$bam->{filepath}\n";
        }
    }

    # inform
    stout_and_sterr "[INFO]\tCreate reads group objects from $bam->{tissue} bam file.\n"
                         ."\t$bam->{filepath}\n";
}

#--- add reads group object(s) to this bam ---
sub add_ReadsGroup_OBs{
    my $bam = shift;
    my %parm = @_;
    $bam->{rgOB}->{$_->get_RGid} = $_ for @{$parm{rgOB_Aref}};
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
    my $samtools = $parm{samtools};

    unless( $region_refseg && $region_st_pos && $region_ed_pos && $outFilePrefix){
        warn_and_exit "bam->get_region_depth( refseg=>?, st_pos=>?, ed_pos=>?, out_prefix=>?)\n";
    }

    unless( defined($samtools) && -e $samtools ){
        warn_and_exit "<ERROR>:\tCannot locate samtools in function bam->get_region_depth().\n";
    }

    my $bamForDepth = $bam->{filepath};
    my $viewRegParm = "$region_refseg:$region_st_pos-$region_ed_pos";

    # filter soft-clip
    if( $no_SoftClip ){
        my $tmp_bam = "$outFilePrefix.noSoftClip.bam";
        my $tmpbam_cmd = "( $samtools view $bam->{filepath} $viewRegParm | awk '\$6 !~ /S/' | $samtools view -S -b -o $tmp_bam - ) && ( $samtools index $tmp_bam )";
        trible_run_for_success($tmpbam_cmd, 'tmpNoSoftClipBam', {esdo_Nvb=>1});
        $bamForDepth = $tmp_bam;
    }

    # get depth
    $bam->{regionDepthFile} = "$outFilePrefix.depth.gz";
    my $depth_cmd = "( $samtools depth -a -r $viewRegParm -q $min_baseQ -Q $min_mapQ $bamForDepth | gzip -c > $bam->{regionDepthFile} )";
    trible_run_for_success($depth_cmd, 'bamDepth', {esdo_Nvb=>1});

    # sweep
    if( $bamForDepth ne $bam->{filepath} ){
        `rm -rf $bamForDepth`;
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
    my $samtools = $parm{samtools};
    my $bcftools = $parm{bcftools};
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

#--- do stat on marker at given position ---
## maker: barcode (10x) or PE-id (Hi-C, MinIon and PacBio)
sub get_pos_marker_stat{

    my $bam = shift;
    my %parm = @_;
    my $chr = $parm{chr};
    my $pos = $parm{pos};
    my $marker = $parm{marker}; # '10x', 'Hi-C', 'MinIon', 'PacBio'
    my $samtools = $parm{samtools};

    my $PmakerStat_Href = {};
    # FLAG: -F 0x100(sd) + 0x400(d) + 0x800(sp) = 0xD00
    open (BAM,"$samtools view -F 0xD00 $bam->{filepath} $chr:$pos-$pos |") || die"fail read bam: $!\n";
    while(<BAM>){
        # load reads object
        my $reads_OB = BioFuse::BioInfo::Objects::Reads_OB->new( ReadsLineText => $_, _rc_optfd => 1, _rc_rdstr => 1 );
        next if( $reads_OB->is_unmap );
        # get marker of pos
        my $marker_str = ( $marker =~ /10x/i ? $reads_OB->get_10x_barc : $reads_OB->get_pid );
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
    my $samtools = $parm{samtools};

    my $AmakerStat_Href = {};
    # FLAG: -F 0x100(sd) + 0x400(d) + 0x800(sp) = 0xD00
    open (BAM,"$samtools view -F 0xD00 $bam->{filepath} $chr:$pos-$pos |") || die"fail read bam: $!\n";
    while(<BAM>){
        # load reads object
        my $reads_OB = BioFuse::BioInfo::Objects::Reads_OB->new( ReadsLineText => $_, _rc_optfd => 1, _rc_rdstr => 1 );
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
            my $marker_str = ( $marker =~ /10x/i ? $reads_OB->get_10x_barc : $reads_OB->get_pid );
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
    my $samtools = $parm{samtools};
    my $mark = $parm{mark} || $bam->get_tag;
    my $viewOpt = $parm{viewOpt} || '';
    my $readsType = $parm{readsType} || 'norm'; # norm/HiC
    my $subrtRef = $parm{subrtRef};
    my $subrtParmAref = $parm{subrtParmAref} || [];
    my $peC_ReportUnit = $parm{peC_ReportUnit} || 1E6;
    my $peOB_AbufferSize = $parm{peOB_AbufferSize} || 1E5;
    my $deal_peOB_pool = $parm{deal_peOB_pool} || 0;
    my $quiet = $parm{quiet} || 0; # not inform PE Count
    my $simpleLoad = $parm{simpleLoad} || 0; # just record reads' basic information

    # objects modules
    my $rdObjSource = {norm => 'Reads_OB',   HiC => 'HicReads_OB'};
    my $rdObjModule = 'BioFuse::BioInfo::Objects::' . $rdObjSource->{$readsType};
    my $peObjSource = {norm => 'PairEnd_OB', HiC => 'HicPairEnd_OB'};
    my $peObjModule = 'BioFuse::BioInfo::Objects::' . $peObjSource->{$readsType};

    # read BAM
    my $pe_Count = 0;
    my $pe_CountInform = $peC_ReportUnit;
    my $last_pid = '';
    my $last_peOB = undef;
    my @peOB_pool = ();
    my $fh = $bam->start_read(samtools => $samtools, viewOpt => $viewOpt);
    while(<$fh>){
        my $reads_OB = $rdObjModule->new( ReadsLineText => $_, _rc_optfd => 1, _rc_rdstr => 1, simpleLoad => $simpleLoad );
        if( $reads_OB->get_pid ne $last_pid ){
            # deal last pe_OB
            $pe_Count++;
            push @peOB_pool, $last_peOB if( defined $last_peOB );
            if( $pe_Count % $peOB_AbufferSize == 0 ){
                if($deal_peOB_pool){
                    &{$subrtRef}(pe_OB_poolAf => \@peOB_pool, @$subrtParmAref);
                }
                else{
                    &{$subrtRef}(pe_OB => $_, @$subrtParmAref) for @peOB_pool;
                }
                # sweep
                @peOB_pool = ();
                # inform
                if($pe_Count >= $pe_CountInform){
                    stout_and_sterr "[INFO]\t".`date`
                                         ."\tload $pe_CountInform PE-reads from $mark bam.\n" unless $quiet;
                    $pe_CountInform += $peC_ReportUnit;
                }
            }
            # create new pe_OB
            $last_peOB = $peObjModule->new;
            # update
            $last_pid = $reads_OB->get_pid;
        }
        # load this reads_OB to pe_OB
        $last_peOB->load_reads_OB(reads_OB => $reads_OB);
    }
    close $fh;
    # deal the last one
    if( defined $last_peOB ){
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
                         ."\tTotally, load $pe_Count PE-reads from $mark bam.\n" unless $quiet;
}

1; ## tell the perl script the successful access of this module.