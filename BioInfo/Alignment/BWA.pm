package BioFuse::BioInfo::Alignment::BWA;

use strict;
use warnings;
use BioFuse::Util::Log qw/ cluck_and_exit stout_and_sterr /;
use BioFuse::Util::Sys qw/ file_exist trible_run_for_success /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              PEfqToSortBam
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Alignment::BWA';
#----- version --------
$VERSION = "0.03";
$DATE = '2022-03-06';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @function_list = qw/
                        PEfqToSortBam
                        FilterPE
                     /;

#--- BWA+samtools pipeline from fq to sort.bam ---
sub PEfqToSortBam{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $ref = $parm{ref};
    my $fqAf = $parm{fqAf}; # BioFuse::BioInfo::Objects::SeqData::FastQ_OB
    my $bam = $parm{bam}; # BioFuse::BioInfo::Objects::SeqData::Bam_OB
    my $bwa = $parm{bwa};
    my $samtools = $parm{samtools};
    # additional
    my $func = $parm{func} || 'mem'; # mem: >70bp, aln: <=70bp
    my $maxIns = $parm{maxIns} || 500; # for 'aln' func, '-a' para
    my $RG_para = $parm{RG_para} || ''; # '-R xxx' / '-r xxx'
    my $view_para = $parm{view_para} || ''; # '-F xxx -f xxx'
    my $discBuPair = $parm{discBuPair} || 0; # discard both unmaped PE-reads
    my $sort_mode = lc($parm{sort_mode} || 'c'); # [c]oordinates or reads-[n]ame
    my $onlyRawBam = $parm{onlyRawBam} || 0; # only get raw bam
    my $maxBiClipRatio = $parm{maxBiClipRatio} || undef; # Soft/Hart-clip part ratio
    my $whClipRatioTimes = $parm{whClipRatioTimes} || 2.5;

    # check
    for my $key (qw/ ref bwa samtools /){
        cluck_and_exit "<ERROR>\tcannot find $key.\n" unless file_exist(filePath => $parm{$key});
    }

    my $bamPath = $bam->filepath;
    my @fqPath = map {$_->filepath} @$fqAf;
    # PEfqToInitBam
    my $raw_bamPath = $bamPath.'.init.bam';
    my $raw_bam = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath => $raw_bamPath, tag => 'rawBam');
    $raw_bam->addTool(samtools => $samtools);
    ## Platform must be recognized one, or GATK cannot run
    ## -Y, soft-clip for supplementary alignments
    ## -p, smart pairing
    ## use awk to filter the bwa standout sam
    ##    awk '
    ##            $0~/^\@/                                    # include SAM header
    ##         || (    !(and($2,0x4) && and($2,0x8))          # Cannot be u and U simultaneously
    ##              && xor(and($2,0x40), and($2,0x80))        # must have and only have one of end_NO 1 and 2
    ##              && and($2,0x1)                            # must be paired sequencing, mean 'p'
    ##            )
    ##        '
    my $cmd;
    my $awk = "awk '\$0~/^\@/ || (!(and(\$2,0x4) && and(\$2,0x8)) && xor(and(\$2,0x40), and(\$2,0x80)) && and(\$2,0x1))'";
    if($func eq 'mem'){
        $cmd =  "$bwa mem -t 8 -Y $RG_para "
               .(@fqPath == 1 ? '-p ' : '')
               ."$ref @fqPath 2>/dev/null "
               .($discBuPair ? "| $awk " : '')
               ."| $samtools view -b -S $view_para -T $ref -o $raw_bamPath - 2>/dev/null";
    }
    else{
        my $sai_1 = $bamPath.'.e1.sai';
        my $sai_2 = $bamPath.'.e2.sai';
        $cmd =  "($bwa aln -L -m 800000 -l 27 -i 10 -k 3 -t 8 -e 10 -f $sai_1 $ref $fqPath[0] 2>/dev/null) && "
               ."($bwa aln -L -m 800000 -l 27 -i 10 -k 3 -t 8 -e 10 -f $sai_2 $ref $fqPath[1] 2>/dev/null) && "
               ."($bwa sampe -a $maxIns $RG_para $ref $sai_1 $sai_2 $fqPath[0] $fqPath[1] 2>/dev/null "
               .($discBuPair ? "| $awk " : '')
               ."| $samtools view -b -S $view_para -T $ref -o $raw_bamPath - 2>/dev/null) && "
               ."(rm -f $sai_1 $sai_2)";
    }
    trible_run_for_success($cmd, 'get raw_bam', {esdo_Nvb=>1});
    # inform
    stout_and_sterr "[INFO]\tBWA $func finished.\n"
                         ."\tRef: $ref\n";
    stout_and_sterr       "\tFq:  $_\n" for @fqPath;

    # only want raw bam
    if($onlyRawBam){
        `mv -f $raw_bamPath $bamPath`;
        return;
    }

    # toNsort
    my $nSort_bamPath = "$bamPath.nSort.bam";
    my $nSort_bam = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath => $nSort_bamPath, tag => "nSortBam");
    $nSort_bam->addTool(samtools => $samtools);
    $raw_bam->toNsort(nSortBam => $nSort_bam);
    $raw_bam->delete_file; # sweep
    # inform
    stout_and_sterr "[INFO]\tSamTools nSort finished.\n";

    # fixmate
    my $fixmate_bamPath = "$bamPath.fixmate.bam";
    my $fixmate_bam = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath => $fixmate_bamPath, tag => 'fixmateBam');
    $fixmate_bam->addTool(samtools => $samtools);
    $nSort_bam->toFixmate(fixmateBam => $fixmate_bam);
    $nSort_bam->delete_file; # sweep
    # inform
    stout_and_sterr "[INFO]\tSamTools fixmate finished.\n";

    # additional filter
    if(defined $maxBiClipRatio){
        my $filter_bamPath = "$bamPath.fixmate.filter.bam";
        my $filter_bam = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath => $filter_bamPath, tag => 'filterBam');
        $filter_bam->addTool(samtools => $samtools);
        # start writing
        $filter_bam->start_write;
        # header
        $filter_bam->write(content=>$_) for @{$fixmate_bam->header_Af};
        # filter reads
        my @Parm = (filterBam=>$filter_bam, maxBiClipRatio=>$maxBiClipRatio, whClipRatioTimes=>$whClipRatioTimes);
        $fixmate_bam->smartBam_PEread(deal_peOB_pool => 1, subrtRef => \&FilterPE, subrtParmAref => \@Parm, quiet => 1);
        # stop writing
        $filter_bam->stop_write;
        # rename to fixmate_bam
        `mv -f $filter_bamPath $fixmate_bamPath`;
    }

    # want nSort
    if($sort_mode eq 'n'){
        `mv -f $fixmate_bamPath $bamPath`;
        return;
    }

    # want cSort
    my $cSort_bam = BioFuse::BioInfo::Objects::SeqData::Bam_OB->new(filepath => $bamPath, tag => "cSortBam");
    $cSort_bam->addTool(samtools => $samtools);
    $fixmate_bam->toCsort(cSortBam => $cSort_bam);
    $cSort_bam->toIndex;
    $fixmate_bam->delete_file; # sweep
    # inform
    stout_and_sterr "[INFO]\tSamTools cSort finished.\n";
}

#--- addtional pe filter ---
sub FilterPE{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pe_OB_poolAf = $parm{pe_OB_poolAf};
    my $filterBam = $parm{filterBam}; # bam object
    my $maxBiClipRatio = $parm{maxBiClipRatio};
       $maxBiClipRatio = 1 unless defined $maxBiClipRatio;
    my $whClipRatioTimes = $parm{whClipRatioTimes} || 2.5;

    for my $pe_OB (@$pe_OB_poolAf){
        my $r1_OB = $pe_OB->rOB_Af(reads_end=>1)->[0];
        my $r2_OB = $pe_OB->rOB_Af(reads_end=>2)->[0];
        my $r1_unmap = $r1_OB->is_unmap;
        my $r2_unmap = $r2_OB->is_unmap;
        # r1 clip
        my $r1_beydClip = 0;
        if(!$r1_unmap){
            my $r1_foreClipLen = $r1_OB->foreClipLen;
            my $r1_hindClipLen = $r1_OB->hindClipLen;
            my $r1_mxCPL = $r1_OB->rlen * $maxBiClipRatio;
               $r1_beydClip   = ($r1_foreClipLen > $r1_mxCPL) && ($r1_hindClipLen > $r1_mxCPL);
               $r1_beydClip ||= ($r1_foreClipLen + $r1_hindClipLen > $r1_mxCPL * $whClipRatioTimes);
        }
        # r2 clip
        my $r2_beydClip = 0;
        if(!$r2_unmap){
            my $r2_foreClipLen = $r2_OB->foreClipLen;
            my $r2_hindClipLen = $r2_OB->hindClipLen;
            my $r2_mxCPL = $r2_OB->rlen * $maxBiClipRatio;
               $r2_beydClip   = ($r2_foreClipLen > $r2_mxCPL) && ($r2_hindClipLen > $r2_mxCPL);
               $r2_beydClip ||= ($r2_foreClipLen + $r2_hindClipLen > $r2_mxCPL * $whClipRatioTimes);
        }
        # clip filter
        if(    (!$r1_unmap && !$r2_unmap && $r1_beydClip && $r2_beydClip)
            || (!$r1_unmap &&  $r2_unmap && $r1_beydClip                )
            || ( $r1_unmap && !$r2_unmap                 && $r2_beydClip)
        ){
            next;
        }
        else{
            for my $end (1,2){
                $filterBam->write(content=>$_->printSAM."\n") for @{$pe_OB->rOB_Af(reads_end=>$end)};
            }
        }
    }
}

1; ## tell the perl script the successful access of this module.
