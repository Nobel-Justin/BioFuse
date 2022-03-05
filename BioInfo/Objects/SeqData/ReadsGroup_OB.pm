package BioFuse::BioInfo::Objects::SeqData::ReadsGroup_OB;

use strict;
use warnings;
use File::Basename qw/ dirname /;
use List::Util qw/ max /;
use BioFuse::Util::Log qw/ cluck_and_exit stout_and_sterr /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Dist::DistStat qw/ engineer_Ntimes_SD_evaluation /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Objects::SeqData::ReadsGroup_OB';
#----- version --------
$VERSION = "0.05";
$DATE = '2019-06-03';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @function_list = qw/
                        new
                        RG_ID
                        LB_ID
                        RG_NO
                        set_RG_NO
                        rID_prefix
                        set_rID_prefix
                        maxRlen_Hf
                        update_maxRlen
                        affix
                        addAffix
                        bam
                        RG_para
                        load_reads_for_ins_evalue
                        evalue_ins
                        test_3p_overlap
                        file_prefix
                        report_format
                        write_report
                        load_report
                        write_insDistLog
                     /;

#--- structure of object
# rg -> RG_ID = $RG_ID
# rg -> RG_NO = $RG_NO
# rg -> LB_ID = $LB_ID
# rg -> bam = $bam
# rg -> tissue = $tissue
# rg -> rID_prefix = $rID_prefix
# rg -> maxRlen = {1=>,2=>}
# rg -> affix = {Cand_nSortBam=>bam_OB, Cand_smartFqGz=>fq_OB};

#--- construction of object ---
sub new{
    my $type = shift;
    my %parm = @_;

    my $rg = {};
    $rg->{RG_ID}  = $parm{RG_ID};
    $rg->{LB_ID}  = $parm{LB_ID};
    $rg->{RG_NO} = undef;
    $rg->{affix} = {};
    # reads
    $rg->{maxRlen} = {1=>0, 2=>0};
    $rg->{rID_prefix} = undef;
    $rg->{rID_prefix_prev} = undef;
    # bam
    if(exists $parm{bam}){
        $rg->{bam} = $parm{bam};
        $rg->{tissue} = $parm{bam}->tissue || 'U'; # Unknown
    }
    # for insert size evalue
    $rg->{Ins_Stat} = {};
    $rg->{pr_forIns_count} = 0;
    $rg->{pr_forIns_count_eng3sd} = 0; # after filtration of the 'eng-3sd'
    $rg->{EnoughReadsBool} = 0; # for stop bam-loading
    $rg->{pr_info_forIns} = {};
    $rg->{pr_maps_forDup} = {};
    $rg->{end_3p_overlap} = 0; # firstly count number, then turns to BOLL
    # insert size info
    $rg->{eligible} = 0;
    $rg->{Ins_mean} = 0;
    $rg->{Ins_sd} = 0;
    $rg->{Ins_peak} = {};
    $rg->{Ins_median} = {};
    # files
    ## initial is postfix
    $rg->{reportPostfix} = 'report';
    $rg->{report} = '';
    $rg->{insDistLogPostfix} = 'insDist.log';
    $rg->{insDistLog} = '';

    bless($rg);
    return $rg;
}

#--- return reads group ID ---
sub RG_ID{
    my $rg = shift;
    return $rg->{RG_ID};
}

#--- return library ID ---
sub LB_ID{
    my $rg = shift;
    return $rg->{LB_ID};
}

#--- return reads group NO ---
sub RG_NO{
    my $rg = shift;
    return $rg->{RG_NO};
}

#--- set RG NO. ---
sub set_RG_NO{
    my $rg = shift;
    my %parm = @_;
    $rg->{RG_NO} = $parm{RG_NO};
}

#--- return rID_prefix ---
sub rID_prefix{
    my $rg = shift;
    return $rg->{rID_prefix};
}

#--- set rID_prefix ---
sub set_rID_prefix{
    my $rg = shift;
    my %parm = @_;
    my $toolTag = $parm{toolTag} || 'BF'; # default as 'BioFuse'
    $rg->{rID_prefix} =  "_$toolTag"
                        .($rg->{tissue}?"ts$rg->{tissue}":'')
                        ."rgNO$rg->{RG_NO}"
                        ."${toolTag}_";
}

#--- return Hash-ref of max rlen ---
sub maxRlen_Hf{
    my $rg = shift;
    return $rg->{maxRlen};
}

#--- update max rlen ---
sub update_maxRlen{
    my $rg = shift;
    my %parm = @_;
    $rg->{maxRlen}->{$parm{endNO}} = max($rg->{maxRlen}->{$parm{endNO}}, $parm{rLen});
}

#--- return affix of given key ---
sub affix{
    my $rg = shift;
    my %parm = @_;
    return $rg->{affix}->{$parm{key}};
}

#--- add affix ---
sub addAffix{
    my $rg = shift;
    my %parm = @_;
    $rg->{affix}->{$_} = $parm{$_} for keys %parm;
}

#--- return bam object ---
sub bam{
    my $rg = shift;
    return $rg->{bam};
}

#--- return RG para used in BWA sampe/mem ---
sub RG_para{
    my $rg = shift;
    my %parm = @_;
    return '@RG\tID:'.$rg->RG_ID.'\tLB:'.$rg->LB_ID.'\tPL:ILLUMINA\tCN:NULL\tPU:'.$rg->RG_ID.'\tSM:'.$parm{sampleID};
}

#--- extract and record new read-id FS prefix ---
sub load_reads_for_ins_evalue{

    my $rg = shift;
    my %parm = @_;
    my $reads_OB       = $parm{reads_OB};
    my $only_SoftClip  = $parm{only_SoftClip} || 0; # at least one end
    my $allow_SoftClip = $parm{allow_SoftClip} || $only_SoftClip;
    my $pr_AimCount    = $parm{pr_AimCount};
    my $pr_DropProb    = $parm{pr_DropProb} || 0;

    # when only softclip is required
    if(     $only_SoftClip
        && !$reads_OB->is_softclip
      ){
        delete $rg->{pr_info_forIns}->{$reads_OB->{pid}}; # just delete
        return;
    }

    # the fw-mapped end has been recorded before
    if( exists $rg->{pr_info_forIns}->{$reads_OB->{pid}} ){
        # check mapping orientation
        unless( $reads_OB->is_rv_map ){
            cluck_and_exit "<ERROR>\tMapping orientation mistake of properly-mapped paired-end $reads_OB->{pid}.\n"
                                ."\t$rg->{bam}->{filepath}\n";
        }
        # for check duplication
        my $preads_OB = $rg->{pr_info_forIns}->{$reads_OB->{pid}};
        my $pr_maps_info = join(';', $preads_OB->{mseg}, $preads_OB->{mpos}, $reads_OB->{mseg}, $reads_OB->{mpos});
        # deal with this end
        if(    $reads_OB->is_good_cigar
            && ( $allow_SoftClip || ! $reads_OB->is_softclip )
            && !exists $rg->{pr_maps_forDup}->{$pr_maps_info}
        ){
            # mark for dedup
            $rg->{pr_maps_forDup}->{$pr_maps_info} = 1;
            # fix and record Tlen
            my $Tlen = $preads_OB->{tlen} + ( ($allow_SoftClip) ? $reads_OB->tlen_FixTlen_with_S : 0 );
            $rg->{Ins_Stat}->{$Tlen} ++;
            # update maximum read length
            $reads_OB->update_rgOB_maxRlen;
            $preads_OB->update_rgOB_maxRlen;
            # overlap at the 3' end ?
            $rg->{end_3p_overlap} += ($Tlen > ($reads_OB->{rlen}+$preads_OB->{rlen})) ? 0 : 1;
            # count in
            if( ++$rg->{pr_forIns_count} == $pr_AimCount ){
                $rg->{EnoughReadsBool} = 1;
                undef $rg->{pr_info_forIns};
                undef $rg->{pr_maps_forDup};
                return;
            }
        }
        # never use it, bye-bye
        delete $rg->{pr_info_forIns}->{$reads_OB->{pid}};
    }
    # to record this fw-mapped end
    elsif(   $reads_OB->is_fw_map
          && $reads_OB->is_good_cigar
          && ( $allow_SoftClip || ! $reads_OB->is_softclip )
          && ( !$pr_DropProb || rand(1) > $pr_DropProb ) # random drop
    ){
        $reads_OB->{tlen} += $reads_OB->tlen_FixTlen_with_S if( $allow_SoftClip ); # count in fore soft-clip of reads
        $rg->{pr_info_forIns}->{$reads_OB->{pid}} = $reads_OB;
    }
}

#--- evalue insert size ---
sub evalue_ins{

    my $rg = shift;
    my %parm = @_;
    my $MinPairForInsEvalue = $parm{MinPairForInsEvalue} || 1E3;

    # check pair count
    if( $rg->{pr_forIns_count} < $MinPairForInsEvalue ){
        $rg->{eligible} = 0;
        stout_and_sterr "[WARN]\tRG $rg->{RG_ID} has paired-end less than $MinPairForInsEvalue for ins evaluation. So ignore it.\n";
        return;
    }
    else{
        $rg->{eligible} = 1;
    }

    # insert size evaluation, default: N=3
    my $Dist_Href = engineer_Ntimes_SD_evaluation( Value2Count_Href => $rg->{Ins_Stat}, ratio_digit => 0 );
    $rg->{Ins_mean} = $Dist_Href->{value_mean};
    $rg->{Ins_sd}   = $Dist_Href->{value_sd};
    $rg->{Ins_peak} = $Dist_Href->{value_peak_Href};
    $rg->{Ins_median} = $Dist_Href->{value_median_Href};
    $rg->{pr_forIns_count_eng3sd} = $Dist_Href->{count_sum};
}

#--- whether certain pr of this RG are PCW ---
sub test_3p_overlap{
    my $rg = shift;
    my %parm = @_;
    my $perct_cutoff = $parm{perct_cutoff} || 0.2;

    $rg->{end_3p_overlap} = (   $rg->{eligible}
                             && $rg->{end_3p_overlap} / $rg->{pr_forIns_count} > $perct_cutoff
                            ) ? 1 : 0;
}

#--- get prefix of stat file of RG ---
sub file_prefix{
    my $rg = shift;
    return "$rg->{tissue}.RG_NO-$rg->{RG_NO}";
}

#--- structure to output report ---
# type = text, just output
# type = bool, yes:no
# type = rbool, no:yes
sub report_format{
    return  [   # tag                type and hierarchical keys
                [ 'Tissue',          { type => 'text', keyAref => ['tissue']                 } ],
                [ 'RG_ID',           { type => 'text', keyAref => ['RG_ID']                  } ],
                [ 'RG_NO',           { type => 'text', keyAref => ['RG_NO']                  } ],
                [ 'Library-ID',      { type => 'text', keyAref => ['LB_ID']                  } ],
                [ 'Eligible',        { type => 'bool', keyAref => ['eligible']               } ],
                [ '#_of_colt_PE',    { type => 'text', keyAref => ['pr_forIns_count']        } ],
                [ '#_of_used_PE',    { type => 'text', keyAref => ['pr_forIns_count_eng3sd'] } ],
                [ 'MaxRlen-end1',    { type => 'text', keyAref => ['maxRlen', '1']          } ],
                [ 'MaxRlen-end2',    { type => 'text', keyAref => ['maxRlen', '2']          } ],
                [ 'Ins-mean',        { type => 'text', keyAref => ['Ins_mean']               } ],
                [ 'Ins-sd',          { type => 'text', keyAref => ['Ins_sd']                 } ],
                [ 'Ins-peakValue',   { type => 'text', keyAref => ['Ins_peak', 'value']      } ],
                [ 'Ins-peakCount',   { type => 'text', keyAref => ['Ins_peak', 'count']      } ],
                [ 'Ins-medianValue', { type => 'text', keyAref => ['Ins_median', 'value']    } ],
                [ 'Ins-medianCount', { type => 'text', keyAref => ['Ins_median', 'count']    } ],
                [ 'rid-prefixCurr',  { type => 'text', keyAref => ['rID_prefix']             } ], 
                [ 'rid-prefixPrev',  { type => 'text', keyAref => ['rID_prefix_prev']        } ], 
                [ 'prime3overlap',   { type => 'bool', keyAref => ['end_3p_overlap']         } ],
                [ 'BamFile',         { type => 'text', keyAref => ['bam', 'filepath']     } ]
            ];
}

#--- write RG report ---
sub write_report{
    my $rg = shift;
    my %parm = @_;
    my $folder = $parm{folder};

    # aim file
    my $rgStatFilePrefix = $rg->file_prefix;
    $rg->{report} = File::Spec->catfile($folder, "$rgStatFilePrefix.$rg->{reportPostfix}");
    # structure of report
    my $structure_Aref = $rg->report_format;
    # output
    open (RGRT, Try_GZ_Write($rg->{report})) || die "fail generate RG report: $!\n";
    for my $lineinfo_Aref (@$structure_Aref){
        my $tag = $lineinfo_Aref->[0];
        my $type = $lineinfo_Aref->[1]->{type};
        my $hkeys_Aref = $lineinfo_Aref->[1]->{keyAref};
        my $value = $rg;
        $value = $value->{ $hkeys_Aref->[$_] } for ( 0 .. scalar(@$hkeys_Aref)-1 );
        $value = $value ? 'yes' : 'no' if( $type eq 'bool' );
        print RGRT join("\t", $tag, $value)."\n";
    }
    close RGRT;
}

#--- load RG report to reconstruct RG_OB ---
sub load_report{
    my $rg = shift;
    my %parm = @_;
    $rg->{report} = $parm{rg_report};

    # read report
    open (RGRT, Try_GZ_Read($rg->{report})) || die "fail read RG report: $!\n";
    my %info_pool = map{ (split /[\s\:]+/)[0,1] } <RGRT>;
    close RGRT;
    # structure of report
    my $structure_Aref = &report_format;
    for my $lineinfo_Aref (@$structure_Aref){
        my $tag = $lineinfo_Aref->[0];
        unless( exists $info_pool{$tag} ){
            cluck_and_exit "<ERROR>\tMiss tag $tag in RG report.\n"
                                ."\t$rg->{report}\n";
        }
        my $value = $info_pool{$tag};
        my $type = $lineinfo_Aref->[1]->{type};
        $value = $value eq 'yes' ? 1 : 0 if( $type eq 'bool' );
        my $hkeys_Aref = $lineinfo_Aref->[1]->{keyAref};
        my $AimAttribute = $rg;
        $AimAttribute = $AimAttribute->{ $hkeys_Aref->[$_] } for ( 0 .. scalar(@$hkeys_Aref)-2 );
        $AimAttribute->{ $hkeys_Aref->[-1] } = $value;
    }
    # remains
    my $rgStatFilePrefix = $rg->file_prefix;
    $rg->{insDistLog} = File::Spec->catfile( dirname($rg->{report}), "$rgStatFilePrefix.$rg->{insDistLogPostfix}" );
    # file check
    unless( -e $rg->{insDistLog} ){
        cluck_and_exit "<ERROR>\tCannot find the insert size distribution file of RG $rg->{RG_ID}\n";
    }
}

#--- write distribution of ins ---
sub write_insDistLog{

    my $rg = shift;
    my %parm = @_;
    my $folder = $parm{folder};

    # aim file
    my $rgStatFilePrefix = $rg->file_prefix;
    $rg->{insDistLog} = File::Spec->catfile($folder, "$rgStatFilePrefix.$rg->{insDistLogPostfix}");
    open (DIST, Try_GZ_Write($rg->{insDistLog})) || die "fail generate ins Dist of RG: $!\n";
    print DIST "#ins\tcount\n";
    print DIST "$_\t$rg->{Ins_Stat}->{$_}\n" for sort {$a<=>$b} keys %{$rg->{Ins_Stat}};
    close DIST;
    # release memory
    undef $rg->{Ins_Stat};
}

1; ## tell the perl script the successful access of this module.
