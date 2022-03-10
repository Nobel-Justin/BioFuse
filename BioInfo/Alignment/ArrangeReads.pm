package BioFuse::BioInfo::Alignment::ArrangeReads;

use strict;
use warnings;
use List::Util qw/ max first any /;
use BioFuse::Util::Log qw/ cluck_and_exit stout_and_sterr /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              ArrangeReadsAmongRefseg
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'ArrangeReads';
#----- version --------
$VERSION = "0.07";
$DATE = '2022-03-10';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @function_list = qw/
                        ArrangeReadsAmongRefseg
                        SelectRefsegForReads
                        MxasKeepRef_clipMatch
                     /;

#--- arrange reads among multiple refseg ---
# !!! the source bam should be Nsort (after merged by multiple refseg alignments)
sub ArrangeReadsAmongRefseg{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $nSortBam = $parm{nSortBam}; # bam object
    my $outBam = $parm{outBam}; # bam object
    my $onlyPropMap = $parm{onlyPropMap} || 0; # only properly mapped
    my $skipSuppMap = $parm{skipSuppMap} || 0; # skip supplement/second alignment
    my $skipDupRead = $parm{skipDupRead} || 0; # skip PCR duplicated reads
    my $skipPEunmap = $parm{skipPEunmap} || 0; # skip PE-reads that both unmap
    # refseg need to keep SQ of bam header
    ## for mapped reads, only keep alignment on this refseg
    my $keepRefsegID = $parm{keepRefsegID} || undef;
    # to allow the keepRefsegID is not the only refseg with max AS
    my $keepRefNsolo = $parm{keepRefNsolo} || 0;
    # keep the keepRefsegID as long as its AS satisfies this [m]in [r]atio [d]ifference from the max-AS
    my $keepRefASmrd = $parm{keepRefASmrd} || 0;
    # swap when PE ends map to diff-refseg respectively
    my $swapAlignEnds = $parm{swapAlignEnds} || 0;
    # !! effective when keepRefsegID is not provided
       $swapAlignEnds = 0 if defined $keepRefsegID;
    # for one end, if has multiple alignment, save remains to supplementary
    my $moreMapToSupp = $parm{moreMapToSupp} || 0;
    # !! effective when keepRefsegID is not provided
       $moreMapToSupp = 0 if defined $keepRefsegID;
    my $toSuppMinMapQ = $parm{toSuppMinMapQ};
       $toSuppMinMapQ = 30 if !defined $toSuppMinMapQ;
    my $toSuppMAminCL = $parm{toSuppMAminCL}; # minimum clipped length of main alignment to have supp
       $toSuppMAminCL = 15 if !defined $toSuppMAminCL;

    # start writing
    $outBam->start_write;
    # header
    if(defined $keepRefsegID){
        $outBam->write(content=>$_) for grep !/^\@SQ/ || /\sSN:$keepRefsegID\s/, @{$nSortBam->header_Af};
    }
    else{
        $outBam->write(content=>$_) for @{$nSortBam->header_Af};
    }
    # arrange reads
    my $viewOpt = '';
       $viewOpt .= "-f 0x2" if $onlyPropMap;
       $viewOpt .= "-F 0x900" if $skipSuppMap; # -F 0x800 + 0x100
       $viewOpt .= "-F 0x400" if $skipDupRead;
    $nSortBam->smartBam_PEread(viewOpt => $viewOpt, deal_peOB_pool => 1, subrtRef => \&SelectRefsegForReads, subrtParmAref => \@_, quiet => 1);
    # stop writing
    $outBam->stop_write;
}

#--- arrange reads among multiple refseg ---
sub SelectRefsegForReads{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $pe_OB_poolAf = $parm{pe_OB_poolAf};
    my $outBam = $parm{outBam}; # bam object
    my $skipPEunmap = $parm{skipPEunmap} || 0; # skip PE-reads that both unmap
    # for mapped reads, only keep alignment on this refseg
    my $keepRefsegID = $parm{keepRefsegID} || undef;
    # to allow the keepRefsegID is not the only refseg with max AS
    my $keepRefNsolo = $parm{keepRefNsolo} || 0;
    # keep the keepRefsegID as long as its AS satisfies this [m]in [r]atio [d]ifference from the max-AS
    my $keepRefASmrd = $parm{keepRefASmrd} || 0;
    # swap when PE ends map to diff-refseg respectively
    my $swapAlignEnds = $parm{swapAlignEnds} || 0;
    # !! effective when keepRefsegID is not provided
       $swapAlignEnds = 0 if defined $keepRefsegID;
    # for one end, if has multiple alignment, save remains to supplementary
    my $moreMapToSupp = $parm{moreMapToSupp} || 0;
    # !! effective when keepRefsegID is not provided
       $moreMapToSupp = 0 if defined $keepRefsegID;
    my $toSuppMinMapQ = $parm{toSuppMinMapQ};
       $toSuppMinMapQ = 30 if !defined $toSuppMinMapQ;
    my $toSuppMAminCL = $parm{toSuppMAminCL}; # minimum clipped length of main alignment to have supp
       $toSuppMAminCL = 15 if !defined $toSuppMAminCL;

    for my $pe_OB (@$pe_OB_poolAf){
        my %rOB_Af = map{ ($_, $pe_OB->rOB_Af(reads_end=>$_)) } (1,2);
        # mapped rOB
        my %map_rOB = map { ($_, [grep !$_->is_unmap, @{$rOB_Af{$_}}]) } (1,2);
        my %endsMap = map { ($_, scalar(@{$map_rOB{$_}})) } (1,2);
        # different scenarios
        if( !($endsMap{1} or $endsMap{2}) ){ # _1 and _2 are all unmap
            next if $skipPEunmap;
            for my $end (1,2){ # unmap, just output one time (use first)
                $outBam->write(content=>$_->printSAM."\n") for first {$_->is_unmap} @{$rOB_Af{$end}};
            }
        }
        else{
            my %output;
            # alignment score of each refseg
            my %refseg2AS;
            $refseg2AS{$_->mseg} += $_->alignScore for map {@{$rOB_Af{$_}}} (1,2);
            # avail refseg with max alignment score
            my $max_AS = max(values %refseg2AS);
            my @max_AS_refseg = grep $refseg2AS{$_} == $max_AS, keys %refseg2AS;
            # if keepRefsegID is set
            if(defined $keepRefsegID){
                next if !exists $refseg2AS{$keepRefsegID} || $refseg2AS{$keepRefsegID} == 0; # no alignments on keepRefsegID
                my $skip = 0;
                # judge keepRefsegID
                if($keepRefNsolo){
                    # allow the keepRefsegID is not the only refseg with max AS
                    $skip = 1 if !(any {$_ eq $keepRefsegID}, @max_AS_refseg);
                }
                else{
                    # the keepRefsegID must be the only refseg with max AS
                    $skip = 1 if @max_AS_refseg != 1 || $max_AS_refseg[0] ne $keepRefsegID;
                }
                # try to SAVE keepRefsegID
                if(    $skip
                    && (   # arbitrarily SAVE keepRefsegID according to the keepRefASmrd
                           ($keepRefASmrd != 0 && $refseg2AS{$keepRefsegID} >= $max_AS * (1-$keepRefASmrd))
                        || # arbitrarily SAVE keepRefsegID if clipped part matches to that on refseg with max AS
                           &MxasKeepRef_clipMatch(map_rOB_Hf=>\%map_rOB, max_AS_refseg=>$max_AS_refseg[0], keepRefsegID=>$keepRefsegID)
                       )
                ){
                    $skip = 0;
                }
                # skip?
                next if $skip;
                # only keep keepRefsegID
                @max_AS_refseg = ($keepRefsegID);
            }
            # go on: different scenarios
            if( $endsMap{1} xor $endsMap{2} ){ # _1/_2, one has map (may have unmap simultaneously), another are all unmap
                # choose the paired r_OBs of the first refseg with max AS
                # $outBam->write(content=>$_->printSAM."\n")
                push @{$output{$_->endNO}}, $_
                    for grep {$_->mseg eq $max_AS_refseg[0] || $_->p_mseg eq $max_AS_refseg[0]}
                        map {@{$rOB_Af{$_}}} (1,2);
            }
            else{ # _1 and _2 both have map (may have unmap simultaneously)
                # firstly, find refseg that both _1 and _2 mapped
                my %end2mseg;
                $end2mseg{1}{$_->mseg} = 1 for @{$map_rOB{1}};
                $end2mseg{2}{$_->mseg} = 1 for @{$map_rOB{2}};
                ## sort by AS
                my @msegBiEnd = sort {$refseg2AS{$b}<=>$refseg2AS{$a}} grep exists $end2mseg{2}{$_}, keys %{$end2mseg{1}};
                # if keepRefsegID is set
                if(defined $keepRefsegID){
                    @msegBiEnd = grep $_ eq $keepRefsegID, @msegBiEnd;
                }
                # take msegBiEnd with most AS
                if(@msegBiEnd != 0){
                    # $outBam->write(content=>$_->printSAM."\n")
                    push @{$output{$_->endNO}}, $_
                        for grep {$_->mseg eq $msegBiEnd[0] || $_->p_mseg eq $msegBiEnd[0]}
                            map {@{$rOB_Af{$_}}} (1,2);
                }
                else{ # each refseg has just one end map, swap ends?
                    if(defined $keepRefsegID || !$swapAlignEnds){ # not allowed when keepRefsegID is set, or disabled
                        # $outBam->write(content=>$_->printSAM."\n")
                        push @{$output{$_->endNO}}, $_
                            for grep {$_->mseg eq $max_AS_refseg[0] || $_->p_mseg eq $max_AS_refseg[0]}
                                map {@{$rOB_Af{$_}}} (1,2);
                    }
                    else{ # do swap!
                        # find each end with max AS
                        my %maxAS_rOB = (1=>undef,2=>undef);
                        for my $end (1,2){
                            for my $rOB (@{$map_rOB{$end}}){
                                $maxAS_rOB{$end} = $rOB if !defined $maxAS_rOB{$end} || $maxAS_rOB{$end}->alignScore < $rOB->alignScore;
                            }
                        }
                        # swap: update rOB attributes, and output
                        for my $tEnd (1,2){ # [t]his end
                            my $pEnd = $tEnd % 2 + 1; # [p]air end
                            my $t_rOB = $maxAS_rOB{$tEnd};
                            my $p_rOB = $maxAS_rOB{$pEnd};
                            # update
                            ## paired end mapped refseg and position
                            my $p_mseg = $t_rOB->mseg eq $p_rOB->mseg ? '=' : $p_rOB->mseg;
                            $t_rOB->update_attr(attr_id=>'p_mseg', value=>$p_mseg);
                            $t_rOB->update_attr(attr_id=>'p_mpos', value=>$p_rOB->mpos);
                            ## set Tlen as zero
                            $t_rOB->update_attr(attr_id=>'tlen', value=>0);
                            ## flags
                            $t_rOB->  set_flag(flag=>0x1 );
                            $t_rOB->  set_flag(flag=>0x20) if $p_rOB->is_rv_map;
                            $t_rOB->unset_flag(flag=>0x20) if $p_rOB->is_fw_map;
                            $t_rOB->unset_flag(flag_Af=>[0x2,0x4,0x8,0x200,0x400,0x800]);
                            ## optfd
                            $t_rOB->update_optfd(tag=>'MQ:i:', value=>$p_rOB->mapQ);  # pair-end MQ
                            $t_rOB->update_optfd(tag=>'MC:Z:', value=>$p_rOB->cigar); # pair-end cigar
                            # output
                            # $outBam->write(content=>$t_rOB->printSAM."\n");
                            push @{$output{$tEnd}}, $t_rOB;
                        }
                    }
                }
            }
            # r_OB mapped to the other refseg (if has), to supplementary?
            if($moreMapToSupp){
                for my $tEnd (1,2){ # [t]his end
                    my $m_rOB = first {!$_->is_suppmap} @{$output{$tEnd}}; # selected [m]ain alignment of [t]his end
                    next if $m_rOB->biClipLen <= $toSuppMAminCL; # clipped len of [m]ain alignment
                    my $pEnd = $tEnd % 2 + 1; # [p]air end
                    my $p_rOB = first {!$_->is_suppmap} @{$output{$pEnd}}; # selected [m]ain alignment of [p]air end, maybe unmap
                    for my $t_rOB (@{$map_rOB{$tEnd}}){ # mapped
                        # should diff refseg of [m]ain alignment
                        next if $t_rOB->mseg eq $m_rOB->mseg;
                        # should same refseg of [p]air alignment when [p]air and [m]ain diff refseg
                        next if $m_rOB->mseg ne $p_rOB->mseg && $t_rOB->mseg ne $p_rOB->mseg;
                        # map qual
                        next if $t_rOB->mapQ != 255 && $t_rOB->mapQ < $toSuppMinMapQ;
                        # update
                        ## paired end mapped refseg and position
                        my $p_mseg = $t_rOB->mseg eq $p_rOB->mseg ? '=' : $p_rOB->mseg;
                        $t_rOB->update_attr(attr_id=>'p_mseg', value=>$p_mseg);
                        $t_rOB->update_attr(attr_id=>'p_mpos', value=>$p_rOB->mpos);
                        ## set Tlen as zero
                        $t_rOB->update_attr(attr_id=>'tlen', value=>0);
                        ## flags
                        $t_rOB->  set_flag(flag=>0x800); # set as supplementary
                        $t_rOB->  set_flag(flag=>0x1  );
                        $t_rOB->  set_flag(flag=>0x20 ) if  $p_rOB->is_rv_map;
                        $t_rOB->unset_flag(flag=>0x20 ) if  $p_rOB->is_fw_map;
                        $t_rOB->  set_flag(flag=>0x8  ) if  $p_rOB->is_unmap;
                        $t_rOB->unset_flag(flag=>0x8  ) if !$p_rOB->is_unmap;
                        $t_rOB->  set_flag(flag=>0x400) if  $m_rOB->is_dup; # copy dup status of [m]ain alignment
                        $t_rOB->unset_flag(flag=>0x400) if !$m_rOB->is_dup;
                        $t_rOB->unset_flag(flag_Af=>[0x2,0x4,0x200]);
                        ## optfd
                        $t_rOB->update_optfd(tag=>'MQ:i:', value=>$p_rOB->mapQ);  # pair-end MQ
                        $t_rOB->update_optfd(tag=>'MC:Z:', value=>$p_rOB->cigar); # pair-end cigar
                        # output
                        push @{$output{$tEnd}}, $t_rOB;
                    }
                }
            }
            # output
            $outBam->write(content=>$_->printSAM."\n") for map {@{$output{$_}}} (1,2);
        }
    }
}

#--- check whether the clipped part of keepRefsegID matches with that of max_AS refseg ---
sub MxasKeepRef_clipMatch{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $map_rOB_Hf = $parm{map_rOB_Hf};
    my $max_AS_refseg = $parm{max_AS_refseg};
    my $keepRefsegID = $parm{keepRefsegID};
    my $CLsumMinRatio = $parm{CLsumMinRatio} || 0.8;

    for my $end (1,2){
        next if scalar(@{$map_rOB_Hf->{$end}}) == 0;
        my $mxasRef_rOB = first {$_->mseg eq $max_AS_refseg} @{$map_rOB_Hf->{$end}};
        my $keepRef_rOB = first {$_->mseg eq $keepRefsegID } @{$map_rOB_Hf->{$end}};
        next if !defined $mxasRef_rOB || !defined $keepRef_rOB;
        return 1 if $mxasRef_rOB->clipMatch(other_rOB=>$keepRef_rOB, CLsumMinRatio=>$CLsumMinRatio);
    }
    # last, not match
    return 0;
}

1; ## tell the perl script the successful access of this module.
