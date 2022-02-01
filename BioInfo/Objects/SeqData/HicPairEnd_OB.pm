package BioFuse::BioInfo::Objects::SeqData::HicPairEnd_OB;

use BioFuse::BioInfo::Objects::SeqData::PairEnd_OB; # inheritance

use strict;
use warnings;
use Data::Dumper;
use BioFuse::Util::Log qw/ warn_and_exit /;
use BioFuse::Util::Array qw/ binarySearch /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter BioFuse::BioInfo::Objects::SeqData::PairEnd_OB);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Objects::SeqData::HicPairEnd_OB';
#----- version --------
$VERSION = "0.08";
$DATE = '2021-08-14';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @function_list = qw/
                        new
                        get_rEndWholeAlignJudge
                        get_rEndWholeSuppHaplo
                        isInValidPair
                        testLinkRefSeg
                        dEndSameHapJudge
                        sEndSoloHapJudge
                        sEndInterHapJudge
                        dEndInterHapJudge
                        addHapIDtoReadsOptfd
                     /;

#--- structure of object
# basis, BioFuse::BioInfo::Objects::SeqData::PairEnd_OB

#--- construction of object ---
sub new{
    my $type = shift;

    my $pe_OB = BioFuse::BioInfo::Objects::SeqData::PairEnd_OB->new( @_ );

    bless($pe_OB);
    return $pe_OB;
}

#--- return whole alignment judge string of given rEnd ---
sub get_rEndWholeAlignJudge{
    my $pe_OB = shift;
    my %parm = @_;

    return join( '', map {
                            ($_->get_AlignJudge)
                         }
                     @{ $pe_OB->rOB_Af( reads_end => $parm{rEnd} ) }
                );
}

#--- return whole HashRef containing supported haplo-id of given rEnd ---
sub get_rEndWholeSuppHaplo{
    my $pe_OB = shift;
    my %parm = @_;

    my %rWholeHap;
    for my $reads_OB ( @{ $pe_OB->rOB_Af( reads_end => $parm{rEnd} ) } ){
        my $rAlignHapHref = $reads_OB->get_SuppHaploHref;
        push @{ $rWholeHap{$_} }, @{ $rAlignHapHref->{$_} } for keys %$rAlignHapHref;
    }
    return \%rWholeHap;
}

#--- test is the paired-ends are invalid Hi-C PE ---
## 1) aligned too close, if set
## 2) aligned to same fragement
sub isInValidPair{
    my $pe_OB = shift;
    my %parm = @_;
    my $chr2enzymePosHf = $parm{chr2enzymePosHf};
    my $maxCloseAlignDist = $parm{maxCloseAlignDist} || 1E3;
    my $skipCloseAlign  = $parm{skipCloseAlign} || 0;

    my $rOB_sortAref = $pe_OB->sorted_rOB_Af(rEndAref => [1,2], onlyMap => 1);
    my $rOB_a = $rOB_sortAref->[0];
    my $rOB_b = $rOB_sortAref->[-1];
    if($rOB_a->mseg ne $rOB_b->mseg){
        return 0;
    }
    else{
        if(      # check close alignment
               (    ! $skipCloseAlign
                 && $rOB_a->is_closeAlign(test_rOB => $rOB_b, distance => $maxCloseAlignDist)
               )
            ||   # check fragment index
               (    binarySearch(query => $rOB_a->mpos, array => $chr2enzymePosHf->{$rOB_a->mseg})
                 == binarySearch(query => $rOB_b->mpos, array => $chr2enzymePosHf->{$rOB_b->mseg})
               )
        ){
            my $fw_a = $rOB_a->is_fw_map;
            my $rv_b = $rOB_b->is_rv_map;
            return 'DanglingEnd'  if ($fw_a && $rv_b); # most
            my $rv_a = $rOB_a->is_rv_map;
            my $fw_b = $rOB_b->is_fw_map;
            return 'SelfCircle'   if ($rv_a && $fw_b);
            return 'DumpedPairFw' if ($fw_a && $fw_b);
            return 'DumpedPairRv' if ($rv_a && $rv_b);
        }
        else{
            return 0;
        }
    }
}

#--- test whether links mapped to required chr-[pair] ---
# return 1 means mapped, return 0 means not mapped.
sub testLinkRefSeg{
    my $pe_OB = shift;
    my %parm = @_;
    my $soloBoth = $parm{soloBoth}; # pe must both match a_refseg if only a_refseg provided
    my $chrSortHref = $parm{chrSortHref} || undef;
    my $chrSortKey  = $parm{chrSortKey} || undef;

    # parameters
    unless(exists $parm{a_refSegHref}){
        warn_and_exit "<ERROR>\tlacks 'a_refSegHref' parameter in func 'testLinkRefSeg' of pe_OB.\n";
    }
    my %refSegHref = (a => $parm{a_refSegHref});
    $refSegHref{b} = $parm{b_refSegHref} if exists $parm{b_refSegHref};
    # test refseg
    my $rOB_Af = $pe_OB->sorted_rOB_Af(chrSortHref => $chrSortHref, chrSortKey  => $chrSortKey);
    my %testResult;
    for my $i (0,-1){
        my $mseg = $rOB_Af->[$i]->mseg;
        for my $ab (keys %refSegHref){
            if(exists $refSegHref{$ab}{$mseg}){
                $testResult{$i}{$ab}{isNeed} = 1;
                $testResult{$i}{$ab}{noNeed} = 0;
            }
            else{
                $testResult{$i}{$ab}{isNeed} = 0;
                $testResult{$i}{$ab}{noNeed} = 1;
            }
        }
    }
    # make judgment
    if(exists $refSegHref{b}){
        if(    $testResult{0}{a}{isNeed} * $testResult{-1}{b}{isNeed} != 0
            || $testResult{0}{b}{isNeed} * $testResult{-1}{a}{isNeed} != 0
        ){
            return 1;
        }
        else{
            return 0;
        }
    }
    else{
        if($soloBoth){
            if($testResult{0}{a}{isNeed} * $testResult{-1}{a}{isNeed} != 0){
                return 1;
            }
            else{
                return 0;
            }
        }
        else{
            if($testResult{0}{a}{isNeed} + $testResult{-1}{a}{isNeed} != 0){
                return 1;
            }
            else{
                return 0;
            }
        }
    }
}

#--- judgement on dEnd-hx PE-reads ---
## 1) at least one end must have 'UK' alignment, then go to 2)
## 2) at least two 'hx' alignments (R1-vs-R2) is not close aligned
## *) here, we'd better also check whether one UK is close to 'hx' alignment of the other end.
## return 0, go to phMut-sEnd-hx
## return 1, keep  phMut-dEnd-hx
sub dEndSameHapJudge{
    my $pe_OB = shift;
    my %parm = @_;
    my $maxCloseAlignDist = $parm{maxCloseAlignDist};

    my $R1_rOB_Aref = $pe_OB->rOB_Af( reads_end => 1 );
    my $R2_rOB_Aref = $pe_OB->rOB_Af( reads_end => 2 );

    # at least one end must have 'UK' alignment
    # else, keep dEnd
    my $R1_hasUK = grep ! $_->has_SuppHaplo, @$R1_rOB_Aref;
    my $R2_hasUK = grep ! $_->has_SuppHaplo, @$R2_rOB_Aref;
    return 1 unless( $R1_hasUK || $R2_hasUK );

    # then, check whether the 'hx' alignment of two ends is close-Align
    for my $r1_OB ( grep $_->has_SuppHaplo, @$R1_rOB_Aref ){ # skip 'UK'
        for my $r2_OB ( grep $_->has_SuppHaplo, @$R2_rOB_Aref ){ # skip 'UK'
            # not close !
            # return to keep dEnd
            if( ! $r1_OB->is_closeAlign(test_rOB => $r2_OB, distance => $maxCloseAlignDist) ){
                return 1;
            }
        }
    }
    # all close!
    return 0;
}

#--- judgement on sEnd-hx PE-reads ----
## 1)     the has-solo-hap rEnd cannot have two alignments supporting the same hapID
## 2) OR, the two alignments supporting the same hapID are not close aligned
## return 0, go to phMut-dEnd-hx
## return 1, keep  phMut-sEnd-hx
sub sEndSoloHapJudge{
    my $pe_OB = shift;
    my %parm = @_;
    my $shEnd = $parm{shEnd}; # solo-haplo
    my $maxCloseAlignDist = $parm{maxCloseAlignDist};

    my @shEnd_HasHap_rOB = grep $_->has_SuppHaplo, @{$pe_OB->rOB_Af(reads_end => $shEnd)};
    if(      @shEnd_HasHap_rOB >= 2
        && ! $shEnd_HasHap_rOB[0]->is_closeAlign(test_rOB => $shEnd_HasHap_rOB[-1], distance => $maxCloseAlignDist)
    ){
        return 0;
    }
    else{
        return 1;
    }
}

#--- judgement on sEnd-hInter PE-reads ---
##  1) have two alignments (i.e., SP) that interChr or large-distance 'heuristic';
##  2) one of two alignments is close to the other end's alignment
## *3) might enzyme in its sequence (not in use).
## return 1 means toSelectOneHap, 0 means toKeepMultiHap
sub sEndInterHapJudge{
    my $pe_OB = shift;
    my %parm = @_;
    my $mhEnd = $parm{mhEnd}; # multi-haplo
    my $uhEnd = $parm{uhEnd}; # unkonwn-haplo
    my $maxCloseAlignDist = $parm{maxCloseAlignDist};

    my $mh_rOB_Aref = $pe_OB->rOB_Af(reads_end => $mhEnd);
    my $uh_rOB_Aref = $pe_OB->rOB_Af(reads_end => $uhEnd); # much likely to have only one alignment
    my $doSelectBool = 0; # might select
    if(    scalar(@$mh_rOB_Aref) > 1 # have SP
        && $mh_rOB_Aref->[0]->has_SuppHaplo # both alignments have supported haplotype
        && $mh_rOB_Aref->[1]->has_SuppHaplo
    ){
        # pairwise alignments' distance
        my %closeAlign  = ( u1m1 => $uh_rOB_Aref->[0]->is_closeAlign(test_rOB => $mh_rOB_Aref->[0], distance => $maxCloseAlignDist),
                            m1m2 => $mh_rOB_Aref->[0]->is_closeAlign(test_rOB => $mh_rOB_Aref->[1], distance => $maxCloseAlignDist),
                            u1m2 => $uh_rOB_Aref->[0]->is_closeAlign(test_rOB => $mh_rOB_Aref->[1], distance => $maxCloseAlignDist)  );
        if(    ! $closeAlign{m1m2} # mh: not close
            && ( $closeAlign{u1m1} || $closeAlign{u1m2} ) # mh-uh: either close
        ){
            # return to keep MultiHap
            return 0;
        }
        else{ # select one haplo
            return 1;
        }
    }
    else{ # select one haplo
        return 1;
    }
}

#--- judgement on dEnd-hInter PE-reads ---
## 1) alignments support same haplo  should be close;
## 2) alignments support diff haplos should be not close;
## return 1 means toSelectOneHap, 0 means toKeepMultiHap
sub dEndInterHapJudge{
    my $pe_OB = shift;
    my %parm = @_;
    my $maxCloseAlignDist = $parm{maxCloseAlignDist};

    # prepare hapID -> reads_OB
    my %hapID2rOB;
    for my $rEnd (1, 2){
        my $rOB_Aref = $pe_OB->rOB_Af(reads_end => $rEnd);
        for my $rOB ( @$rOB_Aref ){
            for my $hapID (keys %{$rOB->get_SuppHaploHref}){ # auto-skip 'UK'
                push @{ $hapID2rOB{$hapID} }, $rOB;
            }
        }
    }
    # judgement
    my @hapID = sort keys %hapID2rOB;
    for my $h_i ( 0 .. $#hapID ){
        my $hapID_i = $hapID[$h_i];
        my $rOB_iA = $hapID2rOB{$hapID_i};
        my $rOB_iC = scalar(@$rOB_iA);
        # same hap should be close
        for my $r_i ( 0 .. $rOB_iC-1 ){
            for my $r_j ( $r_i+1 .. $rOB_iC-1 ){
                # not close !
                unless( $rOB_iA->[$r_i]->is_closeAlign(test_rOB => $rOB_iA->[$r_j], distance => $maxCloseAlignDist) ){
                    return 1;
                }
            }
        }
        # diff haps should not be close
        for my $h_j ( $h_i+1 .. $#hapID ){
            my $hapID_j = $hapID[$h_j];
            my $rOB_jA = $hapID2rOB{$hapID_j};
            my $rOB_jC = scalar(@$rOB_jA);
            for my $r_i ( 0 .. $rOB_iC-1 ){
                for my $r_j ( 0 .. $rOB_jC-1 ){
                    # close !
                    if( $rOB_iA->[$r_i]->is_closeAlign(test_rOB => $rOB_jA->[$r_j], distance => $maxCloseAlignDist) ){
                        return 1;
                    }
                }
            }
        }
    }
    # pass judgement
    # return to keep MultiHap
    return 0;
}

#--- add OWN 'XH:Z:' haplo-id to optfd of given rEnd's reads_OBs ---
sub addHapIDtoReadsOptfd{
    my $pe_OB = shift;
    my %parm = @_;
    $_->addHapIDtoOptfd for @{ $pe_OB->rOB_Af( reads_end => $parm{reads_end} ) };
}

1; ## tell the perl script the successful access of this module.
