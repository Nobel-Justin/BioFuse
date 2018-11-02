package BioFuse::BioInfo::Objects::HicPairEnd_OB;

use BioFuse::BioInfo::Objects::PairEnd_OB; # inheritance

use strict;
use warnings;
use Data::Dumper;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter BioFuse::BioInfo::Objects::PairEnd_OB);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Objects::HicPairEnd_OB';
#----- version --------
$VERSION = "0.04";
$DATE = '2018-11-01';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        get_rEndWholeAlignJudge
                        get_rEndWholeSuppHaplo
                        isInValidPair
                        dEndSameHapJudge
                        sEndSoloHapJudge
                        sEndInterHapJudge
                        dEndInterHapJudge
                        addHapIDtoReadsOptfd
                     /;

#--- structure of object
# basis, BioFuse::BioInfo::Objects::PairEnd_OB

#--- construction of object ---
sub new{
    my $type = shift;

    my $pe_OB = BioFuse::BioInfo::Objects::PairEnd_OB->new( @_ );

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
                     @{ $pe_OB->get_reads_OB( reads_end => $parm{rEnd} ) }
                );
}

#--- return whole HashRef containing supported haplo-id of given rEnd ---
sub get_rEndWholeSuppHaplo{
    my $pe_OB = shift;
    my %parm = @_;

    my %rWholeHap;
    for my $reads_OB ( @{ $pe_OB->get_reads_OB( reads_end => $parm{rEnd} ) } ){
        my $rAlignHapHref = $reads_OB->get_SuppHaploHref;
        push @{ $rWholeHap{$_} }, @{ $rAlignHapHref->{$_} } for keys %$rAlignHapHref;
    }
    return \%rWholeHap;
}

#--- test is the paired-ends are aligned too close, becoming invalid Hi-C PE ---
sub isInValidPair{
    my $pe_OB = shift;
    my %parm = @_;
    my $minSplitReadGap = $parm{minSplitReadGap};

    my $rOB_sortAref = $pe_OB->get_sorted_reads_OB(rEndAref => [1,2], onlyMap => 1);
    return $rOB_sortAref->[0]->is_closeAlign( test_rOB => $rOB_sortAref->[-1], distance => $minSplitReadGap );
}

#--- judgement on dEnd-hx PE-reads ---
## 1) at least one end must have 'UK' alignment, then go to 2)
## 2) the 'hx' alignment of two ends cannot be close
## return 0 means keep dEnd, 1 means change to sEnd
sub dEndSameHapJudge{
    my $pe_OB = shift;
    my %parm = @_;
    my $minSplitReadGap = $parm{minSplitReadGap};

    my $R1_rOB_Aref = $pe_OB->get_reads_OB( reads_end => 1 );
    my $R2_rOB_Aref = $pe_OB->get_reads_OB( reads_end => 2 );

    # at least one end must have 'UK' alignment
    # else, keep dEnd
    my $R1_hasUK = grep ! $_->has_SuppHaplo, @$R1_rOB_Aref;
    my $R2_hasUK = grep ! $_->has_SuppHaplo, @$R2_rOB_Aref;
    return 0 unless( $R1_hasUK || $R2_hasUK );

    # then, check whether the 'hx' alignment of two ends is close-Align
    for my $r1_OB ( grep $_->has_SuppHaplo, @$R1_rOB_Aref ){ # skip 'UK'
        for my $r2_OB ( grep $_->has_SuppHaplo, @$R2_rOB_Aref ){ # skip 'UK'
            # close !
            if( $r1_OB->is_closeAlign( test_rOB => $r2_OB, distance => $minSplitReadGap ) ){
                return 1;
            }
        }
    }
    # not close!
    # return to keep dEnd
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
    my $minSplitReadGap = $parm{minSplitReadGap};

    my @shEnd_HasHap_rOB = grep $_->has_SuppHaplo, @{$pe_OB->get_reads_OB(reads_end => $shEnd)};
    if(      @shEnd_HasHap_rOB >= 2
        && ! $shEnd_HasHap_rOB[0]->is_closeAlign(test_rOB => $shEnd_HasHap_rOB[-1], distance => $minSplitReadGap)
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
    my $minSplitReadGap = $parm{minSplitReadGap};

    my $mh_rOB_Aref = $pe_OB->get_reads_OB(reads_end => $mhEnd);
    my $uh_rOB_Aref = $pe_OB->get_reads_OB(reads_end => $uhEnd); # much likely to have only one alignment
    my $doSelectBool = 0; # might select
    if(    scalar(@$mh_rOB_Aref) > 1 # have SP
        && $mh_rOB_Aref->[0]->has_SuppHaplo # both alignments have supported haplotype
        && $mh_rOB_Aref->[1]->has_SuppHaplo
    ){
        # pairwise alignments' distance
        my %closeAlign  = ( u1m1 => $uh_rOB_Aref->[0]->is_closeAlign(test_rOB => $mh_rOB_Aref->[0], distance => $minSplitReadGap),
                            m1m2 => $mh_rOB_Aref->[0]->is_closeAlign(test_rOB => $mh_rOB_Aref->[1], distance => $minSplitReadGap),
                            u1m2 => $uh_rOB_Aref->[0]->is_closeAlign(test_rOB => $mh_rOB_Aref->[1], distance => $minSplitReadGap)  );
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
    my $minSplitReadGap = $parm{minSplitReadGap};

    # prepare hapID -> reads_OB
    my %hapID2rOB;
    for my $rEnd (1, 2){
        my $rOB_Aref = $pe_OB->get_reads_OB(reads_end => $rEnd);
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
                unless( $rOB_iA->[$r_i]->is_closeAlign(test_rOB => $rOB_iA->[$r_j], distance => $minSplitReadGap) ){
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
                    if( $rOB_iA->[$r_i]->is_closeAlign(test_rOB => $rOB_jA->[$r_j], distance => $minSplitReadGap) ){
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
    $_->addHapIDtoOptfd for @{ $pe_OB->get_reads_OB( reads_end => $parm{reads_end} ) };
}

1; ## tell the perl script the successful access of this module.
