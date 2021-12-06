package BioFuse::BioInfo::Objects::Segment::RefSeg_OB;

use strict;
use warnings;
use Data::Dumper;
use List::Util qw/ first /;
use BioFuse::Util::Log qw/ cluck_and_exit stout_and_sterr /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::Interval qw/ Get_Two_Seg_Olen /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Objects::Segment::RefSeg_OB';
#----- version --------
$VERSION = "0.05";
$DATE = '2021-07-22';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        id
                        set_length
                        length
                        set_note
                        note
                        set_circular
                        is_circular
                        set_extLen
                        extLen
                        has_seqID
                        set_seq
                        seq
                        segSeq
                        addMut
                        mutSeqHf
                        mutToSeq
                        writeMutInfo
                        testMutOvp
                        regMapAf
                        writeRegMapInfo
                        regMapPos
                     /;

#--- structure of object
# refseg -> id = $id
# refseg -> length = $length, this is the original length
# refseg -> note = $note
# virus_OB -> circular = 0(no)/1(yes)
# virus_OB -> extLen = $extLen, this is just the extended path length
# refseg -> seq = {orig=>orig_seq, h1=>h1_seq, .., hx=>hx_seq}
# refseg -> regMap = {h1=>[{oSt=>origStp,oEd=>origEdp,hSt=>hapStp,hEd=>hapEdp},..,], h2=>[]}
# refseg -> mut -> seqID -> stp -> {edp=>S, type=>S, ref=>S, allele=>S, discard=>0/1, comment=>S}
#           note: snv/ins/el/[mnp]  ins:+xx; del:-xx|-\d+

#--- construction of object
sub new{
    my $type = shift;
    my %parm = @_;

    my $refseg = {};
    $refseg->{id} = $parm{id};

    bless($refseg);
    return $refseg;
}

#--- get refseg's id ---
sub id{
    my $refseg = shift;
    return $refseg->{id};
}

#--- set refseg's length ---
sub set_length{
    my $refseg = shift;
    my %parm = @_;
    $refseg->{length} = $parm{length};
}

#--- get refseg's length ---
sub length{
    my $refseg = shift;
    return $refseg->{length};
}

#--- set refseg's notes ---
sub set_note{
    my $refseg = shift;
    my %parm = @_;
    $refseg->{note} = $parm{note};
}

#--- return refseg's notes ---
sub note{
    my $refseg = shift;
    return $refseg->{note} || 'N/A';
}

#--- set circular attribute ---
# use 1/0 currently
sub set_circular{
    my $virus_OB = shift;
    my %parm = @_;
    if($parm{circular} eq 'yes' || $parm{circular} eq '1'){
        $virus_OB->{circular} = 1;
    }
    elsif($parm{circular} eq 'no' || $parm{circular} eq '0'){
        $virus_OB->{circular} = 0;
    }
    else{
        cluck_and_exit "<ERROR>\tcannot recognize circular para $parm{circular}.\n";
    }
}

#--- check whether circular or not ---
sub is_circular{
    my $virus_OB = shift;
    return $virus_OB->{circular} || 0; # default 0
}

#--- set extended length specific for circular ---
sub set_extLen{
    my $virus_OB = shift;
    my %parm = @_;
    $virus_OB->{extLen}  = $parm{extLen};
}

#--- return extended length ---
sub extLen{
    my $virus_OB = shift;
    return $virus_OB->{extLen};
}

#--- check whether has given seqID ---
sub has_seqID{
    my $refseg = shift;
    my %parm = @_;
    my $seqID = $parm{seqID};
    my $attr = $parm{attr} || 'mut'; # mut/seq/regMap
    return (exists $refseg->{$attr}->{$seqID} ? 1 : 0);
}

#--- seq sequence ---
sub set_seq{
    my $refseg = shift;
    my %parm = @_;
    my $seqID = $parm{seqID} || 'orig';
    unless(defined $parm{seq}){
        cluck_and_exit "<ERROR>\trequires seq.\n";
    }
    $refseg->{seq}->{$seqID} = $parm{seq};
}

#--- return sequence ---
sub seq{
    my $refseg = shift;
    my %parm = @_;
    my $seqID = $parm{seqID} || 'orig';
    if(   !$parm{test}
       && !exists $refseg->{seq}->{$seqID}
      ){
        cluck_and_exit "<ERROR>\tcannot find $seqID seq.\n";
    }
    return $refseg->{seq}->{$seqID} || undef;
}

#--- get segment sequence based on pos/prime/strd ---
sub segSeq{
    my $refseg = shift;
    my %parm = @_;
    my $seqID = $parm{seqID} || 'orig';
    my $pos = $parm{pos};
    my $strd = $parm{strd}; # +,-
    my $prime = $parm{prime}; # 5p, 3p
    my $span = $parm{span} || 200;
    my $forceLinear = $parm{forceLinear} || 0;

    my $seq = $refseg->seq(seqID=>$seqID);
    my $seqLen = CORE::length($seq);
    # check span
    if($span > $seqLen){
        cluck_and_exit "<ERROR>\tspan ($span) is longer than reference length (".$seqLen.").\n";
    }

    my $segSeq = '';
    if(    ($strd eq '+' && $prime eq '5p')
        || ($strd eq '-' && $prime eq '3p')
    ){ # left hand
        if($pos < $span){ # out of seq range
            if($refseg->is_circular && !$forceLinear){ # circular
                $segSeq = substr($seq,-1*($span-$pos)).substr($seq,0,$pos);
            }
            else{
                $segSeq = substr($seq,0,$pos);
                stout_and_sterr "<WARN>\tget short sequence from $refseg->{id} $seqID genome.\n"
                                     ."\tpos:$pos, span:$span, strd:$strd, prime:$prime, forceLinear:$forceLinear\n";
            }
        }
        else{
            $segSeq = substr($seq,$pos-$span,$span);
        }
    }
    else{ # right hand
        if($pos+$span-1 > $seqLen){ # out of seq range
            if($refseg->is_circular && !$forceLinear){ # circular
                $segSeq = substr($seq,$pos-1).substr($seq,0,$pos+$span-1-$seqLen);
            }
            else{
                $segSeq = substr($seq,$pos-1);
                stout_and_sterr "<WARN>\tget short sequence from $refseg->{id} $seqID genome.\n"
                                     ."\tpos:$pos, span:$span, strd:$strd, prime:$prime, forceLinear:$forceLinear\n";
            }
        }
        else{
            $segSeq = substr($seq,$pos-1,$span);
        }
    }
    # reverse complementary
    ($segSeq = reverse $segSeq) =~ tr/ACGTacgt/TGCAtgca/ if $strd eq '-';

    return $segSeq;
}

#--- record given mut ---
sub addMut{
    my $refseg = shift;
    my %parm = @_;
    $refseg->{mut}->{$parm{seqID}}->{$parm{stp}} = { edp=>$parm{edp},
                                                     type=>$parm{type}, allele=>$parm{allele}, len=>$parm{len},
                                                     discard=>0, comment=>$parm{comment}||'N/A'
                                                   };
}

#--- return mutHf of given seqID ---
sub mutSeqHf{
    my $refseg = shift;
    my %parm = @_;
    return \%{$refseg->{mut}->{$parm{seqID}}};
}

#--- generate seq and region mapping based on mut of given seqID ---
# note: 1) after this, the del allele becomes the deleted sequence, not deleted len any more.
#       2) mut has the 'ref' attribute
#       3) record regmap and seq to corresponding attributes
sub mutToSeq{
    my $refseg = shift;
    my %parm = @_;
    my $seqID = $parm{seqID};
    my $idFsame = $parm{idFsame} || 0;

    my $refID = $refseg->id;
    my $oLen = $refseg->length;
    my $oSeq = uc $refseg->seq(seqID=>'orig');
    # region mapping
    my $regMapAf = $refseg->regMapAf(seqID=>$seqID);
    # each mut
    my $mutSeqHf = $refseg->mutSeqHf(seqID=>$seqID);
    my $lastMutPos = 0;
    my $seq = '';
    for my $mutPos (sort {$a<=>$b} keys %$mutSeqHf){
        my $mutHf = $mutSeqHf->{$mutPos};
        # check pre-dispositioned mut pos: mut in the deleted seq
        if($mutPos <= $lastMutPos){
            $mutHf->{discard} = 1; # mark
            stout_and_sterr "<WARN>\t$refID $seqID LastMutPos:$lastMutPos MutPos:$mutPos encouters overlapped mutation. skip it.\n".Dumper($mutHf);
            next;
        }
        my $type = $mutHf->{type};
        my $mutLen = $mutHf->{len};
        # start mapping
        # note: when one region starts, the xEd pos has offset (-1) to pretend been included in the seq
        if(@$regMapAf == 0){
            if($type eq 'del' && $mutPos == 1){
                $regMapAf->[0] = { oSt=>$mutPos+$mutLen,
                                   oEd=>$mutPos+$mutLen-1,
                                   hSt=>1,
                                   hEd=>0
                                 };
            }
            else{
                $regMapAf->[0] = { oSt=>1,
                                   oEd=>0,
                                   hSt=>1,
                                   hEd=>0
                                 };
            }
        }
        # meet mut
        my $last_oEd = $regMapAf->[-1]->{oEd};
        if($type eq 'snv'){ # not break mapping
            # seq grows
            $seq .= substr($oSeq, $last_oEd, $mutPos-$last_oEd-1) . lc($mutHf->{allele});
            # update reg mapping
            my $regSpan = $mutPos - $last_oEd;
            $regMapAf->[-1]->{oEd} += $regSpan;
            $regMapAf->[-1]->{hEd} += $regSpan;
            # update mut pos boudary
            $lastMutPos = $mutPos;
            # ref attribute
            $mutHf->{ref} = substr($oSeq, $mutPos-1, 1);
        }
        elsif($type eq 'ins'){
            (my $insStr = $mutHf->{allele}) =~ s/[^ACGT]//gi;
            # avoid ins shift due to flanking identical bases
            if(   ($mutPos<$oLen && uc(substr($oSeq,$mutPos,1))   eq uc(substr($insStr, 0,1))) # 3-prime, assuming the alignment is 5-prime preferred, basically
               || ($idFsame      && uc(substr($oSeq,$mutPos-1,1)) eq uc(substr($insStr,-1,1))) # 5-prime, if allowed
            ){
                $mutHf->{discard} = 1; # mark
                stout_and_sterr "<WARN>\t$refID $seqID MutPos:$mutPos encouters shifting INS. skip it.\n".Dumper($mutHf);
                next;
            }
            # seq grows
            $seq .= substr($oSeq, $last_oEd, $mutPos-$last_oEd) . lc($insStr);
            # update reg mapping
            ## last reg finishes
            my $regSpan = $mutPos - $last_oEd;
            $regMapAf->[-1]->{oEd} += $regSpan;
            $regMapAf->[-1]->{hEd} += $regSpan;
            ## new reg starts
            push @$regMapAf, { oSt=>$regMapAf->[-1]->{oEd}+1,
                               oEd=>$regMapAf->[-1]->{oEd},
                               hSt=>$regMapAf->[-1]->{hEd}+$mutLen+1,
                               hEd=>$regMapAf->[-1]->{hEd}+$mutLen
                             };
            # update mut pos boudary
            $lastMutPos = $mutPos;
            # ref attribute
            $mutHf->{ref} = substr($oSeq, $mutPos-1, 1);
        }
        elsif($type eq 'del'){
            # avoid del shift due to flanking identical bases
            if(   ($mutPos+$mutLen<=$oLen && uc(substr($oSeq,$mutPos+$mutLen-1,1)) eq uc(substr($oSeq,$mutPos-1,1))) # 3-prime, assuming the alignment is 5-prime preferred, basically
               || ($idFsame               && $mutPos != 1 && uc(substr($oSeq,$mutPos+$mutLen-2,1)) eq uc(substr($oSeq,$mutPos-2,1))) # 5-prime, if allowed
            ){
                $mutHf->{discard} = 1; # mark
                stout_and_sterr "<WARN>\t$refID $seqID MutPos:$mutPos encouters shifting DEL. skip it.\n".Dumper($mutHf);
                next;
            }
            # seq grows
            $seq .= substr($oSeq, $last_oEd, $mutPos-$last_oEd-1);
            # update reg mapping
            ## last reg finishes
            my $regSpan = $mutPos - $last_oEd;
            $regMapAf->[-1]->{oEd} += $regSpan - 1;
            $regMapAf->[-1]->{hEd} += $regSpan - 1;
            ## new reg starts
            push @$regMapAf, { oSt=>$regMapAf->[-1]->{oEd}+$mutLen+1,
                               oEd=>$regMapAf->[-1]->{oEd}+$mutLen,
                               hSt=>$regMapAf->[-1]->{hEd}+1,
                               hEd=>$regMapAf->[-1]->{hEd}
                             };
            # update mut pos boudary
            $lastMutPos = $mutPos + $mutLen - 1;
            # ref attribute
            $mutHf->{ref} = substr($oSeq, $mutPos-1, 1);
            # update DEL allele, has '-' prefix
            $mutHf->{allele} = '-'.substr($oSeq, $mutPos-1, $mutLen);
        }
    }
    # to the end of the ref genome
    if($lastMutPos < $oLen){
        my $last_oEd = $regMapAf->[-1]->{oEd};
        my $regSpan = $oLen - $lastMutPos;
        # seq grows
        $seq .= substr($oSeq, $last_oEd);
        # finish the last reg
        $regMapAf->[-1]->{oEd} += $regSpan;
        $regMapAf->[-1]->{hEd} += $regSpan;
    }
    # record seq
    $refseg->set_seq(seqID=>$seqID, seq=>$seq);
}

#--- output mut info ---
# refID seqID pos type refBase mutAllele
sub writeMutInfo{
    my $refseg = shift;
    my %parm = @_;
    my $seqID = $parm{seqID};
    my $fileFH = $parm{fileFH} || undef;

    my $refID = $refseg->id;
    my $mutSeqHf = $refseg->mutSeqHf(seqID=>$seqID);
    for my $mutPos (sort {$a<=>$b} keys %$mutSeqHf){
        my $mutHf = $mutSeqHf->{$mutPos};
        my $info = join("\t", $refID, $seqID, $mutPos, uc $mutHf->{type}, $mutHf->{ref}, $mutHf->{allele});
        if(defined $fileFH){
            print {$fileFH} $info."\n";
        }
        else{
            print $info."\n";
        }
    }
}

#--- test overlap of given pos+dist and mut ---
sub testMutOvp{
    my $refseg = shift;
    my %parm = @_;
    my $seqID = $parm{seqID};
    my $pos = $parm{pos};
    my $dist = $parm{dist} || 5;

    my $mutSeqHf = $refseg->mutSeqHf(seqID=>$seqID);
    my $ovp_mutPos = first { $_ <= $pos+$dist  && $mutSeqHf->{$_}->{edp} >= $pos-$dist
                           } sort {$a<=>$b} keys %$mutSeqHf;
    return $mutSeqHf->{$ovp_mutPos} if defined $ovp_mutPos;
    # not found
    return undef;
}

#--- return regMapAf of given seqID ---
sub regMapAf{
    my $refseg = shift;
    my %parm = @_;
    return \@{$refseg->{regMap}->{$parm{seqID}}};
}

#--- output regMap info ---
# refID seqID orig_stp orig_edp hap_stp hap_edp
sub writeRegMapInfo{
    my $refseg = shift;
    my %parm = @_;
    my $seqID = $parm{seqID};
    my $fileFH = $parm{fileFH} || undef;

    my $refID = $refseg->id;
    my $regMapAf = $refseg->regMapAf(seqID=>$seqID);
    for my $rmpHf (@$regMapAf){
        my $info = join("\t", $refID, $seqID, $rmpHf->{oSt}, $rmpHf->{oEd}, $rmpHf->{hSt}, $rmpHf->{hEd});
        if(defined $fileFH){
            print {$fileFH} $info."\n";
        }
        else{
            print $info."\n";
        }
    }
}

#--- map given pos to new pos (two-way) ---
sub regMapPos{
    my $refseg = shift;
    my %parm = @_;
    my $seqID = $parm{seqID};
    my $origPos = $parm{origPos} || undef;
    my $hapPos  = $parm{hapPos}  || undef;

    my $regMapAf = $refseg->regMapAf(seqID=>$seqID);
    if(defined $origPos){
        my $rmpHf = first {    $_->{oSt}<=$origPos
                            && $_->{oEd}>=$origPos
                          } @$regMapAf;
        cluck_and_exit "<ERROR>\tcannot map origPos ($origPos) to haplotype.\n" unless(defined $rmpHf);
        return $rmpHf->{hSt} + $origPos - $rmpHf->{oSt};
    }
    elsif(defined $hapPos){
        my $rmpHf = first {    $_->{hSt}<=$hapPos
                            && $_->{hEd}>=$hapPos
                          } @$regMapAf;
        cluck_and_exit "<ERROR>\tcannot map hapPos ($hapPos) to haplotype.\n" unless(defined $rmpHf);
        return $rmpHf->{oSt} + $hapPos - $rmpHf->{hSt};
    }
}

1; ## tell the perl script the successful access of this module.
