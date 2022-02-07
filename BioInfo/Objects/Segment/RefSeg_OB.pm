package BioFuse::BioInfo::Objects::Segment::RefSeg_OB;

use strict;
use warnings;
use Data::Dumper;
use List::Util qw/ first any /;
use BioFuse::Util::Log qw/ cluck_and_exit stout_and_sterr /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::Interval qw/ Get_Two_Seg_Olen /;
use BioFuse::BioInfo::FASTA qw/ read_fasta_file /;

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
$VERSION = "0.08";
$DATE = '2022-02-07';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @function_list = qw/
                        new
                        id
                        set_length
                        length
                        foreNlen
                        tailNlen
                        is_inNregion
                        set_note
                        note
                        set_circular
                        is_circular
                        set_extLen
                        extLen
                        has_seqID
                        loadSeqFromFa
                        set_seq
                        seq
                        segSeq
                        init_refpos
                        refpos_mergeExt
                        normInDelMut
                        filter_refpos
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
# refseg -> foreNlen = $foreNlen, this is length of the fore N part in original seq
# refseg -> tailNlen = $foreNlen, this is length of the tail N part in original seq
# refseg -> N_region = [[stp1,edp1],[stp2,edp2],..], this is N part in original seq
# refseg -> note = $note
# refseg -> circular = 0(no)/1(yes)
# refseg -> extLen = $extLen, this is just the extended path length
# refseg -> seq = {orig=>orig_seq, h1=>h1_seq, .., hx=>hx_seq}
# refseg -> regMap = {h1=>[{oSt=>origStp,oEd=>origEdp,hSt=>hapStp,hEd=>hapEdp},..,], h2=>[]}
# refseg -> mut -> seqID -> stp -> {edp=>S, type=>S, ref=>S, allele=>S, discard=>0/1, comment=>S}
#           note: snv/ins/el/[mnp]  ins:+xx; del:-xx|-\d+
# refseg -> refpos -> seqID = { pos => $refpos_OB, ...}

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
    my %parm = @_;
    my $seqID = $parm{seqID};

    if(    defined $seqID
        && defined $refseg->{seq}->{$seqID}
    ){
        return CORE::length($refseg->{seq}->{$seqID});
    }
    else{
        return $refseg->{length};
    }
}

#--- get refseg's foreNlen ---
sub foreNlen{
    my $refseg = shift;
    my %parm = @_;
    my $seqID = $parm{seqID};

    if(    defined $seqID
        && defined $refseg->{seq}->{$seqID}
    ){
        my ($foreNstr) = ($refseg->{seq}->{$seqID} =~ /^(N*)/i);
        return CORE::length($foreNstr);
    }
    else{
        return $refseg->{foreNlen};
    }
}

#--- get refseg's tailNlen ---
sub tailNlen{
    my $refseg = shift;
    my %parm = @_;
    my $seqID = $parm{seqID};

    if(    defined $seqID
        && defined $refseg->{seq}->{$seqID}
    ){
        my ($tailNstr) = ($refseg->{seq}->{$seqID} =~ /(N*)$/i);
        return CORE::length($tailNstr);
    }
    else{
        return $refseg->{tailNlen};
    }
}

#--- check whether given pos is in N-region ---
sub is_inNregion{
    my $refseg = shift;
    my %parm = @_;
    my $pos = $parm{pos};
    my $flk = $parm{flk} || 0;
    return any {$pos>=$_->[0]-$flk && $pos<=$_->[1]+$flk} @{$refseg->{N_region}};
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

#--- load sequence from fasta --- 
sub loadSeqFromFa{
    my $refseg = shift;
    my %parm = @_;
    my $fa = $parm{fa};
    my $segName = $parm{segName} || $refseg->id;
    my $seqID = $parm{seqID} || 'orig';

    # read fasta file
    read_fasta_file( FaFile => $fa,
                     needSeg_Href => {$segName=>1},
                     SrefToSaveSeq => \$refseg->{seq}->{$seqID}
                   );
    # inform
    stout_and_sterr "[INFO]\tload refseg ($segName) seq.\n";
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

#--- init refpos info from seq ---
sub init_refpos{
    my $refseg = shift;
    my %parm = @_;
    my $seqID = $parm{seqID} || 'orig';
    my $refposHf = $parm{refposHf}; # may be other refpos hash
       $refposHf = \%{$refseg->{refpos}->{$seqID}} if !defined $refposHf; # internal refpos hash

    my @seq = split //, $refseg->seq(seqID=>$seqID);
    $refposHf->{$_+1} = BioFuse::BioInfo::Objects::Allele::RefPos_OB->new(pos=>$_+1, refAllele=>$seq[$_]) for 0..$#seq;
    # inform
    stout_and_sterr "[INFO]\tmake refpos objects of seq ($seqID).\n";
}

#--- merge the pos info of the extended part ---
## user could provide other refpos hash belongs to this refpos
## note the 'length' is after extended
sub refpos_mergeExt{
    my $refseg = shift;
    my %parm = @_;
    my $seqID = $parm{seqID} || 'orig';
    my $length = $parm{length} || $refseg->length(seqID=>$seqID);
    my $extLen = $parm{extLen} || $refseg->extLen;
    my $refposHf = $parm{refposHf}; # may be other refpos hash
       $refposHf = $refseg->{refpos}->{$seqID} if !defined $refposHf; # internal refpos hash

    for my $shift (1 .. $extLen){
        my $donor_pos = $shift + $length - $extLen;
        # merge info
        my $accept_refpos_OB = $refposHf->{$shift};
        my $donor_refpos_OB  = $refposHf->{$donor_pos};
        print "$shift\n$donor_pos\n" unless defined $accept_refpos_OB;
        $accept_refpos_OB->merge(donor=>$donor_refpos_OB) if defined $donor_refpos_OB;
        # sweep
        delete $refposHf->{$donor_pos};
    }
}

#--- deal with the floating InDel ---
## similar to vcf norm
sub normInDelMut{
    my $refseg = shift;
    my %parm = @_;
    my $seqID = $parm{seqID} || 'orig';
    my $refposHf = $parm{refposHf}; # may be other refpos hash
       $refposHf = $refseg->{refpos}->{$seqID} if !defined $refposHf; # internal refpos hash

    my $norm_again = 1;
    while($norm_again){
        $norm_again = 0; # set 1 when norm happens
        for my $pos (sort {$a<=>$b} keys %$refposHf){
            my $refpos_OB = $refposHf->{$pos};
            next unless $refpos_OB->has_mutation;
            for my $mut_id (@{$refpos_OB->mutList}){
                my ($mut_type,$mut_seq) = (split /,/,$mut_id)[0,1];
                next if $mut_type eq 'snp'; # only for ins or del
                # find match tail-part between mut_seq and fore_refseq
                my $pos_shift = $mut_type eq 'ins' ? 1 : 0;
                my $match_len = 0;
                for my $i (1 .. CORE::length($mut_seq)){
                    my $mut_base = substr($mut_seq, -1*$i, 1);
                    my $fore_pos = $pos - $i + $pos_shift;
                    my $fore_refpos_OB = $refposHf->{$fore_pos};
                    print "$i\n" if($pos == 32324);
                    print Dumper($fore_refpos_OB) if($pos == 32324);
                    last unless defined $fore_refpos_OB;
                    last if $fore_refpos_OB->has_mutation && $fore_pos != $pos; # donot cross mutation
                    if($fore_refpos_OB->refAllele =~ /^$mut_base$/i){ # match
                        $match_len = $i;
                    }
                    else{
                        last;
                    }
                }
                if($pos == 32324){
                    print "normInDelMut\t$mut_id\t$match_len\n";
                }
                # norm: move mutation
                if($match_len){
                    my $accept_pos = $pos - $match_len;
                    my $accept_refpos_OB = $refposHf->{$accept_pos};
                    my $new_mut_seq = substr($mut_seq, -1*$match_len)
                                     .substr($mut_seq, 0, CORE::length($mut_seq)-$match_len);
                    my $new_mut_id = "$mut_type,$new_mut_seq";
                    # mut depth merge
                    my $mutWhSup = $refpos_OB->mutDepth(mut_id=>$mut_id, type=>'S');
                    my $mutFwSup = $refpos_OB->mutDepth(mut_id=>$mut_id, type=>'F');
                    my $mutRvSup = $refpos_OB->mutDepth(mut_id=>$mut_id, type=>'R');
                    $accept_refpos_OB->addMutDepth(mut_id=>$new_mut_id, add_fw=>$mutFwSup, add_rv=>$mutRvSup, add_cp => -1*$mutWhSup);
                    $accept_refpos_OB->setRefDepth(refdepthAf=>[@{$refpos_OB->RefDepth(type=>'A')}]);
                    # deletion
                    if($mut_type eq 'del'){
                        # decrease accepter's ref allele depth
                        $accept_refpos_OB->addRefDepth(add_fw => -1*$mutFwSup, add_rv => -1*$mutRvSup);
                        # increase deleted part's depth
                        for my $j (1 .. $match_len){
                            my $del_refpos_OB = $refposHf->{$pos+$j-1};
                            next unless defined $del_refpos_OB;
                            $del_refpos_OB->addDepth(add=>$mutWhSup);
                            $del_refpos_OB->addRefDepth(add_fw=>$mutFwSup, add_rv=>$mutRvSup);
                        }
                    }
                    # sweep original mutation
                    $refpos_OB->delete_mut(mut_id=>$mut_id);
                    # mark
                    $norm_again = 1; # check whole refpos and try norm again
                    # inform
                    stout_and_sterr "[INFO]\tnorm mutation: 'pos=$pos;mut_id=$mut_id' to 'pos=$accept_pos;mut_id=$new_mut_id'.\n";
                }
            }
        }
    }
}

#--- filter pos info and get bestMut ---
## note 'find_bestMut' func has deleted the mutation fails basic filtration
sub filter_refpos{
    my $refseg = shift;
    my %parm = @_;
    my $seqID = $parm{seqID} || 'orig';
    my $minDepth = $parm{minDepth} || 10; # min whole depth of the position
    my $minAltRC = $parm{minAltRC} || 3; # min altation-supportting reads count
    my $biStrdRC = $parm{biStrdRC} || 0; # min altation-supportting reads count on both stand
    my $minIDist = $parm{minIDist} || 0; # min distance between adjacent InDels
    my $skipDELp = $parm{skipDELp} || 0; # do not deal with pos in the DEL
    my $method   = $parm{method} || 'Bwa-SAMtools-VCF'; # for tag
    my $onlyMUTp = $parm{onlyMUTp} || 0; # only keep pos with bestMut
    my $refposHf = $parm{refposHf}; # may be other refpos hash
       $refposHf = $refseg->{refpos}->{$seqID} if !defined $refposHf; # internal refpos hash
    my $max_rlen = $parm{max_rlen} || 100; # to judge linear ref bilateral ends
    my $refEndTm = $parm{refEndTm} || 3; # times of min supporting count required for linear ref ends mut or N-region shore
    my $debug    = $parm{debug};

    # refend region for linear ref
    my ($foreRefEnd,$tailRefEnd) =   $refseg->is_circular
                                   ? (0, 0)
                                   : ($refseg->foreNlen+$max_rlen*2, $refseg->length-$refseg->tailNlen-$max_rlen*2);

    # filter mutations
    ## basic criterion; adjacent indel exclusion
    my $last_ID_refpos_OB;
    my $minPos = 1;
    for my $pos (sort {$a<=>$b} keys %$refposHf){
        my $refpos_OB = $refposHf->{$pos};
        # obtain the best mutation
        if(   $refpos_OB->depth < $minDepth   # depth filtration
           || $refpos_OB->refAllele =~ /^N$/i # do not allow mutation at 'N' locus
           || ($skipDELp && $pos < $minPos)   # skip pos when skipDELp
          ){
            # inform
            stout_and_sterr "[INFO]\tdiscard allMut ($method) at pos ".$refpos_OB->pos.", due to basic filtration.\n".Dumper($refpos_OB) if $debug && $refpos_OB->has_mutation;
            # delete it anyway
            $refpos_OB->sweep_mutation;
        }
        elsif( $refpos_OB->has_mutation ){ # has mutation(s), so check whether it is good enough
            # try to find the best qualified mutation
            my $refEnd = (   (!$refseg->is_circular && ($pos<=$foreRefEnd || $pos>=$tailRefEnd)) # linear reference bilateral end
                          || $refseg->is_inNregion(pos=>$pos, flk=>$max_rlen*2) # N-region shore
                         ) ? $refEndTm : 0; # for mut at end of linear ref
            $refpos_OB->find_bestMut(minAltRC=>$minAltRC, biStrdRC=>$biStrdRC, method=>$method, debug=>$debug, refEnd=>$refEnd);
            # no qualified mutation!
            unless($refpos_OB->has_bestMut){
                next;
            }
            # compare adjacent indel supports
            if($refpos_OB->bestMut_type =~ /ins|del/){
                if(   defined $last_ID_refpos_OB
                   && $pos - $last_ID_refpos_OB->pos < $minIDist # too close!
                ){
                    if($last_ID_refpos_OB->bestMut_sumSup <= $refpos_OB->bestMut_sumSup){
                        $last_ID_refpos_OB->remove_bestMut;
                        $last_ID_refpos_OB = $refpos_OB; # update
                        # inform
                        stout_and_sterr "[INFO]\tdiscard bestMut ($method) at pos ".$last_ID_refpos_OB->pos.", due to adjacent indel comparison.\n" if $debug;
                    }
                    else{
                        $refpos_OB->remove_bestMut;
                        # inform
                        stout_and_sterr "[INFO]\tdiscard bestMut ($method) at pos ".$refpos_OB->pos.", due to adjacent indel comparison.\n" if $debug;
                    }
                }
                else{
                    $last_ID_refpos_OB = $refpos_OB; # update
                }
            }
            # skipDELp, update pos filteration
            if($skipDELp && $refpos_OB->bestMut_type eq 'del'){
                $minPos = $pos + CORE::length($refpos_OB->bestMut_seq);
            }
        }
    }
    # sweep pos without bestMut
    if($onlyMUTp){
        delete $refposHf->{$_} for grep !$refposHf->{$_}->has_bestMut, keys %$refposHf;
    }
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
    my $noAlert = $parm{noAlert} || 0;

    # orig?
    return $origPos || $hapPos if $seqID eq 'orig';

    my $regMapAf = $refseg->regMapAf(seqID=>$seqID);
    if(defined $origPos){
        my $rmpHf = first {    $_->{oSt}<=$origPos
                            && $_->{oEd}>=$origPos
                          } @$regMapAf;
        if(!defined $rmpHf){
            if($noAlert){
                return undef;
            }
            else{
                cluck_and_exit "<ERROR>\tcannot map origPos ($origPos) to haplotype.\n";
            }
        }
        return $rmpHf->{hSt} + $origPos - $rmpHf->{oSt};
    }
    elsif(defined $hapPos){
        my $rmpHf = first {    $_->{hSt}<=$hapPos
                            && $_->{hEd}>=$hapPos
                          } @$regMapAf;
        if(!defined $rmpHf){
            if($noAlert){
                return undef;
            }
            else{
                cluck_and_exit "<ERROR>\tcannot map hapPos ($hapPos) to original reference.\n";
            }
        }
        return $rmpHf->{oSt} + $hapPos - $rmpHf->{hSt};
    }
}

1; ## tell the perl script the successful access of this module.
