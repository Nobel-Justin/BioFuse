package BioFuse::BioInfo::Objects::Allele::RefPos_OB;

use strict;
use warnings;
use Data::Dumper;
use List::Util qw/ max /;
use BioFuse::Util::Log qw/ cluck_and_exit /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Objects::Allele::RefPos_OB';
#----- version --------
$VERSION = "0.01";
$DATE = '2021-12-29';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        setDepth
                        addDepth
                        setRefDepth
                        addRefDepth
                        RefDepth
                        mutDepth
                        addMut
                        addMutDepth
                        pos
                        refAllele
                        depth
                        merge
                        has_mutation
                        mutList
                        sweep_mutation
                        delete_mut
                        find_bestMut
                        has_bestMut
                        copy_bestMut
                        remove_bestMut
                        reset_bestMut_seq
                        bestMut_id
                        bestMut_infoStr
                        bestMut_type
                        bestMut_seq
                        bestMut_sumSup
                        bestMut_fwSup
                        bestMut_rvSup
                        bestMut_method
                     /;

#--- structure of object
# refpos_OB -> chr = $chr
# refpos_OB -> pos = $pos
# refpos_OB -> refAllele = $refAllele
# refpos_OB -> depth = $depth
# refpos_OB -> refDepth = [$whole, $Forward, $Reverse]
# refpos_OB -> mutation = {mut_id=>[$mut_wholeDepth, $mut_ForwardDepth, $mut_ReverseDepth, $ref_wholeDepth]}
#   mut_id format: 'TYPE,SEQ', e.g., 'snv,A'; 'ins,TC'; 'del,AGC'
# refpos_OB -> bestMut = {sumSup=>S, fwSup=>S, rvSup=>S, mutType=>S, mutSeq=>S, method=>S}

#--- construction of object ---
sub new{
    my $type = shift;
    my %parm = @_;

    my $refpos_OB = {};
    $refpos_OB->{pos} = $parm{pos};
    $refpos_OB->{depth} = 0;
    $refpos_OB->{refAllele} = $parm{refAllele};
    $refpos_OB->{refDepth} = [0,0,0]; # W,F,R

    bless($refpos_OB);
    return $refpos_OB;
}

#--- set total depth ---
sub setDepth{
    my $refpos_OB = shift;
    my %parm = @_;
    $refpos_OB->{depth} = $parm{depth};
}

#--- add to total depth ---
sub addDepth{
    my $refpos_OB = shift;
    my %parm = @_;
    $refpos_OB->{depth} += $parm{add};
}

#--- set array of ref depth ---
sub setRefDepth{
    my $refpos_OB = shift;
    my %parm = @_;
    $refpos_OB->{refDepth} = $parm{refdepthAf};
}

#--- add to ref depth ---
sub addRefDepth{
    my $refpos_OB = shift;
    my %parm = @_;
    my $add_fw = $parm{add_fw} || 0;
    my $add_rv = $parm{add_rv} || 0;

    $refpos_OB->{refDepth}->[0] += ($add_fw + $add_rv);
    $refpos_OB->{refDepth}->[1] +=  $add_fw;
    $refpos_OB->{refDepth}->[2] +=  $add_rv;
    # minimum should be zero
    $_ = max($_,0) for @{$refpos_OB->{refDepth}};
}

#--- return ref allele depth ---
sub RefDepth{
    my $refpos_OB = shift;
    my %parm = @_;
    my $type = $parm{type} || 'S'; # Sum/Fw/Rv
    return $refpos_OB->{refDepth}      if $type eq 'A'; # return depth array-ref
    return $refpos_OB->{refDepth}->[0] if $type eq 'S';
    return $refpos_OB->{refDepth}->[1] if $type eq 'F';
    return $refpos_OB->{refDepth}->[2] if $type eq 'R';
}

#--- return mutation depth array-ref ---
sub mutDepth{
    my $refpos_OB = shift;
    my %parm = @_;
    my $mut_id = $parm{mut_id};
    my $type = $parm{type} || 'S'; # Sum/Fw/Rv/Compare/Array-ref

    if(!$refpos_OB->has_mutation(mut_id=>$mut_id)){
        cluck_and_exit "<ERROR>\tcan not find mut ($mut_id) in refpos object.\n".Dumper($refpos_OB);
    }
    return $refpos_OB->{mutation}->{$mut_id}      if $type eq 'A'; # return depth array-ref
    return $refpos_OB->{mutation}->{$mut_id}->[0] if $type eq 'S';
    return $refpos_OB->{mutation}->{$mut_id}->[1] if $type eq 'F';
    return $refpos_OB->{mutation}->{$mut_id}->[2] if $type eq 'R';
    return $refpos_OB->{mutation}->{$mut_id}->[3] if $type eq 'C'; # ref allele depth to compare
}

#--- record mutation ---
sub addMut{
    my $refpos_OB = shift;
    my %parm = @_;
    my $mut_id = $parm{mut_id};
    my $depthAf = $parm{depthAf};

    if(!exists $refpos_OB->{mutation}->{$mut_id}){
        $refpos_OB->{mutation}->{$mut_id} = $depthAf;
    }
    else{
        $refpos_OB->{mutation}->{$mut_id}->[$_] += $depthAf->[$_] for 0 .. @$depthAf-1;
    }
}

#--- add depth of mutation ---
sub addMutDepth{
    my $refpos_OB = shift;
    my %parm = @_;
    my $mut_id = $parm{mut_id};
    my $add_fw = $parm{add_fw} || 0;
    my $add_rv = $parm{add_rv} || 0;
    my $add_cp = $parm{add_cp} || 0; # ref allele depth to compare

    $refpos_OB->{mutation}->{$mut_id}->[0] += ($add_fw + $add_rv);
    $refpos_OB->{mutation}->{$mut_id}->[1] +=  $add_fw;
    $refpos_OB->{mutation}->{$mut_id}->[2] +=  $add_rv;
    $refpos_OB->{mutation}->{$mut_id}->[3] +=  $add_cp;
    # minimum should be zero
    $_ = max($_,0) for @{$refpos_OB->{mutation}->{$mut_id}};
}

#--- return pos ---
sub pos{
    my $refpos_OB = shift;
    return $refpos_OB->{pos};
}

#--- return refAllele ---
sub refAllele{
    my $refpos_OB = shift;
    return $refpos_OB->{refAllele};
}

#--- return total depth ---
sub depth{
    my $refpos_OB = shift;
    return $refpos_OB->{depth};
}

#--- merge with given refpos_OB ---
sub merge{
    my $refpos_OB = shift;
    my %parm = @_;
    my $donor_OB = $parm{donor}; # also is refpos OB

    # check
    if($refpos_OB->refAllele ne $donor_OB->refAllele){
        cluck_and_exit "<ERROR>\tdonor refpos_OB refAllele not match to acceptor refpos_OB\n"
                             ."\tdonor refpos_OB:\n".Dumper($donor_OB)
                             ."\tacceptor refpos_OB:\n".Dumper($refpos_OB);
    }
    # depth
    $refpos_OB->addDepth(add=>$donor_OB->depth);
    # refDepth
    $refpos_OB->addRefDepth(add_fw=>$donor_OB->RefDepth(type=>'F'));
    $refpos_OB->addRefDepth(add_rv=>$donor_OB->RefDepth(type=>'R'));
    # mutation
    $refpos_OB->addMut(mut_id=>$_, depthAf=>[@{$donor_OB->mutDepth(mut_id=>$_,type=>'A')}]) for @{$donor_OB->mutList};
}

#--- has mutations ---
sub has_mutation{
    my $refpos_OB = shift;
    my %parm = @_;
    my $mut_id = $parm{mut_id};

    if(    exists $refpos_OB->{mutation}
        && scalar(keys %{$refpos_OB->{mutation}}) != 0
        && (   !defined $mut_id
            || exists $refpos_OB->{mutation}->{$mut_id}
           )
    ){
        return 1;
    }
    else{
        $refpos_OB->sweep_mutation; # try anyway
        return 0;
    }
}

#--- return mutation list ---
sub mutList{
    my $refpos_OB = shift;
    if(exists $refpos_OB->{mutation}){
        return [ sort keys %{$refpos_OB->{mutation}} ];
    }
    else{
        return [];
    }
}

#--- sweep all mutations ---
sub sweep_mutation{
    my $refpos_OB = shift;
    delete $refpos_OB->{mutation};
    delete $refpos_OB->{bestMut};
}

#--- delete given mut_id ---
sub delete_mut{
    my $refpos_OB = shift;
    my %parm = @_;
    my $mut_id = $parm{mut_id};
    delete $refpos_OB->{mutation}->{$mut_id};
    $refpos_OB->sweep_mutation if scalar(keys %{$refpos_OB->{mutation}})==0;
}

#--- find best mutation ---
sub find_bestMut{
    my $refpos_OB = shift;
    my %parm = @_;
    my $minAltRc = $parm{minAltRc} || 3; # min altation-supportting reads count
    my $biStrdRC = $parm{biStrdRC} || 0; # min altation-supportting reads count on both stand
    my $method   = $parm{method} || 'Bwa-SAMtools-VCF'; # for tag

    for my $mut_id (sort keys %{$refpos_OB->{mutation}}){
        my ($mut_type, $mut_seq) = split /,/, $mut_id;
        my ($sumSup, $fwSup, $rvSup, $refSup_toCmp) = @{$refpos_OB->{mutation}->{$mut_id}};
        # filter the mutation
        if(    ( $sumSup < $minAltRc ) # at least so many supported reads
               # if required, both strand supported
            || ( $biStrdRC && ($fwSup < $biStrdRC || $rvSup < $biStrdRC) )
        ){
            $refpos_OB->delete_mut(mut_id=>$mut_id);
            next;
        }
        elsif( # must be dominant against ref-allel
               ( $refSup_toCmp >  0 && $sumSup <= $refSup_toCmp )
            || ( $refSup_toCmp <= 0 && $sumSup <= $refpos_OB->RefDepth(type=>'S') )
                 # must be dominant against ref-allel, for non-del
            #    ( $mut_type !~ /del/ && $sumSup <= $refSup_toCmp )
            #    # must be dominant against ref-allel, for del specifically
            # || ( $mut_type =~ /del/ && $sumSup <= $refpos_OB->RefDepth(type=>'S') )
          ){
            next;
        }
        elsif( !exists $refpos_OB->{bestMut}
            || $sumSup > $refpos_OB->{bestMut}->{sumSup}
        ){  # update
            $refpos_OB->{bestMut} = { sumSup=>$sumSup, fwSup=>$fwSup, rvSup=>$rvSup, mutType=>$mut_type, mutSeq=>$mut_seq, method=>$method };
        }
    }
}

#--- has best mutation ---
sub has_bestMut{
    my $refpos_OB = shift;
    return exists $refpos_OB->{bestMut} ? 1 : 0;
}

#--- set best mutation --
sub copy_bestMut{
    my $refpos_OB = shift;
    my %parm = @_;
    my $donor_OB = $parm{donor}; # another refpos_OB with bestMut

    # check
    if($refpos_OB->has_bestMut){
        cluck_and_exit "<ERROR>\tcannot set bestMut for refpos_OB, as it already has.\n".Dumper($refpos_OB);
    }
    if(!$donor_OB->has_bestMut){
        cluck_and_exit "<ERROR>\tdonor refpos_OB lacks bestMut.\n".Dumper($donor_OB);
    }
    # copy
    $refpos_OB->{bestMut} = $donor_OB->{bestMut};
}

#--- remove the best mutation ---
sub remove_bestMut{
    my $refpos_OB = shift;
    delete $refpos_OB->{bestMut};
}

#--- reset seq of the best mutation --
sub reset_bestMut_seq{
    my $refpos_OB = shift;
    my %parm = @_;
    $refpos_OB->{bestMut}->{mutSeq} = $parm{mut_seq};
}

#--- return type of the best mutation ---
sub bestMut_id{
    my $refpos_OB = shift;
    return join(',', map {("$refpos_OB->{bestMut}->{$_}")} qw/ mutType mutSeq/);
}

#--- return type of the best mutation ---
sub bestMut_infoStr{
    my $refpos_OB = shift;
    return join(';', map {("$_=$refpos_OB->{bestMut}->{$_}")} qw/ mutType mutSeq sumSup fwSup rvSup method/);
}

#--- return type of the best mutation ---
sub bestMut_type{
    my $refpos_OB = shift;
    return $refpos_OB->{bestMut}->{mutType};
}

#--- return seq of the best mutation ---
sub bestMut_seq{
    my $refpos_OB = shift;
    return $refpos_OB->{bestMut}->{mutSeq};
}

#--- return sumSup of the best mutation ---
sub bestMut_sumSup{
    my $refpos_OB = shift;
    return $refpos_OB->{bestMut}->{sumSup};
}

#--- return fwSup of the best mutation ---
sub bestMut_fwSup{
    my $refpos_OB = shift;
    return $refpos_OB->{bestMut}->{fwSup};
}

#--- return rvSup of the best mutation ---
sub bestMut_rvSup{
    my $refpos_OB = shift;
    return $refpos_OB->{bestMut}->{rvSup};
}

#--- return rvSup of the best mutation ---
sub bestMut_method{
    my $refpos_OB = shift;
    return $refpos_OB->{bestMut}->{method};
}

1; ## tell the perl script the successful access of this module.
