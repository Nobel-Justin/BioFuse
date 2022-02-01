package BioFuse::BioInfo::Objects::SeqData::HicReads_OB;

use BioFuse::BioInfo::Objects::SeqData::Reads_OB; # inheritance

use strict;
use warnings;
require Exporter;
use Data::Dumper;
use BioFuse::Util::Log qw/ warn_and_exit stout_and_sterr /;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter BioFuse::BioInfo::Objects::SeqData::Reads_OB);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Objects::SeqData::HicReads_OB';
#----- version --------
$VERSION = "0.03";
$DATE = '2018-11-04';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @function_list = qw/
                        new
                        load_AlignJudge
                        load_SuppHaplo
                        del_SuppHaplo
                        onlyKeep_SuppHaplo
                        get_AlignJudge
                        get_SuppHaploHref
                        get_SuppHaploStr
                        has_SuppHaplo
                        is_fromUnPhasedRegRand
                        addHapIDtoOptfd
                        recover_SuppHaploAttr
                     /;

#--- structure of object
# basis, BioFuse::BioInfo::Objects::SeqData::Reads_OB
# reads_OB -> hapID = {$hapID=>[$allele_OB, ..], ..}
# reads_OB -> alignJudge = $alignJudgeString

#--- construction of object ---
sub new{
    my $type = shift;

    my $reads_OB = BioFuse::BioInfo::Objects::SeqData::Reads_OB->new( @_ );

    bless($reads_OB);
    return $reads_OB;
}

#--- record alignment judgement string ---
sub load_AlignJudge{
    my $reads_OB = shift;
    my %parm = @_;
    my $substitute = $parm{substitute} || 0;

    if( $substitute ){
        $reads_OB->{alignJudge}  = $parm{alignJudge};
    }
    else{
        $reads_OB->{alignJudge} .= $parm{alignJudge};
    }
}

#--- record haplo-id this reads supports ---
# with allele_OB supports this haplo
sub load_SuppHaplo{
    my $reads_OB = shift;
    my %parm = @_;
    push @{ $reads_OB->{hapID}->{ $parm{hapID} } }, $parm{allele_OB};
}

#--- delete haplo-id this reads supports ---
sub del_SuppHaplo{
    my $reads_OB = shift;
    my %parm = @_;
    delete $reads_OB->{hapID}->{ $parm{hapID} };
}

#--- delete all other haplo-id(s) this reads supports ---
sub onlyKeep_SuppHaplo{
    my $reads_OB = shift;
    my %parm = @_;
    delete $reads_OB->{hapID}->{$_} for grep $_ ne $parm{hapID}, keys %{ $reads_OB->get_SuppHaploHref };
}

#--- return alignment judgement string ---
sub get_AlignJudge{
    my $reads_OB = shift;
    my %parm = @_;
    return $reads_OB->{alignJudge} || '';
}

#--- return HashRef containing supported haplo-id ---
sub get_SuppHaploHref{
    my $reads_OB = shift;
    return $reads_OB->{hapID} || {};
}

#--- return HashRef containing supported haplo-id ---
sub get_SuppHaploStr{
    my $reads_OB = shift;
    if( $reads_OB->has_SuppHaplo ){
        return join(',', sort keys %{$reads_OB->get_SuppHaploHref});
    }
    else{
        return 'UK';
    }
}

#--- test whether has haplotype supported ---
sub has_SuppHaplo{
    my $reads_OB = shift;
    return scalar( keys %{ $reads_OB->get_SuppHaploHref } );
}

#--- whether the hapID is set (i.e., h[x]Intra) as flanking region not phased ---
## it and its paired end do not have phased contacts in flanking region
## 'RegionNotPhased'
##  see func 'get_rOBpair_HapLinkCount' in HaploHiC::PhasedHiC::phasedPEtoContact
## 'XU:Z:'
##  see func 'assign_sEndUKend_haplotype' in HaploHiC::PhasedHiC::sEndSoloHapConfirm
##  and func 'assign_dEndUKend_haplotype' in HaploHiC::PhasedHiC::dEndUkHapConfirm
sub is_fromUnPhasedRegRand{
    my $reads_OB = shift;
    return $reads_OB->optfd_has_regex(regex => 'XU:Z:RegionNotPhased');
}

#--- add OWN 'XH:Z:' haplo-id to optfd ---
## if already has, then update
sub addHapIDtoOptfd{
    my $reads_OB = shift;
    my $haploStr = $reads_OB->get_SuppHaploStr;
    if( $reads_OB->{optfd} =~ /XH:Z:\S+/ ){
        $reads_OB->{optfd} =~ s/XH:Z:\S+/XH:Z:$haploStr/;
    }
    else{
        $reads_OB->{optfd} .= "\tXH:Z:$haploStr";
    }
}

#--- recover supported haplo-id from 'XH:Z:' tag ---
sub recover_SuppHaploAttr{
    my $reads_OB = shift;
    if( $reads_OB->{optfd} =~ /XH:Z:(\S+)/ ){
        $reads_OB->load_SuppHaplo( hapID => $_, allele_OB => 'NULL' ) for grep !/UK/i, split /,/, $1;
    }
}

1; ## tell the perl script the successful access of this module.
