package BioFuse::BioInfo::Objects::Sample::Sample_OB;

use strict;
use warnings;
use Data::Dumper;
use BioFuse::Util::Log qw/ cluck_and_exit /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Objects::Sample::Sample_OB';
#----- version --------
$VERSION = "0.01";
$DATE = '2019-05-23';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        setType
                        type
                        setPatient
                        patient
                        setTissue
                        tissue
                        set_note
                        note
                        addBam
                        bam
                        addFQ
                     /;

#--- structure of object
# sample -> id = $id
# sample -> type = 'case'/'control'
# sample -> tissue = $tissue
# sample -> patient = $patient_object
# sample -> note = $note
# sample -> bam = [$bam, $other_bam, ...]
# sample -> fq = [[$fq_PE1,$fq_PE2], [$fq_SE], ...]

#--- construction of object
sub new{
    my $type = shift;
    my %parm = @_;

    my $sample = {};
    $sample->{id} = $parm{id};

    bless($sample);
    return $sample;
}

#--- set sample type ---
## case OR control
sub setType{
    my $sample = shift;
    my %parm = @_;
    $sample->{type} = $parm{type};
}

#--- return sample type ---
## case OR control
sub type{
    my $sample = shift;
    return $sample->{type};
}

#--- set sample's patient ---
sub setPatient{
    my $sample = shift;
    my %parm = @_;
    my $patient = $parm{patient_OB};
    $sample->{patient} = $patient;
    $patient->addSample(sample=>$sample);
}

#--- return sampl'se patient object ---
sub patient{
    my $sample = shift;
    return $sample->{patient};
}

#--- set sample tissue source ---
sub setTissue{
    my $sample = shift;
    my %parm = @_;
    $sample->{tissue} = $parm{tissue};
}

#--- return sample tissue source ---
sub tissue{
    my $sample = shift;
    return $sample->{tissue};
}

#--- set sample's notes ---
sub set_note{
    my $sample = shift;
    my %parm = @_;
    $sample->{note} = $parm{note};
}

#--- return sample's notes ---
sub note{
    my $sample = shift;
    return $sample->{note};
}

#--- add sample's bam ---
sub addBam{
    my $sample = shift;
    my %parm = @_;
    push @{ $sample->{bam} }, $parm{bam};
}

#--- return Array-ref of sample's bam ---
sub bam{
    my $sample = shift;
    return $sample->{bam};
}

#--- add sample's fq ---
sub addFQ{
    my $sample = shift;
    my %parm = @_;
    push @{ $sample->{fq} }, $parm{fq_Af};
}

1; ## tell the perl script the successful access of this module.
