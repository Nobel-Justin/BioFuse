package BioFuse::BioInfo::Objects::Sample::Patient_OB;

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

$MODULE_NAME = 'BioFuse::BioInfo::Objects::Sample::Patient_OB';
#----- version --------
$VERSION = "0.01";
$DATE = '2019-05-23';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        addNote
                        addSample
                     /;

#--- structure of object
# patient -> id = $id
# patient -> sample = {type=>[$sample1,$sample2]}
# patient -> note = $note

#--- construction of object
sub new{
    my $type = shift;
    my %parm = @_;

    my $patient = {};
    $patient->{id} = $parm{id};

    bless($patient);
    return $patient;
}

#--- set patient's notes ---
sub addNote{
    my $patient = shift;
    my %parm = @_;
    $patient->{note} = $parm{note};
}

#--- add sample ---
sub addSample{
    my $patient = shift;
    my %parm = @_;
    my $sample = $parm{sample_OB};
    my $type = $sample->getType;
    unless(defined $type){
        cluck_and_exit "ERROR\tcannot get type of sample\n".Dumper($sample);
    }
    push @{ $patient->{sample}->{$type} }, $sample;
}

1; ## tell the perl script the successful access of this module.
