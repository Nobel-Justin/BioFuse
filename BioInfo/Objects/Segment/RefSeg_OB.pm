package BioFuse::BioInfo::Objects::Segment::RefSeg_OB;

use strict;
use warnings;
use Data::Dumper;
use BioFuse::Util::Log qw/ cluck_and_exit /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;

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
$VERSION = "0.01";
$DATE = '2019-05-24';

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
                     /;

#--- structure of object
# refseg -> id = $id
# refseg -> length = $length
# refseg -> note = $note
# refseg -> tmp = $tmp

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

1; ## tell the perl script the successful access of this module.
