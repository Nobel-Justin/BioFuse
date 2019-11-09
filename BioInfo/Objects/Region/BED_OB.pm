package BioFuse::BioInfo::Objects::Region::BED_OB;

use strict;
use warnings;
use Data::Dumper;
use List::Util qw/ sum /;
use BioFuse::Util::Log qw/ cluck_and_exit /;
use BioFuse::BioInfo::BED qw/ read_bed_file /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Objects::Region::BED_OB';
#----- version --------
$VERSION = "0.01";
$DATE = '2019-11-10';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        tag
                        length
                        calc_length
                        set_note
                        note
                     /;

#--- structure of object
# bed -> tag = $tag
# bed -> filepath = $filepath
# bed -> length = $length
# bed -> note = $note

#--- construction of object
sub new{
    my $type = shift;
    my %parm = @_;

    my $bed = {};
    $bed->{filepath} = $parm{filepath} || undef;
    $bed->{tag} = $parm{tag} || undef;

    bless($bed);
    return $bed;
}

#--- return file path ---
sub filepath{
    my $bed = shift;
    return $bed->{filepath};
}

#--- return bed tag ---
sub tag{
    my $bed = shift;
    return $bed->{tag};
}

#--- get bed's length ---
sub length{
    my $bed = shift;
    $bed->calc_length unless defined $bed->{length};
    return $bed->{length};
}

#--- calc length ---
sub calc_length{
    my $bed = shift;
    my $ItvHref = read_bed_file(bedFile=>$bed->{filepath}, nonName=>1);
    $bed->{length} = sum(map {sum(map {$_->[1]-$_->[0]} @{$ItvHref->{$_}})} keys %$ItvHref);
}

#--- set bed's notes ---
sub set_note{
    my $bed = shift;
    my %parm = @_;
    $bed->{note} = $parm{note};
}

#--- return bed's notes ---
sub note{
    my $bed = shift;
    return $bed->{note} || 'N/A';
}

1; ## tell the perl script the successful access of this module.
