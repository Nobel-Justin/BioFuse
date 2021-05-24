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
$VERSION = "0.02";
$DATE = '2021-05-21';

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
                        set_seq
                        seq
                     /;

#--- structure of object
# refseg -> id = $id
# refseg -> length = $length
# refseg -> note = $note
# refseg -> seq = {orig=>orig_seq, h1=>h1_seq, .., hx=>hx_seq}

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

#--- seq sequence ---
sub set_seq{
    my $refseg = shift;
    my %parm = @_;
    unless(defined $parm{seqID} && defined $parm{seq}){
        cluck_and_exit "<ERROR>\trequires both seqID and seq.\n";
    }
    $refseg->{seq}->{$parm{seqID}} = $parm{seq};
}

#--- return sequence ---
sub seq{
    my $refseg = shift;
    my %parm = @_;
    unless(defined $parm{seqID} && exists $refseg->{seq}->{$parm{seqID}}){
        cluck_and_exit "<ERROR>\tcannot find $parm{seqID} seq.\n";
    }
    return $refseg->{seq}->{$parm{seqID}};
}

1; ## tell the perl script the successful access of this module.
