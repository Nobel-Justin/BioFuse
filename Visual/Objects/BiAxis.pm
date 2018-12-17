package BioFuse::Visual::Objects::BiAxis;

use strict;
use warnings;
use Data::Dumper;
use List::Util qw/ min sum /;
use BioFuse::Util::Log qw/ warn_and_exit /;
use BioFuse::Visual::SVG_Util::RadSys qw/ $deg2rad get_coordinate_on_circle /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::Visual::Objects::BiAxis';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-12-16';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        calc_orig
                        get_origX
                        get_origY
                        get_axis
                        valueToSVGcoord
                        draw
                     /;

#--- structure of object
# bi_axis -> axis = [axis_1, axis_2]
# bi_axis -> orig = [orig_X, orig_Y];
# bi_axis -> drawn = 0/1

#--- construction of object ---
sub new{
    shift;
    my %parm = @_;

    my $bi_axis = {};
    # location
    $bi_axis->{axis} = [$parm{axis_1}, $parm{axis_2}];
    bless($bi_axis);
    $bi_axis->calc_orig;
    return $bi_axis;
}

#--- calculate the SVG X-Y coordinate of orig point of bi_axis ---
sub calc_orig{
    my $bi_axis = shift;

    my $rad1 = $bi_axis->{axis}->[0]->get_headRad;
    my $rad2 = $bi_axis->{axis}->[1]->get_headRad;
    # check
    if(abs(cos($rad1-$rad2)) > 0.8){
        warn_and_exit "<ERROR>\tAxises are close to be parallel in BiAxis object.\n".Dumper($bi_axis);
    }
    my $sin1 = sin $rad1;
    my $cos1 = cos $rad1;
    my $sin2 = sin $rad2;
    my $cos2 = cos $rad2;
    # axis locations
    my $x1 = $bi_axis->{axis}->[0]->get_origX;
    my $y1 = $bi_axis->{axis}->[0]->get_origY;
    my $x2 = $bi_axis->{axis}->[1]->get_origX;
    my $y2 = $bi_axis->{axis}->[1]->get_origY;

    my $x = sprintf "%.3f",
            (   $sin1 * $sin2 * ($y2 - $y1)
              + $sin1 * $cos2 * $x2
              - $cos1 * $sin2 * $x1
            ) /
            (   $sin1 * $cos2
              - $cos1 * $sin2
            );
    my $y = sprintf "%.3f",
            (   $cos1 * $cos2 * ($x1 - $x2)
              + $sin1 * $cos2 * $y1
              - $cos1 * $sin2 * $y2
            ) /
            (   $sin1 * $cos2
              - $cos1 * $sin2
            );
    $bi_axis->{orig} = [$x, $y];
}

#--- return SVG X coordinate of orig point ---
sub get_origX{
    my $bi_axis = shift;
    return $bi_axis->{orig}->[0];
}

#--- return SVG Y coordinate of orig point ---
sub get_origY{
    my $bi_axis = shift;
    return $bi_axis->{orig}->[1];
}

#--- return the axis object cooresponding to given NO. ---
sub get_axis{
    my $bi_axis = shift;
    my %parm = @_;
    my $i = (($parm{no} || 1) - 1) % 2;
    return $bi_axis->{axis}->[$i];
}

#--- return SVG X-Y coordinate of given value ---
sub valueToSVGcoord{
    my $bi_axis = shift;
    my %parm = @_;
    my $valueAf = $parm{valueAf};
    my $adjustOuterAf = $parm{adjustOuterAf} || [0,0];

    my ($x1,$y1) = $bi_axis->{axis}->[0]->valueToSVGcoord(value => $valueAf->[0], adjustOuter => $adjustOuterAf->[0]);
    my ($x2,$y2) = $bi_axis->{axis}->[1]->valueToSVGcoord(value => $valueAf->[1], adjustOuter => $adjustOuterAf->[1]);
    my $x = $x1 + $x2 - $bi_axis->{orig}->[0];
    my $y = $y1 + $y2 - $bi_axis->{orig}->[1];
    return ($x, $y);
}

#--- draw bi-axis ---
sub draw{
    my $bi_axis = shift;
    my %parm = @_;
    my $svg_obj = $parm{svg_obj};

    # drawn before?
    unless($bi_axis->{drawn}){
        $bi_axis->get_axis(no=>$_)->draw_quick(svg_obj => $svg_obj) for (1,2);
    }
    # mark
    $bi_axis->{drawn} = 1;
}

1; ## tell the perl script the successful access of this module.
