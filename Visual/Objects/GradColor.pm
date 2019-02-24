package BioFuse::Visual::Objects::GradColor;

use strict;
use warnings;
use Data::Dumper;
use List::Util qw/max min/;
use BioFuse::Util::Log qw/ warn_and_exit /;
use BioFuse::Visual::SVG_Util::RectSysEle qw/ draw_a_parallelogram /;
use BioFuse::Visual::Objects::Axis;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::Visual::Objects::GradColor';
#----- version --------
$VERSION = "0.01";
$DATE = '2019-02-24';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        set_axis
                        get_stRGB
                        get_edRGB
                        get_stValue
                        get_edValue
                        valueToRGB
                        draw_quick
                        draw_body
                     /;

#--- structure of object
# gc -> rgb = {st=>[r,g,b], ed=>[r,g,b], gap=>[]}
# gc -> range = {st=>$stValue, ed=>$edValue, method=>'linear'}
# gc -> anchor = {x=>$anchorX, y=>$anchorY, anchor=>'middle/start/end'}
# gc -> shape = {type=>'rect', height=>$height, width=>$width}
# gc -> drawn = {body=>0/1, stub=>0/1, label=>0/1}
# gc -> axis = $gc_axis

#--- construction of object ---
sub new{
    shift;
    my %parm = @_;

    my $gc = {};
    # color
    $gc->{rgb}->{st} = [split /,/, $parm{stRGB} || '255,255,255'];
    $gc->{rgb}->{ed} = [split /,/, $parm{edRGB} || '255,0,0'];
    $gc->{rgb}->{gap}= [map {($gc->{rgb}->{ed}->[$_] - $gc->{rgb}->{st}->[$_])} (0..2)];
    # range
    $gc->{range}->{st} = $parm{stValue};
    $gc->{range}->{ed} = $parm{edValue};
    # anchor
    $gc->{anchor}->{x} = $parm{anchorX} || 100;
    $gc->{anchor}->{y} = $parm{anchorY} || 100;
    $gc->{anchor}->{anchor} = $parm{anchor} || 'start';
    # shape
    $gc->{shape}->{type} = 'rect';
    $gc->{shape}->{height} = $parm{height} || 10;
    $gc->{shape}->{width} = $parm{width} || 40;
    # draw sign
    $gc->{drawn} = {body=>0, axis=>0};

    bless($gc);
    return $gc;
}

#--- set axis ---
sub set_axis{
    my $gc = shift;
    my %parm = @_;
    my $clockwise = $parm{clockwise} || 0;
    my $label = $parm{label} || '';
    my $labelFZ = $parm{labelFZ} || 10;
    my $stubText = $parm{stubText} || '';
    # axis
    my $ratio = $gc->{anchor}->{anchor} eq 'end'    ?   1 :
                $gc->{anchor}->{anchor} eq 'middle' ? 0.5 : 0;
    my $axis_X = $gc->{anchor}->{x} - $gc->{shape}->{width} * $ratio;
    my $axis_Y = $gc->{anchor}->{y} + $gc->{shape}->{height} * 0.75 * ($clockwise ? 1 : -1);
    $gc->{axis} = BioFuse::Visual::Objects::Axis->new(origP_X=>$axis_X, origP_Y=>$axis_Y, axisLen=>$gc->{shape}->{width}, headAng=>90);
    $gc->{axis}->set_tic(len=>3, width=>1, clockwise=>$clockwise);
    $gc->{axis}->add_resol(minValue=>$gc->{range}->{st}, maxValue=>$gc->{range}->{ed});
    $gc->{axis}->set_stub(stubStep=>$gc->{range}->{ed} - $gc->{range}->{st}, stubText=>$stubText);
    $gc->{axis}->set_label(text=>$label, fontsize=>$labelFZ);
}

#--- return stRGB ---
sub get_stRGB{
    my $gc = shift;
    return $gc->{rgb}->{st};
}

#--- return edRGB ---
sub get_edRGB{
    my $gc = shift;
    return $gc->{rgb}->{ed};
}

#--- return stValue ---
sub get_stValue{
    my $gc = shift;
    return $gc->{range}->{st};
}

#--- return edValue ---
sub get_edValue{
    my $gc = shift;
    return $gc->{range}->{ed};
}

#--- return rgb of given value ---
sub valueToRGB{
    my $gc = shift;
    my %parm = @_;
    my $value = $parm{value};
    # check
    if(    !$parm{allowOutOfRange}
        && !($value >= $gc->{range}->{st} && $value <= $gc->{range}->{ed})
    ){
        warn_and_exit "<ERROR>\tvalue $value out of range of GradColor object.\n"
                            ."\t".Dumper($gc);
    }
    # calculate
    my $ratio = ($value - $gc->{range}->{st}) / ($gc->{range}->{ed} - $gc->{range}->{st});
    return [ map {min(max(int($gc->{rgb}->{st}->[$_] + $ratio * $gc->{rgb}->{gap}->[$_]),0),255)} (0..2) ];
}

#--- quick draw the gradColor ---
## embeded display of different parts
sub draw_quick{
    my $gc = shift;
    $gc->draw_body(@_);
    $gc->{axis}->draw_quick(@_, stubSciNum=>3);
}

#--- draw the body of gradColor ---
sub draw_body{
    my $gc = shift;
    my %parm = @_;
    my $svg_obj = $parm{svg_obj};

    # drawn before?
    return if $gc->{drawn}->{body};

    draw_a_parallelogram( svg_obj => $svg_obj,
                          x => $gc->{axis}->get_origX + $gc->{shape}->{width} / 2,
                          y => $gc->{anchor}->{y},
                          colorGradOrit => 'h',
                          colorGradStColor => 'rgb('.join(',',@{$gc->{rgb}->{st}}).')',
                          colorGradStOpacity => 1,
                          colorGradEdColor => 'rgb('.join(',',@{$gc->{rgb}->{ed}}).')',
                          colorGradEdOpacity => 1,
                          head_bottom_side_len => $gc->{shape}->{width},
                          left_right_side_len => $gc->{shape}->{height},
                          boundary_color => 'none'
                        );
    # mark
    $gc->{drawn}->{body} = 1;
}

1; ## tell the perl script the successful access of this module.
