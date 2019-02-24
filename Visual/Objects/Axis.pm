package BioFuse::Visual::Objects::Axis;

use strict;
use warnings;
use List::Util qw/ min sum first any /;
use POSIX qw/ ceil /;
use Data::Dumper;
use BioFuse::Util::Log qw/ warn_and_exit /;
use BioFuse::Util::Interval qw/ Get_Two_Seg_Olen /;
use BioFuse::Visual::SVG_Util::RadSys qw/ $PI $deg2rad get_coordinate_on_circle /;
use BioFuse::Visual::SVG_Util::Font qw/ show_text_in_line /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::Visual::Objects::Axis';
#----- version --------
$VERSION = "0.02";
$DATE = '2019-02-24';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        add_resol
                        set_label
                        set_stub
                        add_stub
                        set_tic
                        set_resolItvGap
                        validate_resol
                        update_axisLen
                        get_origX
                        get_origY
                        get_axisLen
                        get_resol
                        get_stub
                        get_headAng
                        get_headRad
                        get_origValue
                        valueToAxisDist
                        valueToSVGcoord
                        extendCoord
                        draw_quick
                        draw_axisBodyLine
                        draw_stub
                        draw_label
                     /;

#--- structure of object
# definition of attributes: http://ploticus.sourceforge.net/doc/axis.html
# axis -> axisLen = $axisLen
# axis -> origPos = {x=>$origP_X, y=>$origP_Y, value=>$origP_value, dist_offset=>$dist_offset}
# axis -> headAng = $headAng ## it's clock-wise angle (0-360) from zero o'clock
# axis -> headRad = $headRad ## radian of $headAng
# axis -> resol = [{min=>$min, max=>$max, method=>'linear', resol=>}, {}, ..]
# axis -> label = {text=>$text, gapToBodyLine=>$gapToBodyLine}
# axis -> stub = [$stub1, $stub2, ..]
# axis -> stubText = $stubText
# axis -> stubAttr = {ticGap=>$stubTicGap, vtToBL=>$verticalToBL, plToSX=>$parallelToSX, sciNum=>$stubSciNum}
# axis -> tic = {len=>$len, width=>$width, clockwise=>0/1}
# axis -> resolItvGap = $resolItvGap
# axis -> drawn = {bodyline=>0/1, stub=>0/1, legend=>0/1, label=>0/1}
# axis -> font = {family=>$family, size=>{label=>$xx,stub=>$xx}}
# axis -> stroke = {color=>$strokeCol, width=>$strokeWid}

#--- construction of object ---
sub new{
    shift;
    my %parm = @_;

    my $axis = {};
    # location
    $axis->{origPos} = {x=>$parm{origP_X}||100, y=>$parm{origP_Y}||100, value=>$parm{origP_value}||0, dist_offset=>undef};
    $axis->{headAng} = $parm{headAng} || 0;
    $axis->{headRad} = $parm{headAng} * $deg2rad;
    # basic
    $axis->{axisLen} = $parm{axisLen} || 100;
    $axis->{resol} = $parm{resol} || [];
    $axis->{resolItvGap} = $parm{resolItvGap} || 0;
    $axis->{stub} = $parm{stub} || [];
    $axis->{drawn} = {bodyline=>0, stub=>0, legend=>0, label=>0};
    # global appearance
    $axis->{font} = {family=>$parm{fontfam}||'Arial', size=>{}};
    $axis->{stroke} = {color=>$parm{strokeCol}||'black', width=>$parm{strokeWid}||1};

    # check
    unless($axis->{headAng} >= 0 && $axis->{headAng} < 360){
        warn_and_exit "<ERROR>\tThe head angle must be [0,360).\n";
    }

    bless($axis);
    return $axis;
}

#--- set resolution ---
sub add_resol{
    my $axis = shift;
    my %parm = @_;
    my $minValue = $parm{minValue} || 0;
    my $maxValue = $parm{maxValue} || 10;
    my $pixelLen = $parm{pixelLen} || $axis->{axisLen};
    my $resol = $parm{resol} || (($maxValue - $minValue) / $pixelLen);
    # record
    push @{$axis->{resol}}, {min => $minValue, max => $maxValue, method => 'linear', resol => $resol};
    # validate resol
    $axis->validate_resol;
    # update axis length
    $axis->update_axisLen;
}

#--- set label of axis ---
sub set_label{
    my $axis = shift;
    my %parm = @_;
    my $fontsize = $parm{fontsize} || 12;
    $axis->{label}->{text} = $parm{text};
    $axis->{font}->{size}->{label} = $parm{fontsize} || 14;
    $axis->{label}->{gapToBodyLine} = $parm{gapToBodyLine} || ($axis->{font}->{size}->{label} * 1.5);
}

#--- set stubs on axis ---
sub set_stub{
    my $axis = shift;
    my %parm = @_;
    my $stubStep = $parm{stubStep} || 0;
    my $stubPixSpan = $parm{stubPixSpan} || 20;
    my $updateMax = $parm{updateMax} || 0;
    my $skipItvBoundary = $parm{skipItvBoundary} || 0;
    my $skipEndBoundary = $parm{skipEndBoundary} || 0;
    my $stubText = $parm{stubText} || undef; # e.g., depth unit 'X'
    my $fontsize = $parm{fontsize} || 8;
    my $verticalToBL = $parm{verticalToBL} || 0; # vertical to axis body line
    my $parallelToSX = $parm{parallelToSX} || 0; # parallel to SVG X-axis, prior to verticalToBL
    my $stubTicGap = $parm{stubTicGap} || 2;

    # candidate stub step
    my @candStubStep;
    for my $r (0..9) { push @candStubStep, $_ * (10 ** $r) for (1..6,8) }
    # set stubs of each resol-interval
    my $resolItvCount = scalar @{$axis->{resol}};
    for my $i (0 .. $resolItvCount-1){
        my $hf = $axis->{resol}->[$i];
        # allocate proper stub step
        unless($stubStep){
            $stubStep = $hf->{resol} * $stubPixSpan;
            my $idx = first {$candStubStep[$_] <= $stubStep && $candStubStep[$_+1] >= $stubStep} (0 .. $#candStubStep-1);
            $stubStep = defined $idx ? $candStubStep[$idx] : ( $stubStep > 10 ? $candStubStep[-1] : $candStubStep[0] );
        }
        # update max value if set
        $hf->{max} = ceil($hf->{max} / $stubStep) * $stubStep if $updateMax;
        # record stub
        for(my $stub = $hf->{min}; $stub <= $hf->{max}; $stub += $stubStep){
            $axis->add_stub(stub => $stub) unless $skipItvBoundary && ($stub == $hf->{min} || $stub == $hf->{max});
        }
    }
    # add skipEndBoundary
    unless($skipEndBoundary){
        $axis->add_stub(stub => $axis->{resol}->[0]->{min});
        $axis->add_stub(stub => $axis->{resol}->[-1]->{max});
    }
    # validate resol
    $axis->validate_resol;
    # update axis length
    $axis->update_axisLen;
    # record text shown at the last stub
    $axis->{stubText} = $stubText;
    $axis->{font}->{size}->{stub} = $fontsize;
    # stub attributes
    $verticalToBL ||= $parallelToSX;
    $axis->{stubAttr} = {ticGap => $stubTicGap, vtToBL => $verticalToBL, plToSX => $parallelToSX};
}

#--- manually add stub ---
sub add_stub{
    my $axis = shift;
    my %parm = @_;
    my $stub = $parm{stub};
    push $axis->{stub}, $stub unless any {$_==$stub} @{$axis->{stub}};
    @{$axis->{stub}} = sort {$a<=>$b} @{$axis->{stub}};
}

#--- set tics on axis ---
sub set_tic{
    my $axis = shift;
    my %parm = @_;
    $axis->{tic} = {len => $parm{len}||3, width => $parm{width}||1, clockwise => $parm{clockwise}||0};
}

#--- set gap between resol intervals ---
sub set_resolItvGap{
    my $axis = shift;
    my %parm = @_;
    $axis->{resolItvGap} = $parm{resolItvGap} if defined $parm{resolItvGap};
}

#--- check the resol-interval has no overlap with each other ---
sub validate_resol{
    my $axis = shift;
    my %parm = @_;
    my $offset = $parm{offset} || 11E-5;
    # sort
    @{$axis->{resol}} = sort {$a->{min}<=>$b->{min}} @{$axis->{resol}};
    # pairwise comparison
    my $resolItvCount = scalar @{$axis->{resol}};
    for my $i (0 .. $resolItvCount-1){
        my $i_hf = $axis->{resol}->[$i];
        for my $j ($i+1 .. $resolItvCount-1){
            my $j_hf = $axis->{resol}->[$j];
            warn_and_exit "<ERROR>\tfind overlaps between resol interval.\n".Dumper($i_hf).Dumper($j_hf)
              if Get_Two_Seg_Olen($i_hf->{min}, $i_hf->{max}, $j_hf->{min}, $j_hf->{max}, $offset);
        }
    }
}

#--- update axis length according to resol interval(s) ---
sub update_axisLen{
    my $axis = shift;
    $axis->{axisLen} = sum( map{ (($_->{max} - $_->{min}) / $_->{resol} + $axis->{resolItvGap}) } @{$axis->{resol}} );
    $axis->{axisLen} -= $axis->{resolItvGap};
    $axis->{axisLen} = ceil $axis->{axisLen};
}

#--- return SVG X coordinate of orig point ---
sub get_origX{
    my $axis = shift;
    return $axis->{origPos}->{x};
}

#--- return SVG Y coordinate of orig point ---
sub get_origY{
    my $axis = shift;
    return $axis->{origPos}->{y};
}

#--- return axis length ---
sub get_axisLen{
    my $axis = shift;
    return $axis->{axisLen};
}

#--- return resolution ---
sub get_resol{
    my $axis = shift;
    return $axis->{resol};
}

#--- return stub ---
sub get_stub{
    my $axis = shift;
    return $axis->{stub};
}

#--- return axis head angle ---
sub get_headAng{
    my $axis = shift;
    return $axis->{headAng};
}

#--- return axis head radian ---
sub get_headRad{
    my $axis = shift;
    return $axis->{headRad};
}

#--- return value cooresponding to axis orig point ---
sub get_origValue{
    my $axis = shift;
    return $axis->{origPos}->{value};
}

#--- get given value's shift-distance along the axis ---
sub valueToAxisDist{
    my $axis = shift;
    my %parm = @_;
    my $value = $parm{value};
    my $adjustOuter = $parm{adjustOuter} || 0;
    my $absolutDist = $parm{absolutDist} || 0;

    if(    !$absolutDist
        && !defined $axis->{origPos}->{dist_offset}
    ){
        $axis->{origPos}->{dist_offset} = $axis->valueToAxisDist(value => $axis->{origPos}->{value}, absolutDist => 1);
    }

    # calculate value's shift-distance along the axis
    my $dist = $axis->{resolItvGap} * -1;
    my $findItv = 0;
    my $resolItvCount = scalar @{$axis->{resol}};
    my $i;
    for ($i = 0; $i < $resolItvCount; $i++){
        my $hf = $axis->{resol}->[$i];
        last if $value < $hf->{min};
        $findItv = 1 if $value <= $hf->{max};
        $dist += $axis->{resolItvGap};
        $dist += (min($hf->{max}, $value) - $hf->{min}) / $hf->{resol};
        last if $findItv;
    }
    # outer-value or gapped
    unless($findItv){
        if(    $adjustOuter
            && ($i == 0 || $i == $resolItvCount)
        ){
            $dist = $dist > 0 ? $axis->{axisLen} : 1;
        }
        else{
            warn_and_exit "<ERROR>\tvalue ($value) is out of resol-interval\n".Dumper($axis->{resol});
        }
    }
    # rative distance
    return $absolutDist ? $dist : $dist - $axis->{origPos}->{dist_offset};
}

#--- get given value's SVG X-Y coordinates along the axis ---
sub valueToSVGcoord{
    my $axis = shift;
    return get_coordinate_on_circle(cx => $axis->{origPos}->{x}, cy => $axis->{origPos}->{y}, rad => $axis->{headRad}, radius => $axis->valueToAxisDist(@_));
}

#--- return SVG X-Y coordinates ---
## extending given length to given angle from given position
## the angle is clock-wise angle (0-360) from zero o'clock perpendicular to this axis
sub extendCoord{
    my $axis = shift;
    my %parm = @_;
    my $coordAf = $parm{coordAf};
    my $relateAng = $parm{relateAng} || 0;
    my $distance = $parm{distance} || 3;

    my $radian = $relateAng * $deg2rad;
    my $radCal = $radian + $axis->{headRad} - $PI / 2;
    my $x = $coordAf->[0] + $distance * sin $radCal;
    my $y = $coordAf->[1] - $distance * cos $radCal;
    return ($x, $y);
}

#--- quick draw the axis ---
## embeded display of different parts
sub draw_quick{
    my $axis = shift;
    $axis->draw_axisBodyLine(@_);
    $axis->draw_stub(@_);
    $axis->draw_label(@_);
}

#--- draw the body line of axis ---
sub draw_axisBodyLine{
    my $axis = shift;
    my %parm = @_;
    my $svg_obj = $parm{svg_obj};

    # drawn before?
    return if $axis->{drawn}->{bodyline};

    # if set, update stroke attributes
    $axis->{stroke}->{color} = $parm{strokeCol} if defined $parm{strokeCol};
    $axis->{stroke}->{width} = $parm{strokeWid} if defined $parm{strokeWid};

    # draw
    my @strokeParm = (stroke => $axis->{stroke}->{color}, 'stroke-width' => $axis->{stroke}->{width}, 'stroke-linecap' => 'round');
    my $last_edXY = undef;
    for my $hf (@{$axis->{resol}}){
        my ($x_st, $y_st) = $axis->valueToSVGcoord(value => $hf->{min});
        my ($x_ed, $y_ed) = $axis->valueToSVGcoord(value => $hf->{max});
        # line
        $svg_obj->line(x1 => $x_st, y1 => $y_st, x2 => $x_ed, y2 => $y_ed, @strokeParm);
        # gap circle sign
        if(    $axis->{resolItvGap} > $axis->{stroke}->{width} * 6
            && defined $last_edXY
        ){
            my $time = 3;
            my $xGapU = ($x_st - $last_edXY->[0]) / $time;
            my $yGapU = ($y_st - $last_edXY->[1]) / $time;
            for my $t (1 .. $time-1){
                my $cx = $last_edXY->[0] + $xGapU * $t;
                my $cy = $last_edXY->[1] + $yGapU * $t;
                $svg_obj->circle(cx => $cx, cy => $cy, r => $axis->{stroke}->{width}, fill => $axis->{stroke}->{color});
            }
        }
        # update
        $last_edXY = [$x_ed, $y_ed];
    }
    # mark
    $axis->{drawn}->{bodyline} = 1;
}

#--- draw stub along the axis ---
sub draw_stub{
    my $axis = shift;
    my %parm = @_;
    my $svg_obj = $parm{svg_obj};

    # drawn before?
    return if $axis->{drawn}->{stub};

    # if set, update
    $axis->{stubAttr}->{vtToBL} = $parm{verticalToBL} if defined $parm{verticalToBL};
    $axis->{stubAttr}->{plToSX} = $parm{parallelToSX} if defined $parm{parallelToSX};
    $axis->{stubAttr}->{ticGap} = $parm{stubTicGap} if defined $parm{stubTicGap};
    $axis->{stubAttr}->{sciNum} = $parm{stubSciNum} if defined $parm{stubSciNum};

    # draw
    my $verticalToBL = $axis->{stubAttr}->{vtToBL};
    my $parallelToSX = $axis->{stubAttr}->{plToSX};
    # stub SVG-text-anchor-point Y coordinate to adjust
    my $height_adjust;
    if($axis->{headAng} <= 180){
        $height_adjust = $verticalToBL ? 1 : $axis->{tic}->{clockwise} ? 2 : 0;
    }
    else{
        $height_adjust = $verticalToBL ? 1 : $axis->{tic}->{clockwise} ? 0 : 2;
    }
    # stub SVG-text to rotate
    my $rotate_degree = $parallelToSX ? 0 : ($axis->{headAng} - 90 * ($verticalToBL ? 0 : 1));
    while($rotate_degree >  180){ $rotate_degree -= 360}
    while($rotate_degree < -180){ $rotate_degree += 360}
    while($rotate_degree >   90){ $rotate_degree -= 180}
    while($rotate_degree <  -90){ $rotate_degree += 180}
    # stub-tic angle related to body line
    my $relateAng = $axis->{tic}->{clockwise} ? 180 : 0;
       $relateAng += 180 * (($axis->{headAng}>270 || abs($axis->{headAng}-135)<=45) ? 1 : 0) - $axis->{headAng} if $parallelToSX;
    # stub SVG-text anchor
    my $text_anchor;
    if($axis->{headAng} <= 90 || $axis->{headAng} >= 270){
        $text_anchor = $verticalToBL ? ($axis->{tic}->{clockwise} ? 'start' : 'end') : 'middle';
    }
    else{
        $text_anchor = $verticalToBL ? ($axis->{tic}->{clockwise} ? 'end' : 'start') : 'middle';
    }
    # draw each stub
    my $stubCount = @{$axis->{stub}};
    for my $i (0 .. $stubCount-1){
        my $stub = $axis->{stub}->[$i];
        my ($stub_x, $stub_y) = $axis->valueToSVGcoord(value => $stub);
        # tic
        my ($ticE_x, $ticE_y) = $axis->extendCoord(coordAf => [$stub_x, $stub_y], relateAng => $relateAng, distance => $axis->{tic}->{len});
        $svg_obj->line( x1 => $stub_x, y1 => $stub_y, x2 => $ticE_x, y2 => $ticE_y,
                        stroke => $axis->{stroke}->{color}, 'stroke-width' => $axis->{tic}->{width}, 'stroke-linecap' => 'round'
                      );
        # stub
        ($stub = sprintf ("%.$axis->{stubAttr}->{sciNum}e",$stub)) =~ s/\.?0*e\+0?/E/ if $axis->{stubAttr}->{sciNum} && $stub > 10**$axis->{stubAttr}->{sciNum};
        my ($text_x, $text_y) = $axis->extendCoord(coordAf => [$stub_x, $stub_y], relateAng => $relateAng, distance => $axis->{tic}->{len}+$axis->{stubAttr}->{ticGap});
        $stub = $text_anchor eq 'end' ? "$axis->{stubText} $stub" : "$stub $axis->{stubText}" if $i == $stubCount-1 && defined $axis->{stubText};
        show_text_in_line( svg_obj => $svg_obj, text_x => $text_x, text_y => $text_y, text => $stub, text_anchor => $text_anchor,
                           font_family => $axis->{font}->{family}, font_size => $axis->{font}->{size}->{stub},
                           height_adjust => $height_adjust, rotate_degree => $rotate_degree
                         );
    }
    # mark
    $axis->{drawn}->{stub} = 1;
}

#--- draw label along the axis ---
sub draw_label{
    my $axis = shift;
    my %parm = @_;
    my $svg_obj = $parm{svg_obj};

    # drawn before?
    return if $axis->{drawn}->{label} || !$axis->{label}->{text};

    # draw
    my $verticalToBL = 0;
    # label SVG-text-anchor-point Y coordinate to adjust
    my $height_adjust;
    if($axis->{headAng} <= 180){
        $height_adjust = $verticalToBL ? 1 : $axis->{tic}->{clockwise} ? 2 : 0;
    }
    else{
        $height_adjust = $verticalToBL ? 1 : $axis->{tic}->{clockwise} ? 0 : 2;
    }
    # label SVG-text to rotate
    my $rotate_degree = $axis->{headAng} - 90 * ($verticalToBL ? 0 : 1);
    while($rotate_degree >  180){ $rotate_degree -= 360}
    while($rotate_degree < -180){ $rotate_degree += 360}
    while($rotate_degree >   90){ $rotate_degree -= 180}
    while($rotate_degree <  -90){ $rotate_degree += 180}
    # stub-tic angle related to body line
    my $relateAng = $axis->{tic}->{clockwise} ? 180 : 0;
    # label SVG-text anchor
    my $text_anchor;
    if($axis->{headAng} <= 90 || $axis->{headAng} >= 270){
        $text_anchor = $verticalToBL ? ($axis->{tic}->{clockwise} ? 'start' : 'end') : 'middle';
    }
    else{
        $text_anchor = $verticalToBL ? ($axis->{tic}->{clockwise} ? 'end' : 'start') : 'middle';
    }

    # draw label
    my $stubCount = @{$axis->{stub}};
    my $label_stub =   $stubCount % 2
                     ? $axis->{stub}->[int($stubCount/2)]
                     : ($axis->{stub}->[$stubCount/2] + $axis->{stub}->[$stubCount/2-1]) / 2;
    my ($label_stub_x, $label_stub_y) = $axis->valueToSVGcoord(value => $label_stub);
    my ($label_text_x, $label_text_y) = $axis->extendCoord(coordAf => [$label_stub_x, $label_stub_y], relateAng => $relateAng, distance => $axis->{label}->{gapToBodyLine});
    show_text_in_line( svg_obj => $svg_obj, text_x => $label_text_x, text_y => $label_text_y, text => $axis->{label}->{text}, text_anchor => $text_anchor,
                       font_family => $axis->{font}->{family}, font_size => $axis->{font}->{size}->{label},
                       height_adjust => $height_adjust, rotate_degree => $rotate_degree
                     );
    # mark
    $axis->{drawn}->{label} = 1;
}

1; ## tell the perl script the successful access of this module.
