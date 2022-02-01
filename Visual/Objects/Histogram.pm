package BioFuse::Visual::Objects::Histogram;

use strict;
use warnings;
use Data::Dumper;
use BioFuse::Visual::SVG_Util::RectSysEle qw/ draw_a_parallelogram /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::Visual::Objects::Histogram';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-12-16';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @function_list = qw/
                        new
                        load_data
                        set_attr
                        draw
                     /;

#--- structure of object
# histogram -> bi_axis = $bi_axis
# histogram -> isStack = 0/1
# histogram -> pillarWid = $pillarWid
# histogram -> fill = [$color1, $color2, ..]
# histogram -> stroke = {color=>[$color1, $color2, ..], width=>[$width1, $width2, ..]}
# histogram -> data = [{vA=>$vA, vB=>[$vB1,$vB2,..], fill=>{i=>$color_i, ..}, strokeCol=>{i=>$color_i, ..}, strokeWid=>{i=>$width_i}}]

#--- construction of object ---
sub new{
    shift;
    my %parm = @_;

    my $histogram = {};
    $histogram->{bi_axis} = $parm{bi_axis};
    $histogram->{isStack} = $parm{isStack} || 0;
    $histogram->{pillarWid} = $parm{pillarWid} || 5;
    $histogram->{data} = [];
    $histogram->{fill} = ['none'];
    $histogram->{stroke} = {color => ['black'], width => [0.5]};
    bless($histogram);
    return $histogram;
}

#--- load data ---
sub load_data{
    my $histogram = shift;
    my %parm = @_;
    push @{$histogram->{data}}, { vA => $parm{vA},
                                  vB => $parm{vB} || [],
                                  fill => $parm{fill} || {},
                                  strokeCol => $parm{strokeCol} || {},
                                  strokeWid => $parm{strokeWid} || {}
                                };
}

#--- set attribute ---
sub set_attr{
    my $histogram = shift;
    my %parm = @_;
    $histogram->{fill} = $parm{fillCol} if defined $parm{fillCol};
    $histogram->{stroke}->{color} = $parm{strokeCol} if defined $parm{strokeCol};
    $histogram->{stroke}->{width} = $parm{strokeWid} if defined $parm{strokeWid};
}

#--- draw histogram ---
sub draw{
    my $histogram = shift;
    my %parm = @_;
    my $svg_obj = $parm{svg_obj};

    # draw bi-axis
    my $biAxis = $histogram->{bi_axis};
    $biAxis->draw(svg_obj => $svg_obj);

    # draw histogram
    my $vA_axis = $biAxis->get_axis(no=>1);
    my $vB_axis = $biAxis->get_axis(no=>2);
    my $vB_orig = $vB_axis->get_origValue;
    my $upleft_angle = 180 - ($vA_axis->get_headAng - $vB_axis->get_headAng);
    my $rotate_degree = $vB_axis->get_headAng + 90 - $upleft_angle;
    for my $Hf (sort {$a->{vA}<=>$b->{vA}} @{$histogram->{data}}){
        my $vA = $Hf->{vA};
        my $vB_count = scalar @{$Hf->{vB}};
        my @i = (0 .. $vB_count-1);
           @i = sort {$Hf->{vB}->[$b] <=> $Hf->{vB}->[$a]} @i unless $histogram->{isStack};
        my $stackSum = 0;
        for my $i (@i){
            my $vtop = $Hf->{vB}->[$i] + $stackSum;
            my $vbom = $vB_orig + $stackSum;
            $stackSum = $vtop if $histogram->{isStack};
            my ($xtop,$ytop) = $biAxis->valueToSVGcoord(valueAf=>[$vA,$vtop], adjustOuterAf=>[1,0]);
            my ($xbom,$ybom) = $biAxis->valueToSVGcoord(valueAf=>[$vA,$vbom], adjustOuterAf=>[1,0]);
            my $cx = ($xtop + $xbom) / 2;
            my $cy = ($ytop + $ybom) / 2;
            my $fill = $Hf->{fill}->{$i} || $histogram->{fill}->[$i] || $histogram->{fill}->[-1];
            my $strokeCol = $Hf->{strokeCol}->{$i} || $histogram->{stroke}->{color}->[$i] || $histogram->{stroke}->{color}->[-1];
            my $strokeWid = $Hf->{strokeWid}->{$i} || $histogram->{stroke}->{width}->[$i] || $histogram->{stroke}->{width}->[-1];
            my $height = $vB_axis->valueToAxisDist(value => $vtop) - $vB_axis->valueToAxisDist(value => $vbom);
               $height += 0.001 unless $height;
            draw_a_parallelogram( svg_obj => $svg_obj,
                                  x => $cx,
                                  y => $cy,
                                  fill_color => $fill,
                                  boundary_color => $strokeCol,
                                  boundary_width => $strokeWid,
                                  head_bottom_side_len => $histogram->{pillarWid},
                                  left_right_side_len => $height,
                                  upleft_angle => $upleft_angle,
                                  rotate_degree => $rotate_degree
                                );
        }
    }
}

1; ## tell the perl script the successful access of this module.
