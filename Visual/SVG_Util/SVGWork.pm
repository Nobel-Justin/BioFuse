package BioFuse::Visual::SVG_Util::SVGWork;

use strict;
use warnings;
use SVG;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              initialize_SVG_obj
              output_SVG_file
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Visual::SVG_Util::SVGWork';
#----- version --------
$VERSION = "0.02";
$DATE = '2021-01-11';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @function_list = qw/
                        initialize_SVG_obj
                        output_SVG_file
                     /;

#--- initialize SVG object ---
sub initialize_SVG_obj{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $bg_width  = $parm{bg_width} || '800';
    my $bg_height = $parm{bg_height} || '600';
    my $bg_col = $parm{bg_col} || 'white';
    my $bg_stroke_col = $parm{bg_stroke_col} || 'none';
    my $AUTHOR = $parm{AUTHOR} || 'anonymous';
    my $EMAIL = $parm{EMAIL} || 'test@test.com';
    my $skip_bg = $parm{skip_bg} || 0;

    my $svg_obj = SVG->new( width=>$bg_width, height=>$bg_height, author=>$AUTHOR, 'author-mail'=>$EMAIL);
    $svg_obj->rect(x=>0, y=>0, width=>$bg_width, height=>$bg_height, fill=>$bg_col, stroke=>$bg_stroke_col) unless $skip_bg;
    return $svg_obj;
}

#--- output SVG file ---
sub output_SVG_file{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $svg_obj = $parm{svg_obj};
    my $svg_file = $parm{svg_file};

    open (SVGFILE,">$svg_file") || die"fail write SVG file: $!\n";
    print SVGFILE $svg_obj->xmlify;
    close SVGFILE;
}

#--- 
1; ## tell the perl script the successful access of this module.
