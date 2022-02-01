package BioFuse::Visual::BioInfo::Depth;

use strict;
use warnings;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              draw_depth_spectrum
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Visual::BioInfo::Depth';
#----- version --------
$VERSION = "0.05";
$DATE = '2018-12-03';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @function_list = qw/
                        draw_depth_spectrum
                     /;

# PLS also check BioFuse::BioInfo::Depth
# assume winDepth_Href has such structure:
# winDepth_Href = {}
# winDepth_Href -> $tissue -> $window_NO -> 'pos2depth' -> pos = orig_depth # soon delete
# winDepth_Href -> $tissue -> $window_NO -> 'pos2adjdepth' -> pos = adjust_depth
# winDepth_Href -> $tissue -> $window_NO -> 'mean_adjdepth' = adjust_depth_mean
# winDepth_Href -> $tissue -> $window_NO -> 'showdepth' = smooth_depth_to_show

#--- draw depth spectrum of each tissue ---
sub draw_depth_spectrum{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $svg_obj = $parm{svg_obj};
    my $winDepth_Href = $parm{winDepth_Href}; # next key is tissue
    my $smoDepthKey = $parm{smoDepthKey} || 'showdepth';
    my $maxWinNOtoDraw = $parm{maxWinNOtoDraw};
    my $axisZX = $parm{axisZX};
    my $axisZY = $parm{axisZY};
    my $yResol = $parm{yResol};
    my $SpmCol_Href = $parm{SpmCol_Href}; # next key is tissue
    my $maxDraw1st = $parm{maxDraw1st} || 0;

    # draw spectrum with large one at the backend
    for my $window_NO ( 1 .. $maxWinNOtoDraw ){
        # first to draw larger one
        my @sort_tis =  $maxDraw1st
                       # max comes first
                       ? ( sort { $winDepth_Href->{$b}->{$window_NO}->{$smoDepthKey}
                                  <=>
                                  $winDepth_Href->{$a}->{$window_NO}->{$smoDepthKey}
                                } keys %$winDepth_Href )
                       # case -> ctrl
                       : ( sort keys %$winDepth_Href );
        # draw spectrum
        for my $tis_idx ( 0 .. $#sort_tis ){
            my $tissue = $sort_tis[$tis_idx];
            my $opacity = ($tis_idx == 0) ? 1 : 0.8;
            my $smoDepth = $winDepth_Href->{$tissue}->{$window_NO}->{$smoDepthKey};
            next if( !defined $smoDepth || $smoDepth == 0 );
            $svg_obj->line( x1 => $axisZX + $window_NO,
                            y1 => $axisZY,
                            x2 => $axisZX + $window_NO,
                            y2 => $axisZY - $smoDepth / $yResol,
                            stroke => $SpmCol_Href->{$tissue},
                            'stroke-width' => 1,
                            'stroke-linecap' => 'round',
                            opacity => $opacity
                          );
        }
    }
}

#---
1; ## tell the perl script the successful access of this module.
