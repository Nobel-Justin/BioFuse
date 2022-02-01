package BioFuse::Visual::SVG_Util::Color;

use strict;
use warnings;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              %COLOR_DB
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Visual::SVG_Util::Color';
#----- version --------
$VERSION = "0.01";
$DATE = '2019-02-24';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

# url: http://www.december.com/html/spec/colorsvg.html
our %COLOR_DB=(
                1=>'blue',    2=>'lime',     3=>'gold',      4=>'tomato',     5=>'purple',
                6=>'brown',   7=>'green',    8=>'orange',    9=>'skyblue',   10=>'pink',
               11=>'orchid', 12=>'darkred', 13=>'seagreen', 14=>'chocolate', 15=>'tan',

               dark_col=>'blue|brown|green|purple|chocolate|darkred'
             );

#--------- functions in this pm --------#
my @function_list = qw//;


1; ## tell the perl script the successful access of this module.
