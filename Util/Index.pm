package BioFuse::Util::Index;

use strict;
use warnings;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              Pos2Idx
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Util::Index';
#----- version --------
$VERSION = "0.03";
$DATE = '2018-10-31';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        Pos2Idx
                     /;

#--- get posWinIdx for given pos ---
sub Pos2Idx{
    # options
    shift if ($_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $winSize = $parm{winSize} || 1000;
    return int( $parm{pos} / $winSize );
}

#--- 
1; ## tell the perl script the successful access of this module.
