package BioFuse::Util::Array;

use strict;
use warnings;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              binarySearch
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Util::Array';
#----- version --------
$VERSION = "0.33";
$DATE = '2018-12-03';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @functoion_list = qw/
                        binarySearch
                     /;

#--- binary search in a given numeric ascending array ---
# it always assigns the closest lower number
sub binarySearch{
    shift if(@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $query = $parm{query};
    my $array = $parm{array}; # ref of array

    my ($L_i, $U_i) = (0, @$array - 1); # lower, upper boundary index of searching interval
    my $i;
    while($L_i <= $U_i){
        $i = int(($L_i + $U_i) / 2);
        if($array->[$i] < $query){
            $L_i = $i + 1;
        }
        elsif($array->[$i] > $query){
            $U_i = $i - 1;
        }
        else{
            return $i + 1;
        }
    }
    return $L_i;
}

1; ## tell the perl script the successful access of this module.
