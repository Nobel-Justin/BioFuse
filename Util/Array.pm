package BioFuse::Util::Array;

use strict;
use warnings;
use List::Util qw/ min max /;
use BioFuse::Util::Log qw/ warn_and_exit /;
use Data::Dumper;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              mergeOverlap
              Get_Two_Seg_Olen
              binarySearch
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Util::Array';
#----- version --------
$VERSION = "0.32";
$DATE = '2018-11-15';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @functoion_list = qw/
                        mergeOverlap
                        Get_Two_Seg_Olen
                        binarySearch
                     /;

#--- Get the overlapped length of two intervals ---
# Note: only works for intervals with integer boundaries
# if want to use it for demical boundaries, inter the fifth option
sub Get_Two_Seg_Olen{
    shift if(@_ && $_[0] =~ /$MODULE_NAME/);
    my ($s1,$e1,$s2,$e2,$offset) = @_;
    $offset = 1 unless defined($offset);
    my $overlap_len = min($e1,$e2) - max($s1,$s2) + $offset;
    return ($overlap_len<0) ? 0 : $overlap_len;
}

#--- binary search in a given numeric ascending array ---
# it always assigns the closest lower number
sub binarySearch {
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

#--- merge overlapped region ---
# default max_dist is 0, which means will not merge adjacent region, such as [10,19],[20,30]
# when set mergeAdjacent, max_dist is set as at least one, then adjacent regions will be merged
sub mergeOverlap{
    shift if(@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $regionAref = $parm{regionAref}; # [ [st1,ed1],[st2,ed2],[st3,ed3],... ]
    my $max_dist = $parm{max_dist} || 0;
    my $mergeAdjacent = $parm{mergeAdjacent} || 0;

    # set max_dist as at least one
    $max_dist = max(1, $max_dist) if $mergeAdjacent;
    # sort by end-pos
    @$regionAref = sort {$a->[1] <=> $b->[1]} @$regionAref;
    # merge overlap
    for (my $i = scalar(@$regionAref)-1; $i > 0; $i--){
        # check credible interval
        if($regionAref->[$i]->[0] > $regionAref->[$i]->[1]){
            warn_and_exit "<ERROR>\tindex $i has reversed interval: [@{$regionAref->[$i]}]\n";
        }
        # check overlap
        if( $regionAref->[$i]->[0] <= $regionAref->[$i-1]->[1] + $max_dist ){
            my $new_st = min($regionAref->[$i]->[0],$regionAref->[$i-1]->[0]);
            splice(@$regionAref, $i-1, 2, [$new_st, $regionAref->[$i]->[1]]);
        }
    }
    # sort by st-pos
    @$regionAref = sort {$a->[0] <=> $b->[0]} @$regionAref;
}

1; ## tell the perl script the successful access of this module.
