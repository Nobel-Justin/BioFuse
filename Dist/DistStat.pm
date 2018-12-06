package BioFuse::Dist::DistStat;

use strict;
use warnings;
use List::Util qw/ min max sum /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              get_value_mean
              get_trimmed_mean
              engineer_Ntimes_SD_evaluation
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Dist::DistStat';
#----- version --------
$VERSION = "0.02";
$DATE = '2018-01-09';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        get_value_mean
                        get_trimmed_mean
                        engineer_Ntimes_SD_evaluation
                     /;

#--- return mean value ---
# default to use TMM (ratio=0), means to accept all records.
sub get_value_mean{

    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $value_Aref    = $parm{value_Aref};
    my $Eng_SD_Ntimes = $parm{Eng_SD_Ntimes} || 0;
    my $TrimRatio     = $parm{TrimRatio} || 0;
    my $ratio_digit   = $parm{ratio_digit} || 2;

    if( $Eng_SD_Ntimes ){
        my %Value2Count;
        $Value2Count{$_}++ for @$value_Aref;
        my $Dist_Href = &engineer_Ntimes_SD_evaluation(
                            Value2Count_Href => \%Value2Count,
                            SD_Ntimes => $Eng_SD_Ntimes,
                            ratio_digit => $ratio_digit
                        );
        return $Dist_Href->{value_mean};
    }
    else{
        return &get_trimmed_mean(
                    array_Aref => $value_Aref,
                    TrimRatio => $TrimRatio,
                    ratio_digit => $ratio_digit
                );
    }
}

#--- TMM mean method ---
sub get_trimmed_mean{

    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $array_Aref  = $parm{array_Aref};
    my $TrimRatio   = $parm{TrimRatio} || 0;
    my $revs_sort   = $parm{revs_sort} || 0;
    my $ratio_digit = $parm{ratio_digit} || 2;

    # check available
    if( abs($TrimRatio) >= 1 ){
        die "<ERROR>\tCannot use TrimRatio ($TrimRatio) for TMM operation.\n";
    }

    my $remain_count;
    my $remain_value_sum;

    if( $TrimRatio == 0 ){
        $remain_count = scalar(@$array_Aref);
        $remain_value_sum = sum(@$array_Aref);
    }
    else{
        my $count = scalar(@$array_Aref);
        my $thereshold = int(abs($TrimRatio) * $count);
        my @numbers = $revs_sort ? (sort {$b<=>$a} @$array_Aref) : (sort {$a<=>$b} @$array_Aref);
        # only left trimming
        if( $TrimRatio < 0 ){
            $remain_count = $count - $thereshold;
            $remain_value_sum = sum( @numbers[ $thereshold .. $#numbers ] );
        }
        # bilaterally trimming
        else{
            $remain_count = $count - 2 * $thereshold;
            $remain_value_sum = sum( @numbers[ $thereshold .. $#numbers-$thereshold ] );
        }
        # check available
        if( $remain_count <= 0 ){
            die "<ERROR>\tNo value left after TMM operation.\n";
        }
    }

    # return
    return sprintf "%.${ratio_digit}f", $remain_value_sum / $remain_count;
}

#--- engineering N-times-SD filtering ---
# it will modify source Value2Count Hash
sub engineer_Ntimes_SD_evaluation{

    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $Value2Count_Href = $parm{Value2Count_Href};
    my $SD_Ntimes   = exists $parm{SD_Ntimes} ? $parm{SD_Ntimes} : 3;
    my $ratio_digit = exists $parm{ratio_digit} ? $parm{ratio_digit} : 2;

    # engineering N-times-SD filtering
    ENG_FILTER_INFO:{
        ## initial
        my %value_peak = ( value => 0, count => 0 );
        my $value_mean = 0;
        my $value_sd = 0;
        my %value_median = ( value => 0, count => 0 );
        my $value_sum = 0;
        my $count_sum = 0;
        ## peak and sum of value
        for my $value (sort {$a<=>$b} keys %$Value2Count_Href){
            my $count = $Value2Count_Href->{$value};
            %value_peak = ( value => $value, count => $count ) if( $count > $value_peak{value} );
            $value_sum += $value * $count;
            $count_sum += $count;
        }
        ## value mean
        $value_mean = $value_sum / $count_sum;
        ## sd of value
        $value_sd += (($_ - $value_mean) ** 2) * $Value2Count_Href->{$_} for sort {$a<=>$b} keys %$Value2Count_Href;
        $value_sd = sqrt( $value_sd / $count_sum );
        ## N-times-SD filtering
        my $redo_bool = 0;
        for my $abnormal_value (grep abs( $_ - $value_mean ) > $SD_Ntimes * $value_sd , keys %$Value2Count_Href){
            # discard the abnormals
            delete $Value2Count_Href->{$abnormal_value};
            $redo_bool = 1;
        }
        # redo
        redo ENG_FILTER_INFO if($redo_bool);
        # median value
        my $count_accum = 0;
        for my $value (sort {$a<=>$b} keys %$Value2Count_Href){
            $count_accum += $Value2Count_Href->{$value};
            if( $count_accum >= $count_sum / 2 ){
                %value_median = ( value => $value, count => $Value2Count_Href->{$value} );
                last;
            }
        }
        # refine value_mean
        $value_mean = sprintf "%.${ratio_digit}f", $value_mean;
        $value_sd   = sprintf "%.${ratio_digit}f", $value_sd;
        # return
        return  {
                    value_peak_Href => \%value_peak,
                    value_median_Href => \%value_median,
                    value_mean => $value_mean,
                    value_sd => $value_sd,
                    count_sum => $count_sum
                };
    }
}

#--- 
1; ## tell the perl script the successful access of this module.
