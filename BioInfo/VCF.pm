package BioFuse::BioInfo::VCF;

use strict;
use warnings;
use List::Util qw/ sum /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              check_I16_alt
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::VCF';
#----- version --------
$VERSION = "0.01";
$DATE = '2021-12-16';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @function_list = qw/
                        check_I16_alt
                     /;

#--- determine whether alt is ok based on I16 information ---
#--- this provides merged situation of the non-ref Q13 bases ---
sub check_I16_alt{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $I16Af = $parm{I16Af};
    my $min_mTD = $parm{min_mTD}; # min mean Tail Dist
    my $min_mBQ = $parm{min_mBQ}; # min mean BaseQ
    my $min_mMQ = $parm{min_mMQ}; # min mean MapQ
    my $only_TD = $parm{only_TD}; # only judge read Tail Dist

    # these are all based on Q13
    my $alt_readC = sum(@$I16Af[2,3]); # alt-supported reads count
    my $alt_mBQ = $I16Af->[6]  / $alt_readC; # non-ref mean baseQ
    my $alt_mMQ = $I16Af->[10] / $alt_readC; # non-ref mean MapQ
    my $alt_mTD = $I16Af->[14] / $alt_readC; # non-ref mean reads tail Dist

    # only judge the read tail distance
    if($only_TD){
        if($alt_mTD < $min_mTD){
            return 0;
        }
        else{
            return 1;
        }
    }
    # baseQ && mapQ && tail-distance
    if(    $alt_mBQ < $min_mBQ
        || $alt_mMQ < $min_mMQ
        || $alt_mTD < $min_mTD
    ){
        return 0;
    }
    else{
        return 1;
    }
}

1; ## tell the perl script the successful access of this module.
