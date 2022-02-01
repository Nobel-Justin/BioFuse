package BioFuse::Stat::ConfInt;

use strict;
use warnings;
use BioFuse::Util::Log qw/ warn_and_exit /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              CIL2Z
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::Stat::ConfInt';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-12-06';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--- functions in this pm ---#
my @function_list = qw/
                        CIL2Z
                     /;

#--- variables ---#
my %TwoSideCI2Z = ( 0.99 => 2.576,
                    0.98 => 2.326,
                    0.95 => 1.96,
                    0.90 => 1.645
                  );

#--- get z matched to given Confidence Interval level ---
sub CIL2Z{
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $CIL = $parm{CIL} || 0.95;

    unless(exists $TwoSideCI2Z{$CIL}){
        warn_and_exit "<ERROR>\tcannot find z matched to given Confidence Interval level ($CIL).\n";
    }
    return $TwoSideCI2Z{$CIL};
}

1; ## tell the perl script the successful access of this module.
