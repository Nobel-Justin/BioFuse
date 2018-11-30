package BioFuse::Dist::Statistics::FisherTest;

use strict;
use warnings;
use List::Util qw/ sum /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              $E_Value
              fisher_P
              log_sum
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::Dist::Statistics::FisherTest';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-11-26';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        fisher_P
                        log_sum
                     /;

my $E_Value = 2.718281828;

#--- get four values Fisher P
sub fisher_P{
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
	my ($av,$bv,$cv,$dv) = @_;
	my $ab_sum = $av + $bv;
	my $ac_sum = $av + $cv;
	my $bd_sum = $bv + $dv;
	my $cd_sum = $cv + $dv;
	my $sum = sum($av,$bv,$cv,$dv);
	my $P_1 = &log_sum($ab_sum) + &log_sum($ac_sum) + &log_sum($bd_sum) + &log_sum($cd_sum);
	my $P_2 = &log_sum($av) + &log_sum($bv) + &log_sum($cv) + &log_sum($dv) + &log_sum($sum);
	my $P = $E_Value ** ($P_1 - $P_2);
	return $P;
}

#--- log sum for factorial
sub log_sum{
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
	return $_[0] < 1 ? 0 : sum( map {log($_)} (1 .. $_[0]) );
}

1; ## tell the perl script the successful access of this module.
