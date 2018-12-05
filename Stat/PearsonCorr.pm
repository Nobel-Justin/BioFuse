package BioFuse::Stat::PearsonCorr;

use strict;
use warnings;
use List::Util qw/ min max sum /;
use BioFuse::Util::Log qw/ warn_and_exit /;
use BioFuse::Dist::DistStat qw/ get_value_mean /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  get_pearson_corr
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Stat::PearsonCorr';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-05-01';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						get_pearson_corr
					 /;

#--- return pearson correlation coefficient ---
sub get_pearson_corr{

	# options
	shift if (@_ && $_[0] =~ /$MODULE_NAME/);
	my %parm = @_;
	my $list_a_Aref = $parm{list_a_Aref} || [];
	my $list_b_Aref = $parm{list_b_Aref} || [];
	my $ratio_digit = $parm{ratio_digit} || 4;
	my $error_getNA = $parm{error_getNA} || 0;

	my $func_id = (caller(0))[3];
	if( scalar(@$list_a_Aref) <= 1 ){
		warn "<ERROR>\tList must have at least two elements for $func_id func.\n";
		if( $error_getNA ){ return 'N/A'; } else{ exit(1); }
	}
	if( scalar(@$list_a_Aref) != scalar(@$list_b_Aref) ){
		warn "<ERROR>\tLists have unequal number of elements in $func_id func.\n";
		if( $error_getNA ){ return 'N/A'; } else{ exit(1); }
	}

	# mean
	my $list_a_mean = get_value_mean( value_Aref => $list_a_Aref );
	my $list_b_mean = get_value_mean( value_Aref => $list_b_Aref );
	# Deviation
	my @list_a_dev = map{ $_ - $list_a_mean } @$list_a_Aref;
	my @list_b_dev = map{ $_ - $list_b_mean } @$list_b_Aref;
	# sum of mulplied diff-mean
	my $sum_ProdDev = 0;
	$sum_ProdDev += $list_a_dev[$_] * $list_b_dev[$_] for (0 .. $#list_a_dev);
	# sum of diff-mean square
	my $list_a_SumSQdev = 0;
	$list_a_SumSQdev += $list_a_dev[$_] ** 2 for (0 .. $#list_a_dev);
	my $list_b_SumSQdev = 0;
	$list_b_SumSQdev += $list_b_dev[$_] ** 2 for (0 .. $#list_b_dev);

	if( $list_a_SumSQdev * $list_b_SumSQdev == 0 ){
		warn "<ERROR>\tsum_sq_mdiff is zero in $func_id func.\n";
		if( $error_getNA ){ return 'N/A'; } else{ exit(1); }
	}

	my $pearson_corr = $sum_ProdDev / ( sqrt($list_a_SumSQdev) * sqrt($list_b_SumSQdev) );
	return sprintf "%.${ratio_digit}f", $pearson_corr;
}

#---
1; ## tell the perl script the successful access of this module.
