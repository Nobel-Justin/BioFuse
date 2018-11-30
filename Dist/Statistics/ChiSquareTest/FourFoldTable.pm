package BioFuse::Dist::Statistics::ChiSquareTest::FourFoldTable;

use strict;
use warnings;
use List::Util qw/ min sum /;
use Statistics::Distributions;
use BioFuse::Util::Log qw/ warn_and_exit /;
use BioFuse::Dist::Statistics::FisherTest qw/ fisher_P  /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::Dist::Statistics::ChiSquareTest::FourFoldTable';
#----- version --------
$VERSION = "0.02";
$DATE = '2018-11-26';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        get_odd_ratio
                        get_risk_ratio
                        get_theroy_value
                        get_P_value
                        get_chi_square
                        get_method
                     /;

#--- structure of object
# ffTableOB -> real_value = { 1=>{1=>,2=>}, 2=>{1=>,2=>} }
# ffTableOB -> sum  = $sum
# ffTableOB -> row_sum = {1=>,2=>}
# ffTableOB -> col_sum = {1=>,2=>}
# ffTableOB -> theroy_value -> { 1=>{1=>,2=>}, 2=>{1=>,2=>} }
# ffTableOB -> P_cal -> {method=>, chi_square=>, value=>$P}

#--- construction of object
sub new{
	shift;
    my %parm = @_;
    my $Value_Href = $parm{Value_Href};
    my $method = $parm{method} || undef; # basal_fomula or adjusted_fomula or fisher_test

	my $ffTableOB = {};
	# real values
	$ffTableOB->{real_value}->{1}->{1} = $Value_Href->{r1c1}; # 1st row, 1st column
	$ffTableOB->{real_value}->{1}->{2} = $Value_Href->{r1c2}; # 1st row, 2nd column
	$ffTableOB->{real_value}->{2}->{1} = $Value_Href->{r2c1}; # 2nd row, 1st column
	$ffTableOB->{real_value}->{2}->{2} = $Value_Href->{r2c2}; # 2nd row, 2nd column
	# sum of all values
	$ffTableOB->{sum} = sum(values %$Value_Href);
	# sum of column and row
	$ffTableOB->{row_sum}->{1} = $ffTableOB->{real_value}->{1}->{1} + $ffTableOB->{real_value}->{1}->{2};
	$ffTableOB->{row_sum}->{2} = $ffTableOB->{real_value}->{2}->{1} + $ffTableOB->{real_value}->{2}->{2};
	$ffTableOB->{col_sum}->{1} = $ffTableOB->{real_value}->{1}->{1} + $ffTableOB->{real_value}->{2}->{1};
	$ffTableOB->{col_sum}->{2} = $ffTableOB->{real_value}->{1}->{2} + $ffTableOB->{real_value}->{2}->{2};
	# theroy values
	my @theroy_value;
	for my $row_NO (1 , 2){
		for my $col_NO (1 , 2){
			$ffTableOB->{theroy_value}->{$row_NO}->{$col_NO} = ($ffTableOB->{row_sum}->{$row_NO} * $ffTableOB->{col_sum}->{$col_NO}) / $ffTableOB->{sum};
			push @theroy_value , $ffTableOB->{theroy_value}->{$row_NO}->{$col_NO};
		}
	}
	# judge the P value calculation
	if(    defined $method
		&& $method =~ /^basal_fomula|adjusted_fomula|fisher_test$/
	){
		$ffTableOB->{P_cal}->{method} = $method;
	}
	elsif($ffTableOB->{sum} >= 40 && min(@theroy_value) >= 5){ # normal x2-test
		$ffTableOB->{P_cal}->{method} = 'basal_fomula';
	}
	elsif($ffTableOB->{sum} >= 40 && min(@theroy_value) < 5 && min(@theroy_value) >= 1){
		$ffTableOB->{P_cal}->{method} = 'adjusted_fomula';
	}
	else{
		$ffTableOB->{P_cal}->{method} = 'fisher_test';
	}
	bless($ffTableOB);
	return $ffTableOB;
}

#--- return odd ratio of given respect and aim
sub get_odd_ratio{
	my $ffTableOB = shift;
    my %parm = @_;
    my $aim = $parm{aim} || 'c1'; # r1/r2/c1/c2
    my $wrt = $parm{wrt} || 'r2'; # r1/r2/c1/c2 with repect to

    my ($aim_rc,$aim_no) = ($aim =~ /([rc])([12])/);
    my ($wrt_rc,$wrt_no) = ($wrt =~ /([rc])([12])/);
    if($aim_rc eq $wrt_rc){
    	warn_and_exit "<ERROR>\twrong input for odd_ratio\n"
    						."\taim=$aim, wrt=>$wrt\n";
    }
    # calculate odd ratio
    my $aim_ot = ($aim_no == 1 ? 2 : 1);
    my $wrt_ot = ($wrt_no == 1 ? 2 : 1);
    if($aim_rc eq 'r'){
    	return    ($ffTableOB->{real_value}->{$aim_no}->{$wrt_ot} * $ffTableOB->{real_value}->{$aim_ot}->{$wrt_no})
    		   / (($ffTableOB->{real_value}->{$aim_no}->{$wrt_no} * $ffTableOB->{real_value}->{$aim_ot}->{$wrt_ot}) || 1);
    }
    else{
    	return    ($ffTableOB->{real_value}->{$wrt_ot}->{$aim_no} * $ffTableOB->{real_value}->{$wrt_no}->{$aim_ot})
    		   / (($ffTableOB->{real_value}->{$wrt_no}->{$aim_no} * $ffTableOB->{real_value}->{$wrt_ot}->{$aim_ot}) || 1);
    }
}

#--- return risk ratio of given respect and aim
sub get_risk_ratio{
	my $ffTableOB = shift;
    my %parm = @_;
    my $aim = $parm{aim} || 'c1'; # r1/r2/c1/c2
    my $wrt = $parm{wrt} || 'r2'; # r1/r2/c1/c2 with repect to

    my ($aim_rc,$aim_no) = ($aim =~ /([rc])([12])/);
    my ($wrt_rc,$wrt_no) = ($wrt =~ /([rc])([12])/);
    if($aim_rc eq $wrt_rc){
    	warn_and_exit "<ERROR>\twrong input for odd_ratio\n"
    						."\taim=$aim, wrt=>$wrt\n";
    }
    # calculate risk ratio
    my $aim_ot = ($aim_no == 1 ? 2 : 1);
    my $wrt_ot = ($wrt_no == 1 ? 2 : 1);
    if($aim_rc eq 'r'){
    	return    ($ffTableOB->{real_value}->{$aim_no}->{$wrt_ot} * $ffTableOB->{col_sum}->{$wrt_no})
    		   / (($ffTableOB->{real_value}->{$aim_no}->{$wrt_no} * $ffTableOB->{col_sum}->{$wrt_ot}) || 1);
    }
    else{
    	return    ($ffTableOB->{real_value}->{$wrt_ot}->{$aim_no} * $ffTableOB->{row_sum}->{$wrt_no})
    		   / (($ffTableOB->{real_value}->{$wrt_no}->{$aim_no} * $ffTableOB->{row_sum}->{$wrt_ot}) || 1);
    }
}

#--- return the theroy values Aref
sub get_theroy_value{ # a, b, c, d
	my $ffTableOB = shift;
	my @theroy_value;
	for my $row_NO (1 , 2){
		for my $col_NO (1 , 2){
			push @theroy_value , $ffTableOB->{theroy_value}->{$row_NO}->{$col_NO};
		}
	}
	return \@theroy_value;
}

#--- return x^2-test P value
sub get_P_value{
	my $ffTableOB = shift;
	# calculated before
	return $ffTableOB->{P_cal}->{value} if defined $ffTableOB->{P_cal}->{value};
	# calculate
	my $method = $ffTableOB->{P_cal}->{method};
	if($method =~ /fomula/){ # basal_fomula or adjusted_fomula
		my $adjust = ($method eq 'basal_fomula')?0:0.5;
		my $chi_square = 0;
		for my $row_NO (1 , 2){
			for my $col_NO (1 , 2){
				my $real_value = $ffTableOB->{real_value}->{$row_NO}->{$col_NO};
				my $theroy_value = $ffTableOB->{theroy_value}->{$row_NO}->{$col_NO};
				$chi_square += (abs($real_value-$theroy_value) - $adjust) ** 2 / $theroy_value;
			}
		}
		$ffTableOB->{P_cal}->{chi_square} = $chi_square;
		$ffTableOB->{P_cal}->{value} = Statistics::Distributions::chisqrprob (1,$chi_square);
		return $ffTableOB->{P_cal}->{value}; # freedom degree is (2-1)*(2-1)
	}
	else{ # fisher test, two sides
		## original P value
		my $av = $ffTableOB->{real_value}->{1}->{1}; # 1st row, 1st column
		my $bv = $ffTableOB->{real_value}->{1}->{2}; # 1st row, 2nd column
		my $cv = $ffTableOB->{real_value}->{2}->{1}; # 2nd row, 1st column
		my $dv = $ffTableOB->{real_value}->{2}->{2}; # 2nd row, 2nd column
		my $orig_P = fisher_P($av,$bv,$cv,$dv);
		## get the minimum grid, row and col
		my $min_row_NO = ($ffTableOB->{row_sum}->{1}<$ffTableOB->{row_sum}->{2})?1:2;
		my $min_col_NO = ($ffTableOB->{col_sum}->{1}<$ffTableOB->{col_sum}->{2})?1:2;
		my $min_sum = min($ffTableOB->{row_sum}->{$min_row_NO} , $ffTableOB->{col_sum}->{$min_col_NO});
		## all combination for Fisher test
		my $Fisher_P = 0;
		for (0 .. $min_sum){
			$ffTableOB->{real_value}->{$min_row_NO}->{$min_col_NO} = $_;
			$ffTableOB->{real_value}->{$min_row_NO}->{$min_col_NO%2+1} = $ffTableOB->{row_sum}->{$min_row_NO} - $_;
			$ffTableOB->{real_value}->{$min_row_NO%2+1}->{$min_col_NO} = $ffTableOB->{col_sum}->{$min_col_NO} - $_;
			$ffTableOB->{real_value}->{$min_row_NO%2+1}->{$min_col_NO%2+1} = $ffTableOB->{col_sum}->{$min_col_NO%2+1} - $ffTableOB->{real_value}->{$min_row_NO}->{$min_col_NO%2+1};
			my $av = $ffTableOB->{real_value}->{1}->{1}; # 1st row, 1st column
			my $bv = $ffTableOB->{real_value}->{1}->{2}; # 1st row, 2nd column
			my $cv = $ffTableOB->{real_value}->{2}->{1}; # 2nd row, 1st column
			my $dv = $ffTableOB->{real_value}->{2}->{2}; # 2nd row, 2nd column
			my $Pv = fisher_P($av,$bv,$cv,$dv);
			$Fisher_P += $Pv if($Pv <= $orig_P);
		}
		## reset
		$ffTableOB->{real_value}->{1}->{1} = $av; # 1st row, 1st column
		$ffTableOB->{real_value}->{1}->{2} = $bv; # 1st row, 2nd column
		$ffTableOB->{real_value}->{2}->{1} = $cv; # 2nd row, 1st column
		$ffTableOB->{real_value}->{2}->{2} = $dv; # 2nd row, 2nd column
		$ffTableOB->{P_cal}->{chi_square} = 'NA_Fisher';
		$ffTableOB->{P_cal}->{value} = $Fisher_P;
		return $ffTableOB->{P_cal}->{value}; # freedom degree is (2-1)*(2-1)
	}
}

#--- return chi_square value 
sub get_chi_square{
	my $ffTableOB = shift;
	return $ffTableOB->{P_cal}->{chi_square};
}

#--- return test method
sub get_method{
	my $ffTableOB = shift;
	return $ffTableOB->{P_cal}->{method};
}

1; ## tell the perl script the successful access of this module.
