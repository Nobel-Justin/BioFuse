package BioFuse::Dist::Statistics::SpearmanCorr;

use strict;
use warnings;
use List::Util qw/ min max sum /;
use BioFuse::Util::Log qw/ warn_and_exit /;
use BioFuse::Dist::DistStat qw/ get_value_mean /;
use BioFuse::Dist::Statistics::PearsonCorr qw/ get_pearson_corr /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
			  get_spearman_corr
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Dist::Statistics::SpearmanCorr';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-05-25';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
						get_spearman_corr
					 /;

#--- return spearman correlation coefficient ---
sub get_spearman_corr{

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

	# sort a_A, assign orig order
	my %a_v2i;
	my $a_o = 0;
	my %a_i2o = map { my $value = $list_a_Aref->[$_];
					  push @{$a_v2i{$value}}, $_;
					  ( $_, { orig_o => ++$a_o } )
					}
				sort { $list_a_Aref->[$a] <=> $list_a_Aref->[$b] }
				( 0 .. scalar(@$list_a_Aref)-1 );
	# assign norm order
	for my $i_Aref ( values %a_v2i){
		if( scalar(@$i_Aref) > 1 ){
			my @orig_o = map { $a_i2o{$_}->{orig_o} } @$i_Aref;
			my $mean_o = get_value_mean( value_Aref => \@orig_o );
			$a_i2o{$_}->{norm_o} = $mean_o for @$i_Aref;
		}
		else{
			$a_i2o{$_}->{norm_o} = $a_i2o{$_}->{orig_o} for @$i_Aref;
		}
	}

	# sort b_A, assign orig order
	my %b_v2i;
	my $b_o = 0;
	my %b_i2o = map { my $value = $list_b_Aref->[$_];
					  push @{$b_v2i{$value}}, $_;
					  ( $_, { orig_o => ++$b_o } )
					}
				sort { $list_b_Aref->[$a] <=> $list_b_Aref->[$b] }
				( 0 .. scalar(@$list_b_Aref)-1 );
	# assign norm order
	for my $i_Aref (values %b_v2i){
		if( scalar(@$i_Aref) > 1 ){
			my @orig_o = map { $b_i2o{$_}->{orig_o} } @$i_Aref;
			my $mean_o = get_value_mean( value_Aref => \@orig_o );
			$b_i2o{$_}->{norm_o} = $mean_o for @$i_Aref;
		}
		else{
			$b_i2o{$_}->{norm_o} = $b_i2o{$_}->{orig_o} for @$i_Aref;
		}
	}

	# prepare order array for PearsonCorr
	my @a_orders = map { $a_i2o{$_}->{norm_o} } sort {$a<=>$b} keys %a_i2o;
	my @b_orders = map { $b_i2o{$_}->{norm_o} } sort {$a<=>$b} keys %b_i2o;
	return get_pearson_corr(
							 list_a_Aref => \@a_orders,
							 list_b_Aref => \@b_orders,
							 ratio_digit => $ratio_digit,
							 error_getNA => $error_getNA
							);
}

#---
1; ## tell the perl script the successful access of this module.
