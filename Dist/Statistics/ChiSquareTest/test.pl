#!/usr/bin/perl -w
use strict;
use BioFuse::Dist::Statistics::ChiSquareTest::FourFoldTable;



my %value = (r1c1=>99, r1c2=>15, r2c1=>750, r2c2=>21);
# my @value = @ARGV;
# my @value = (12, 2, 10, 4);
# my @value = (12, 0, 10, 0);
#my @value = @ARGV;

my $FourFoldTable = BioFuse::Dist::Statistics::ChiSquareTest::FourFoldTable->new(Value_Href=>\%value, method=>'basal_fomula');
my $T = $FourFoldTable->get_theroy_value;
print "@$T\n";
print $FourFoldTable->get_method."\n";

my $P = $FourFoldTable -> get_P_value;
print "$P\n";

print $FourFoldTable->get_chi_square."\n";

print $FourFoldTable->get_odd_ratio."\n";

print $FourFoldTable->get_risk_ratio."\n";
