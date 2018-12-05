#!/usr/bin/perl -w
use strict;
use BioFuse::Stat::ChiSquareTest::FourFoldTable;



my %value = (r1c1=>83, r1c2=>276, r2c1=>180, r2c2=>41);
# my @value = @ARGV;
# my @value = (12, 2, 10, 4);
# my @value = (12, 0, 10, 0);
#my @value = @ARGV;

my $FourFoldTable = BioFuse::Stat::ChiSquareTest::FourFoldTable->new(Value_Href=>\%value);
my $T = $FourFoldTable->get_theroy_value;
print "@$T\n";
print $FourFoldTable->get_method."\n";

my $P = $FourFoldTable -> get_P_value;
print "$P\n";

print $FourFoldTable->get_chi_square."\n";

print $FourFoldTable->get_odd_ratio."\n";

print $FourFoldTable->get_risk_ratio."\n";
