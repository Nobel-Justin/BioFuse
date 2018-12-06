#!/usr/bin/perl -w
use strict;
use BioFuse::Stat::ChiSquareTest::FourFoldTable;
use Data::Dumper;


# my %value = (r1c1=>83, r1c2=>276, r2c1=>180, r2c2=>41); # for x^2
# my %value = (r1c1=>992, r1c2=>2260, r2c1=>165, r2c2=>1017); # for OR
# my %value = (r1c1=>9, r1c2=>41, r2c1=>20, r2c2=>29); # for RR
my %value = (r1c1=>298, r1c2=>2757, r2c1=>81, r2c2=>663); # for RD
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

print Dumper($FourFoldTable->get_odds_ratio);

print Dumper($FourFoldTable->get_risk_ratio);

print Dumper($FourFoldTable->get_ratio_diff(aim=>'c1', wrt=>'r1'));

