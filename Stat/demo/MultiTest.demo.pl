use BioFuse::Stat::MultiTest;
use Data::Dumper;

my $mt = BioFuse::Stat::MultiTest->new;

# $mt->load_test(testID=>1, Pvalue=>0.005, testOB=>9);
# $mt->load_test(testID=>2, Pvalue=>0.009);
# $mt->load_test(testID=>3, Pvalue=>0.05);
# $mt->load_test(testID=>4, Pvalue=>0.1);
# $mt->load_test(testID=>5, Pvalue=>0.2);
# $mt->load_test(testID=>6, Pvalue=>0.3);


my @a = (
	['Total_calories', 0.0006 ],
	['Olive_oil', 0.008 ],
	['Whole_milk', 0.039 ],
	['White_meat', 0.041 ],
	['Proteins', 0.042 ],
	['Nuts', 0.06 ],
	['Cereals_and_pasta', 0.074 ],
	['White_fish', 0.205 ],
	['Butter', 0.212 ],
	['Vegetables', 0.216 ],
	['Skimmed_milk', 0.222 ],
	['Red_meat', 0.251 ],
	['Fruit', 0.269 ],
	['Eggs', 0.275 ],
	['Blue_fish', 0.34 ],
	['Legumes', 0.341 ],
	['Carbohydrates', 0.384 ],
	['Potatoes', 0.569 ],
	['Bread', 0.594 ],
	['Fats', 0.696 ],
	['Sweets', 0.762 ],
	['Dairy_products', 0.94 ],
	['Semi-skimmed_milk', 0.942 ],
	['Total_meat', 0.975 ],
	['Processed_meat', 0.986 ]
);

$mt->load_test(testID=>$_->[0], Pvalue=>$_->[1]) for @a;


$mt->p_adjust(method => 'BH', orig_Q => 0);
# $mt->p_adjust(method => 'BH', orig_Q => 1);
# $mt->p_adjust(method => 'BY', orig_Q => 0, no_arb => 0);
# $mt->p_adjust(method => 'BY', orig_Q => 0, no_arb => 1);
# $mt->p_adjust(method => 'BY', orig_Q => 1);
# $mt->p_adjust(method => 'BF');

$mt->set_P_sig_under_FDR(FDR=>0.25);

# print Dumper($mt->get_all_test);
print Dumper($mt->get_FDR_sig_test);
