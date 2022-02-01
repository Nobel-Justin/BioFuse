package BioFuse::Stat::ChiSquareTest::FourFoldTable;

use strict;
use warnings;
use List::Util qw/ min sum /;
use Statistics::Distributions;
use BioFuse::Util::Log qw/ warn_and_exit /;
use BioFuse::Stat::FisherTest qw/ fisher_P $E_Value /;
use BioFuse::Stat::ConfInt qw/ CIL2Z /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::Stat::ChiSquareTest::FourFoldTable';
#----- version --------
$VERSION = "0.03";
$DATE = '2018-12-06';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @function_list = qw/
                        new
                        get_odds_ratio
                        get_risk_ratio
                        get_ratio_diff
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
# ffTableOB -> OR -> $OR_id = {OR => $odds_ratio, CI => [$CI_lower, $CI_upper], CIL => $CIL, significant => 0/1}
# ffTableOB -> RR -> $RR_id = {RR => $risk_ratio, CI => [$CI_lower, $CI_upper], CIL => $CIL, significant => 0/1}
# ffTableOB -> RD -> $RD_id = {RD => $ratio_diff, CI => [$CI_lower, $CI_upper], CIL => $CIL, significant => 0/1}

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

#--- return odds ratio of given respect and aim
sub get_odds_ratio{
    my $ffTableOB = shift;
    my %parm = @_;
    my $aim = $parm{aim} || 'c1'; # r1/r2/c1/c2
    my $wrt = $parm{wrt} || 'r2'; # r1/r2/c1/c2 with repect to
    my $CIL = $parm{CIL} || 0.95;
    my $ratio_digit = $parm{ratio_digit} || 4;

    my ($aim_rc,$aim_no) = ($aim =~ /([rc])([12])/);
    my ($wrt_rc,$wrt_no) = ($wrt =~ /([rc])([12])/);
    if($aim_rc eq $wrt_rc){
        warn_and_exit "<ERROR>\twrong input for odds_ratio\n"
                            ."\taim=$aim, wrt=>$wrt\n";
    }
    # check id
    my $OR_id = "$aim$wrt";
    return $ffTableOB->{OR}->{$OR_id} if exists $ffTableOB->{OR}->{$OR_id} && $ffTableOB->{OR}->{$OR_id}->{CIL} == $CIL;
    # OR = (Yaim_Nwrt / Naim_Nwrt) / (Yaim_Ywrt / Naim_Ywrt)
    my ($Yaim_Nwrt, $Naim_Nwrt, $Yaim_Ywrt, $Naim_Ywrt);
    # fill number
    my $aim_ot = ($aim_no == 1 ? 2 : 1);
    my $wrt_ot = ($wrt_no == 1 ? 2 : 1);
    if($aim_rc eq 'r'){
        $Yaim_Nwrt = $ffTableOB->{real_value}->{$aim_no}->{$wrt_ot};
        $Naim_Nwrt = $ffTableOB->{real_value}->{$aim_ot}->{$wrt_ot};
        $Yaim_Ywrt = $ffTableOB->{real_value}->{$aim_no}->{$wrt_no};
        $Naim_Ywrt = $ffTableOB->{real_value}->{$aim_ot}->{$wrt_no};
    }
    else{
        $Yaim_Nwrt = $ffTableOB->{real_value}->{$wrt_ot}->{$aim_no};
        $Naim_Nwrt = $ffTableOB->{real_value}->{$wrt_ot}->{$aim_ot};
        $Yaim_Ywrt = $ffTableOB->{real_value}->{$wrt_no}->{$aim_no};
        $Naim_Ywrt = $ffTableOB->{real_value}->{$wrt_no}->{$aim_ot};
    }
    # Haldane-Anscombe correction
    if($Yaim_Nwrt * $Naim_Nwrt * $Yaim_Ywrt * $Naim_Ywrt == 0){
        $_ += 0.5 for ($Yaim_Nwrt, $Naim_Nwrt, $Yaim_Ywrt, $Naim_Ywrt);
    }
    # calculate odds ratio
    my $odds_ratio = ($Yaim_Nwrt / $Naim_Nwrt) / ($Yaim_Ywrt / $Naim_Ywrt);
    # calculate Confidence Interval
    my $z = CIL2Z(CIL => $CIL);
    my $ln_OR = log($odds_ratio);
    my $SE_ln_OR = sqrt(sum(map {(1/$_)} ($Yaim_Nwrt, $Naim_Nwrt, $Yaim_Ywrt, $Naim_Ywrt)));
    my $CI_lower = sprintf "%.${ratio_digit}f", $E_Value ** ($ln_OR - $z * $SE_ln_OR),
    my $CI_upper = sprintf "%.${ratio_digit}f", $E_Value ** ($ln_OR + $z * $SE_ln_OR);
    my $sigificant = ($CI_lower <= 1 && $CI_upper >= 1) ? 0 : 1;
    # record
    $odds_ratio = sprintf "%.${ratio_digit}f", $odds_ratio;
    $ffTableOB->{OR}->{$OR_id} = {OR => $odds_ratio, CI => [$CI_lower, $CI_upper], CIL => $CIL, sigificant => $sigificant};
    # return
    return $ffTableOB->{OR}->{$OR_id};
}

#--- return risk ratio of given respect and aim
sub get_risk_ratio{
    my $ffTableOB = shift;
    my %parm = @_;
    my $aim = $parm{aim} || 'c1'; # r1/r2/c1/c2
    my $wrt = $parm{wrt} || 'r2'; # r1/r2/c1/c2 with repect to
    my $CIL = $parm{CIL} || 0.95;
    my $ratio_digit = $parm{ratio_digit} || 4;

    my ($aim_rc,$aim_no) = ($aim =~ /([rc])([12])/);
    my ($wrt_rc,$wrt_no) = ($wrt =~ /([rc])([12])/);
    if($aim_rc eq $wrt_rc){
        warn_and_exit "<ERROR>\twrong input for odds_ratio\n"
                            ."\taim=$aim, wrt=>$wrt\n";
    }
    # check id
    my $RR_id = "$aim$wrt";
    return $ffTableOB->{RR}->{$RR_id} if exists $ffTableOB->{RR}->{$RR_id} && $ffTableOB->{RR}->{$RR_id}->{CIL} == $CIL;
    # RR = (Yaim_Nwrt / Sum_Nwrt) / (Yaim_Ywrt / Sum_Ywrt)
    my ($Yaim_Nwrt, $Sum_Nwrt, $Yaim_Ywrt, $Sum_Ywrt);
    # fill number
    my $aim_ot = ($aim_no == 1 ? 2 : 1);
    my $wrt_ot = ($wrt_no == 1 ? 2 : 1);
    if($aim_rc eq 'r'){
        $Yaim_Nwrt = $ffTableOB->{real_value}->{$aim_no}->{$wrt_ot};
        $Sum_Nwrt  = $ffTableOB->{col_sum}->{$wrt_ot};
        $Yaim_Ywrt = $ffTableOB->{real_value}->{$aim_no}->{$wrt_no};
        $Sum_Ywrt  = $ffTableOB->{col_sum}->{$wrt_no};
    }
    else{
        $Yaim_Nwrt = $ffTableOB->{real_value}->{$wrt_ot}->{$aim_no};
        $Sum_Nwrt  = $ffTableOB->{row_sum}->{$wrt_ot};
        $Yaim_Ywrt = $ffTableOB->{real_value}->{$wrt_no}->{$aim_no};
        $Sum_Ywrt  = $ffTableOB->{row_sum}->{$wrt_no};
    }
    # Haldane-Anscombe correction
    if($Yaim_Nwrt * $Sum_Nwrt * $Yaim_Ywrt * $Sum_Ywrt == 0){
        $_ += 0.5 for ($Yaim_Nwrt, $Yaim_Ywrt);
        $_ += 1   for ($Sum_Nwrt,  $Sum_Ywrt);
    }
    # calculate risk ratio
    my $risk_ratio = ($Yaim_Nwrt / $Sum_Nwrt) / ($Yaim_Ywrt / $Sum_Ywrt);
    # calculate Confidence Interval
    my $z = CIL2Z(CIL => $CIL);
    my $ln_RR = log($risk_ratio);
    my $SE_ln_RR = sqrt(($Sum_Nwrt - $Yaim_Nwrt) / $Yaim_Nwrt / $Sum_Nwrt + ($Sum_Ywrt - $Yaim_Ywrt) / $Yaim_Ywrt / $Sum_Ywrt);
    my $CI_lower = sprintf "%.${ratio_digit}f", $E_Value ** ($ln_RR - $z * $SE_ln_RR),
    my $CI_upper = sprintf "%.${ratio_digit}f", $E_Value ** ($ln_RR + $z * $SE_ln_RR);
    my $sigificant = ($CI_lower <= 1 && $CI_upper >= 1) ? 0 : 1;
    # record
    $risk_ratio = sprintf "%.${ratio_digit}f", $risk_ratio;
    $ffTableOB->{RR}->{$RR_id} = {RR => $risk_ratio, CI => [$CI_lower, $CI_upper], CIL => $CIL, sigificant => $sigificant};
    # return
    return $ffTableOB->{RR}->{$RR_id};
}

#--- return ratio difference of given respect and aim
sub get_ratio_diff{
    my $ffTableOB = shift;
    my %parm = @_;
    my $aim = $parm{aim} || 'c1'; # r1/r2/c1/c2
    my $wrt = $parm{wrt} || 'r2'; # r1/r2/c1/c2 with repect to
    my $CIL = $parm{CIL} || 0.95;
    my $ratio_digit = $parm{ratio_digit} || 4;

    my ($aim_rc,$aim_no) = ($aim =~ /([rc])([12])/);
    my ($wrt_rc,$wrt_no) = ($wrt =~ /([rc])([12])/);
    if($aim_rc eq $wrt_rc){
        warn_and_exit "<ERROR>\twrong input for odds_ratio\n"
                            ."\taim=$aim, wrt=>$wrt\n";
    }
    # check id
    my $RD_id = "$aim$wrt";
    return $ffTableOB->{RD}->{$RD_id} if exists $ffTableOB->{RD}->{$RD_id} && $ffTableOB->{RD}->{$RD_id}->{CIL} == $CIL;
    # RD = (Yaim_Nwrt / Sum_Nwrt) - (Yaim_Ywrt / Sum_Ywrt)
    my ($Yaim_Nwrt, $Sum_Nwrt, $Yaim_Ywrt, $Sum_Ywrt);
    # fill number
    my $aim_ot = ($aim_no == 1 ? 2 : 1);
    my $wrt_ot = ($wrt_no == 1 ? 2 : 1);
    if($aim_rc eq 'r'){
        $Yaim_Nwrt = $ffTableOB->{real_value}->{$aim_no}->{$wrt_ot};
        $Sum_Nwrt  = $ffTableOB->{col_sum}->{$wrt_ot};
        $Yaim_Ywrt = $ffTableOB->{real_value}->{$aim_no}->{$wrt_no};
        $Sum_Ywrt  = $ffTableOB->{col_sum}->{$wrt_no};
    }
    else{
        $Yaim_Nwrt = $ffTableOB->{real_value}->{$wrt_ot}->{$aim_no};
        $Sum_Nwrt  = $ffTableOB->{row_sum}->{$wrt_ot};
        $Yaim_Ywrt = $ffTableOB->{real_value}->{$wrt_no}->{$aim_no};
        $Sum_Ywrt  = $ffTableOB->{row_sum}->{$wrt_no};
    }
    # Haldane-Anscombe correction
    if($Sum_Nwrt * $Sum_Ywrt == 0){
        $_ += 0.5 for ($Yaim_Nwrt, $Yaim_Ywrt);
        $_ += 1   for ($Sum_Nwrt,  $Sum_Ywrt);
    }
    # calculate risk difference
    my $ratio_Nwrt = $Yaim_Nwrt / $Sum_Nwrt;
    my $ratio_Ywrt = $Yaim_Ywrt / $Sum_Ywrt;
    my $ratio_diff = $ratio_Nwrt - $ratio_Ywrt;
    # calculate Confidence Interval
    my $z = CIL2Z(CIL => $CIL);
    my $SE_RD = sqrt($ratio_Nwrt * (1 - $ratio_Nwrt) / $Sum_Nwrt + $ratio_Ywrt * (1 - $ratio_Ywrt) / $Sum_Ywrt);
    my $CI_lower = sprintf "%.${ratio_digit}f", $ratio_diff - $z * $SE_RD;
    my $CI_upper = sprintf "%.${ratio_digit}f", $ratio_diff + $z * $SE_RD;
    my $sigificant = ($CI_lower <= 0 && $CI_upper >= 0) ? 0 : 1;
    # record
    $ratio_diff = sprintf "%.${ratio_digit}f", $ratio_diff;
    $ffTableOB->{RD}->{$RD_id} = {RD => $ratio_diff, CI => [$CI_lower, $CI_upper], CIL => $CIL, sigificant => $sigificant};
    # return
    return $ffTableOB->{RD}->{$RD_id};
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
