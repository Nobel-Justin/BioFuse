package BioFuse::Stat::MultiTest;

use strict;
use warnings;
use Data::Dumper;
use List::Util qw/ min sum first /;
use BioFuse::Util::Log qw/ stout_and_sterr warn_and_exit /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::Stat::MultiTest';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-12-05';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        load_test
                        get_all_test
                        get_FDR_sig_test
                        get_p_adjust_method
                        p_adjust
                        set_P_sig_under_FDR
                     /;

#--- structure of object
# multi_test -> {test} = [ {ID=>, Pvalue=>, Qvalue=>, significant=>} ]
# multi_test -> {p_adjust_method} = BH/BY/BF
# multi_test -> {FDR} = $FDR

#--- construction of object ---
sub new{
    my $type = shift;
    my %parm = @_;

    my $multi_test = {};
    ## basic attributes
    $multi_test->{test} = [];
    $multi_test->{p_adjust_method} = 'none';
    $multi_test->{FDR} = undef;

    bless($multi_test);
    return $multi_test;
}

#--- load one test ---
sub load_test{
    my $multi_test = shift;
    my %parm = @_;
    my $testID = $parm{testID};
    my $Pvalue = $parm{Pvalue};
    my $testOB = $parm{testOB}; # optional

    # check
    if($multi_test->get_p_adjust_method ne 'none'){
        $multi_test->{p_adjust_method} = 'none'; # reset
        $multi_test->{FDR} = undef;
        stout_and_sterr "<WARN>\tReset pre-set p.adjust method when load test.\n";
    }

    # record
    push @{$multi_test->{test}}, { ID => $testID,
                                   Pvalue => $Pvalue,
                                   Qvalue => undef,
                                   significant => undef
                                 };
    # specific for 'testOB'
    $multi_test->{test}->[-1]->{testOB} = $testOB if defined $testOB;
}

#--- return Aref of all test ---
sub get_all_test{
    my $multi_test = shift;
    return $multi_test->{test};
}

#--- return Aref of significant test ---
sub get_FDR_sig_test{
    my $multi_test = shift;
    return [ map {($_)}
             grep {defined $_->{significant} && $_->{significant} == 1}
             @{$multi_test->{test}}
           ];
}

#--- return method of p.adjust ---
sub get_p_adjust_method{
    my $multi_test = shift;
    return $multi_test->{p_adjust_method};
}

#--- do p.adjust ---
## BH: Benjamini–Hochberg procedure (default)
## BY: Benjamini–Hochberg–Yekutieli procedure
## BF: Bonferroni procedure
sub p_adjust{
    my $multi_test = shift;
    my %parm = @_;
    my $method = $parm{method} || 'BH';
    my $orig_Q = $parm{orig_Q} || 0; # to disable cummin(Q) in 'BH' and 'BY'
    my $no_arb = $parm{no_arb} || 0; # to declear not arbitrary dependence in 'BY'

    my $testAf = $multi_test->{test};
    my $n = scalar @$testAf;
    if(    $method eq 'BH'
        || $method eq 'BY'
    ){
        # decreasing sort by P
        @$testAf = sort {$b->{Pvalue} <=> $a->{Pvalue}} @$testAf;
        # Q = P * n * c(m) / i;
        ## BH: c(m) = 1
        ## BY: c(m) = 1 when not arbitrary dependence
        my $i = $n;
        my $min = 1;
        my $cm = ($method eq 'BH' || $no_arb) ? 1 : sum( map {(1/$_)} (1 .. $n) );
        for (@$testAf){
            $_->{Qvalue} = min(1, $min, $_->{Pvalue} * $n * $cm / $i--);
            $min = min($min, $_->{Qvalue}) unless $orig_Q;
        }
        # increasing sort by Q
        @$testAf = sort {$a->{Pvalue} <=> $b->{Pvalue}} @$testAf;
    }
    elsif($method eq 'BF'){
        # Q = P * n
        $_->{Qvalue} = min(1, $n * $_->{Pvalue}) for @$testAf;
    }
    else{
        warn_and_exit "<ERROR>\tcannot recognize p_adjust method: $method\n";
    }

    # record method
    $multi_test->{p_adjust_method} = $method;
}

#--- set FDR and check P significant ---
sub set_P_sig_under_FDR{
    my $multi_test = shift;
    my %parm = @_;
    my $FDR = $parm{FDR} || 0.1;

    # check Q method
    if($multi_test->get_p_adjust_method eq 'none'){
        warn_and_exit "<ERROR>\tP values have not been adjusted.\n";
    }

    my $testAf = $multi_test->{test};
    my $n = scalar @$testAf;
    # increasing sort by P
    @$testAf = sort {$a->{Pvalue} <=> $b->{Pvalue}} @$testAf;
    # find max_i that has Q < FDR
    my $max_i = first {$testAf->[$_]->{Qvalue} <= $FDR} reverse (0 .. $n-1);
       $max_i = -1 unless defined $max_i;
    # set significant
    $testAf->[$_]->{significant} = 1 for (0 .. $max_i);
    $testAf->[$_]->{significant} = 0 for ($max_i+1 .. $n-1);

    # record FDR
    $multi_test->{FDR} = $FDR;
}

1; ## tell the perl script the successful access of this module.
