package BioFuse::Util::Random;

use strict;
use warnings;
use BioFuse::Util::Log qw/ stout_and_sterr warn_and_exit /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              GetRandBool
              PickRandAele
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Util::Random';
#----- version --------
$VERSION = "0.31";
$DATE = '2018-10-30';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @function_list = qw/
                        GetRandBool
                        PickRandAele
                     /;

#--- get a random probability ---
sub GetRandBool{
    shift @_ if(@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $Aref = $parm{Aref};
    my $RandCmp = $parm{RandCmp} || 0.5;
    my $verbose = $parm{verbose} || 0;

    $RandCmp = &PickRandAele( Aref => $Aref ) if( defined $Aref );
    my $RandProb = rand(1);
    my $RandBool = ( ($RandProb < $RandCmp) || 0 );

    if( $verbose ){
        $RandProb = sprintf "%.3f", $RandProb;
        stout_and_sterr "[RAND]\tRandProb: $RandProb, RandCmp: $RandCmp, RandBool: $RandBool\n";
    }

    return $RandBool;
}

#--- randomly pick one element from given Array (ref) ---
sub PickRandAele{
    shift @_ if(@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $Aref = $parm{Aref};

    my $Acnt = scalar(@$Aref);
    if( $Acnt == 0 ){
        warn_and_exit "<ERROR>\tempty array in PickRandAele func.\n";
    }
    elsif( $Acnt == 1 ){
        return $Aref->[0];
    }
    else{
        return $Aref->[int(rand($Acnt))];
    }
}

1; ## tell the perl script the successful access of this module.
