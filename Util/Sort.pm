package BioFuse::Util::Sort;

use strict;
use warnings;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              sortByStrAndSubNum
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Util::Sort';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-11-01';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        sortByStrAndSubNum
                     /;

#--- sort with string and its sub-number ---
## such chr-pos
## ascending sort
sub sortByStrAndSubNum{
    # options
    shift if ($_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $Str_a = $parm{Str_a};
    my $Num_a = $parm{Num_a};
    my $Str_b = $parm{Str_b};
    my $Num_b = $parm{Num_b};
    my $SortHref = $parm{SortHref} || undef;
    my $SortKey  = $parm{SortKey}  || undef;

    return    (  (defined $SortHref && defined $SortKey)
               ? (($SortHref->{$Str_a}->{$SortKey} || 1E10) # when such mseg doesn't exist,
                  <=>
                  ($SortHref->{$Str_b}->{$SortKey} || 1E10) # use LARGE number as its turn
                 )
               : ( $Str_a
                  cmp
                   $Str_b
                 )
              )
           || $Num_a <=> $Num_b;
}

#--- 
1; ## tell the perl script the successful access of this module.
