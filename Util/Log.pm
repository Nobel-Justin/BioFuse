package BioFuse::Util::Log;

use strict;
use warnings;
use Data::Dumper;
use Carp qw/ cluck /;
use List::Util qw/ max /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              warn_and_exit
              cluck_and_exit
              stout_and_sterr
              alignDisplay
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Util::Log';
#----- version --------
$VERSION = "0.34";
$DATE = '2021-12-28';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        warn_and_exit
                        cluck_and_exit
                        stout_and_sterr
                        alignDisplay
                     /;

#----------- warn out the content and exit -----------
sub warn_and_exit{
    my ($warn_content, $exit_signal) = @_;
    $exit_signal = 1 unless defined($exit_signal);
    warn "$warn_content";
    exit($exit_signal);
}

#----------- cluck the content and exit -----------
sub cluck_and_exit{
    my ($warn_content, $exit_signal) = @_;
    $exit_signal = 1 unless defined($exit_signal);
    cluck "$warn_content";
    exit($exit_signal);
}

#--- warn out the content and exit ---
sub stout_and_sterr{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my $content = shift @_;
    my %parm = @_;
    my $stderr = $parm{stderr} || 0;

    $| = 1; # no buffer
    # default STDOUT sololy
    print STDOUT "$content";
    # STDERR additionally, when stated
    warn "$content" if( $stderr );
}

#--- align display of given content ---
## contentAf = [ [row1col1,row1col2,..], [row2col1,row2col2,..], .. ]
sub alignDisplay{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $contentAf = $parm{contentAf};
    my $prefix = $parm{prefix} || '';
    my $colGap = $parm{colGap} || 3;
    my $fillEmpty = $parm{fillEmpty};
    my $alignTo = $parm{alignTo} || 'H'; # Head/Tail

    # width of each column
    my @colWidth;
    for my $rowAf (@$contentAf){
        $colWidth[$_] =   defined $colWidth[$_]
                        ? max($colWidth[$_], length($rowAf->[$_]))
                        : length($rowAf->[$_])
                        for 0 .. scalar(@$rowAf)-1;
    }
    # add column gap
    $_+=$colGap for @colWidth;
    $colWidth[0] -= $colGap if $alignTo eq 'T';
    # display
    my @display;
    my $alignSign = $alignTo eq 'T' ? '' : '-';
    for my $rowAf (@$contentAf){
        my $rowStr = $prefix;
        for my $idx (0 .. $#colWidth){
            # check
            my $value;
            if(defined $rowAf->[$idx]){
                $value = $rowAf->[$idx];
            }
            else{
                &cluck_and_exit ("<ERROR>\tlack ".($idx+1)." column in row:\n".Dumper($rowAf)) if !defined $fillEmpty;
                $value = $fillEmpty;
            }
            
            $rowStr .= sprintf("%$alignSign$colWidth[$idx]s", $value);
        }
        push @display, $rowStr;
    }
    # return
    return join("\n", @display)."\n";
}

1; ## tell the perl script the successful access of this module.
