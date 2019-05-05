package BioFuse::Util::Log;

use strict;
use warnings;
use Carp qw/ cluck /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              warn_and_exit
              cluck_and_exit
              stout_and_sterr
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Util::Log';
#----- version --------
$VERSION = "0.32";
$DATE = '2019-05-05';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @functoion_list = qw/
                        warn_and_exit
                        cluck_and_exit
                        stout_and_sterr
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

1; ## tell the perl script the successful access of this module.
