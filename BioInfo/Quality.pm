package BioFuse::BioInfo::Quality;

use strict;
use warnings;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              baseQ_char2score
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Quality';
#----- version --------
$VERSION = "0.31";
$DATE = '2018-10-30';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @functoion_list = qw/
                        baseQ_char2score
                     /;

#--- convert quality char to score ---
sub baseQ_char2score{
    shift @_ if(@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $Q_char = $parm{Q_char};
    # set 33 for Sanger and Illumina 1.8+.
    # set 64 for Solexa, Illumina 1.3+ and 1.5+.
    my $Q_offset = $parm{Q_offset} || 33;

    return ( ord($Q_char) - $Q_offset );
}

1; ## tell the perl script the successful access of this module.
