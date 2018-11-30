package BioFuse::Util::FileHeader;

use strict;
use warnings;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              getHeaderTag
              lineInfoToHash
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Util::FileHeader';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-11-24';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        getHeaderTag
                        lineInfoToHash
                     /;

#--- get tags in file header (prefixed by '#') ---
## use \s+ to split
sub getHeaderTag{
    # options
    shift if ($_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $header = $parm{header};

    $header =~ s/^#+//;
    my @head_tag = split /\s+/, lc($header);
    return \@head_tag;
}

#--- convert line info to headTag dict ---
## use \s+ to split
sub lineInfoToHash{
    # options
    shift if ($_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $headTagAf = $parm{headTagAf};
    my $lineInfo = $parm{lineInfo};

    my @info = split /\s+/, $lineInfo;
    my %info = map{ ($headTagAf->[$_],$info[$_]) } (0 .. scalar(@$headTagAf)-1);
    return \%info;
}

#--- 
1; ## tell the perl script the successful access of this module.
