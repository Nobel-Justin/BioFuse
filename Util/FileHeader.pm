package BioFuse::Util::FileHeader;

use strict;
use warnings;
use BioFuse::Util::GZfile qw/ Try_GZ_Read /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              getHeaderTag
              lineInfoToHash
              read_headed_list
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Util::FileHeader';
#----- version --------
$VERSION = "0.02";
$DATE = '2019-11-25';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @function_list = qw/
                        getHeaderTag
                        lineInfoToHash
                        read_headed_list
                     /;

#--- get tags in file header (prefixed by '#') ---
## use \s+ to split
sub getHeaderTag{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
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
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $headTagAf = $parm{headTagAf};
    my $lineInfo = $parm{lineInfo};

    my @info = split /\s+/, $lineInfo;
    my %info = map{ ($headTagAf->[$_],$info[$_]) } (0 .. scalar(@$headTagAf)-1);
    return \%info;
}

#--- load headed list and do ---
## the list must have a header
sub read_headed_list{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $list = $parm{list};
    my $subrtRef = $parm{subrtRef};
    my $subrtParmAref = $parm{subrtParmAref};
    my $tab_div = $parm{tab_div} || 0;

    open (LIST,Try_GZ_Read($list)) || die "fail read list: $!\n";
    my @tag = map {s/^#//; ($_)} split /\s+/, lc(<LIST>);
    while(<LIST>){
        next if(/^#/);
        chomp;
        my @info = $tab_div ? split /\t+/ : split /\s+/;
        my %tagmap = map{ ($tag[$_], $info[$_]) } (0 .. $#info);
        # run sub-routine
        if (defined $subrtRef){
            my @parm = (tagmap => \%tagmap);
            push @parm, @$subrtParmAref if (defined $subrtParmAref);
            &{$subrtRef}(@parm);
        }
    }
    close LIST;
}

#--- 
1; ## tell the perl script the successful access of this module.
