package BioFuse::Util::Sys;

use strict;
use warnings;
use Cwd qw/ abs_path /;
use BioFuse::Util::Log qw/ warn_and_exit /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              file_exist
              trible_run_for_success
              check_java_version
              reset_folder
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Util::Sys';
#----- version --------
$VERSION = "0.33";
$DATE = '2019-05-07';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        file_exist
                        trible_run_for_success
                        check_java_version
                        reset_folder
                     /;

#--- verify existence of file/softlink
sub file_exist{
    shift @_ if(@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $filePath = $parm{filePath};
    my $alert = $parm{alert} || 0;

    if(    !defined $filePath
        || length($filePath) == 0
        || !-e abs_path( $filePath )
        || `file $filePath` =~ /broken symbolic link/
    ){
        warn_and_exit "<ERROR>\tFile does not exist:\n"
                            ."\t$filePath\n" if( $alert );
        return 0;
    }
    else{
        return 1;
    }
}

#--- run shell command trible times ---
# verbose_Href->{cmd_Nvb, esdo_Nvb, log_vb}
sub trible_run_for_success{
    shift if(@_ && $_[0] =~ /$MODULE_NAME/);
    my ($command,$type,$verbose_Href) = @_;

    $command = 'set -o pipefail; ' . $command;
    if(!defined($verbose_Href) || !($verbose_Href->{cmd_Nvb}||0)){
        warn "$command\n";
    }

    my $dev_null = '2>/dev/null';
    if(!defined($verbose_Href) || !($verbose_Href->{esdo_Nvb}||0)){
        $dev_null = '';
    }

    my $run_time = 0;
    RUN: {
        chomp(my $run_log = `( $command $dev_null ) && ( echo $type-ok )`);
        $run_time++; # record operation time
        if($run_log !~ /$type-ok$/){
            redo RUN if ($run_time < 3); # maximum is three times
        }
        else{
            if(defined($verbose_Href) && $verbose_Href->{log_vb}){
                $run_log =~ s/$type-ok$//;
                print "$run_log\n";
            }
            return 1;
        }
    }

    # if reach here, fail
    warn_and_exit "<ERROR>\t$type command fails three times.\n"
                        ."\t$command\n";
}

#--- check java version ---
sub check_java_version{
    shift @_ if(@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $javaPath = $parm{javaPath};
    my $minVer = $parm{minVer} || 1.8;

    my ($java_version) = (`$javaPath -Xmx1m -version 2>&1` =~ /\s+version\D+(\d+\.\d+)/);
    if($java_version < $minVer){
        warn_and_exit "<ERROR>\tThe version of java must be at least $minVer, but yours is $java_version\n";
    }
}

#--- try to `rm -rf` and then `mkdir -p` given folder path ---
sub reset_folder{
    shift @_ if(@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $folder = $parm{folder};
    `rm -rf $folder`;
    `mkdir -p $folder`;
}

1; ## tell the perl script the successful access of this module.
