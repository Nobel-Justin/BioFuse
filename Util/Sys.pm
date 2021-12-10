package BioFuse::Util::Sys;

use strict;
use warnings;
use Cwd qw/ abs_path /;
use Data::Dumper;
use BioFuse::Util::Log qw/ cluck_and_exit /;
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
              remove_folder
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Util::Sys';
#----- version --------
$VERSION = "0.37";
$DATE = '2021-12-07';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        file_exist
                        trible_run_for_success
                        check_java_version
                        reset_folder
                        remove_folder
                     /;

#--- verify existence of file/softlink
sub file_exist{
    shift @_ if(@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $filePath = $parm{filePath}; # scalar(filepath) or ref of hash/array/scalar
    my $alert = $parm{alert} || 0;

    my $exist = 1; # initial
    if( !defined $filePath ){
        $exist = 0;
    }
    elsif( ref($filePath) ){ # ref: ARRAY HASH SCALAR
        if( ref($filePath) eq 'ARRAY' ){
            $exist &&= &file_exist(filePath=>$_,alert=>$alert) for @$filePath;
        }
        elsif( ref($filePath) eq 'HASH' ){
            $exist &&= &file_exist(filePath=>$_,alert=>$alert) for keys %$filePath;
        }
        elsif( ref($filePath) eq 'SCALAR' ){
            $exist &&= &file_exist(filePath=>$$filePath,alert=>$alert);
        }
        else{ # wrong type
            $exist = 0;
        }
    }
    elsif( length($filePath) == 0
        || !-e abs_path( $filePath )
        || `file $filePath` =~ /broken symbolic link/
    ){
        $exist = 0;
    }

    cluck_and_exit "<ERROR>\tFile does not exist:\n".Dumper($filePath) if( $alert && $exist == 0 );
    return $exist;
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
    cluck_and_exit "<ERROR>\t$type command fails three times.\n"
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
        cluck_and_exit "<ERROR>\tThe version of java must be at least $minVer, but yours is $java_version\n";
    }
}

#--- try to `rm -rf` and then `mkdir -p` given folder path ---
sub reset_folder{
    shift @_ if(@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $folder = $parm{folder};
    &remove_folder(folder => $folder);
    `mkdir -p $folder`;
}

#--- try to `rm -rf` given folder path ---
sub remove_folder{
    shift @_ if(@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $folder = $parm{folder};
    `rm -rf $folder` if -e $folder;
}

1; ## tell the perl script the successful access of this module.
