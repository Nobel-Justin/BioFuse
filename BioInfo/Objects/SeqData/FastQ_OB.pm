package BioFuse::BioInfo::Objects::SeqData::FastQ_OB;

use strict;
use warnings;
use Cwd qw/ abs_path /;
use Data::Dumper;
use BioFuse::Util::Log qw/ cluck_and_exit /;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::Objects::SeqData::FastQ_OB';
#----- version --------
$VERSION = "0.01";
$DATE = '2019-05-25';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        new
                        filepath
                        tag
                        rgOB_Hf
                        start_read
                        start_write
                        stop_write
                        write
                        add_rgOB
                     /;

#--- structure of object
# fq -> filepath = $filepath
# fq -> tag = $tag
# fq -> tissue = $tissue
# fq -> write_fh = $write_fh
# fq -> rgOB = { rgID->$rgOB }, check BioFuse::BioInfo::Objects::SeqData::ReadsGroup_OB

#--- construction of object ---
sub new{
    my $type = shift;
    my %parm = @_;

    my $fq = {};
    $fq->{filepath} = $parm{filepath} || undef;
    $fq->{tag} = $parm{tag} || undef;
    $fq->{tissue} = $parm{tissue} || undef;
    $fq->{rgOB} = {};

    # $fq->{filepath} = abs_path $fq->{filepath} if defined $fq->{filepath};

    bless($fq);
    return $fq;
}

#--- return file path ---
sub filepath{
    my $fq = shift;
    return $fq->{filepath};
}

#--- return fq tag ---
sub tag{
    my $fq = shift;
    return $fq->{tag};
}

#--- return Hash-ref fq's rgOB ---
sub rgOB_Hf{
    my $fq = shift;
    return $fq->{rgOB};
}

#--- open file-handle to start reading ---
sub start_read{
    my $fq = shift;
    my %parm = @_;

    # check existence
    unless(file_exist(filePath => $fq->filepath)){
        cluck_and_exit "<ERROR>\tCannot find fq\n".Dumper($fq);
    }
    # open fh
    open (my $readFH, Try_GZ_Read($fq->filepath)) || die "fail reading: $!\n".Dumper($fq);

    return $readFH;
}

#--- open file-handle to start writing ---
sub start_write{
    my $fq = shift;
    my %parm = @_;

    # check
    if(exists $fq->{write_fh}){
        cluck_and_exit "<ERROR>\tthe write file-handle of fq_OB already exists.\n".Dumper($fq);
    }
    # open fh
    open ($fq->{write_fh}, Try_GZ_Write($fq->filepath)) || die "fail writing: $!\n".Dumper($fq);
}

#--- close file-handle to stop writing ---
sub stop_write{
    my $fq = shift;
    close  $fq->{write_fh} if exists $fq->{write_fh};
    delete $fq->{write_fh};
}

#--- write contect to fq ---
sub write{
    my $fq = shift;
    my %parm = @_;
    my $content = $parm{content};
    print {$fq->{write_fh}} $content;
}

#--- add reads group object(s) to this fq ---
sub add_rgOB{
    my $fq = shift;
    my %parm = @_;
    $fq->{rgOB}->{$parm{rgOB}->RG_ID} = $parm{rgOB};
}

1; ## tell the perl script the successful access of this module.
