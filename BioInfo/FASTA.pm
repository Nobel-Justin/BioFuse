package BioFuse::BioInfo::FASTA;

use strict;
use warnings;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::Log qw/ stout_and_sterr cluck_and_exit /;
use BioFuse::Util::Sys qw/ file_exist trible_run_for_success /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              read_fasta_file
              write_fasta_file
              BWA_index_fasta
              Faidx_Dict_fasta
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::FASTA';
#----- version --------
$VERSION = "0.32";
$DATE = '2019-05-02';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @functoion_list = qw/
                        read_fasta_file
                        write_fasta_file
                        BWA_index_fasta
                        Faidx_Dict_fasta
                     /;

#--- read fasta file and do something ---
sub read_fasta_file{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $FaFile = $parm{FaFile};
    my $needSeg_Href = $parm{needSeg_Href};
    my $subrtRef = $parm{subrtRef};
    my $subrtParmAref = $parm{subrtParmAref};

    # check fasta file existence
    cluck_and_exit "<ERROR>\tno fasta file supplied. (read_fasta_file)\n" unless defined $FaFile;
    file_exist(filePath=>$FaFile, alert=>1);

    # read fasta file and run sub-routine
    open (FASTA,Try_GZ_Read($FaFile))||die "fail read FaFile: $!\n";
    $/=">";<FASTA>;$/="\n"; # remove '>' prefix
    while(<FASTA>){
        chomp(my $segName = $_); # remove the last "\n"
        $/=">";
        chomp(my $segSeq = <FASTA>); # remove the last '>'
        $/="\n";
        # skip if not needed
        if (defined $needSeg_Href && !exists($needSeg_Href->{$segName})){
            next;
        }
        # remove all blanks
        $segSeq =~ s/\s+//g;
        # run sub-routine
        if (defined $subrtRef){
            my @parm = (segName=>$segName, segSeq_Sref=>\$segSeq);
            push @parm, @$subrtParmAref if (defined $subrtParmAref);
            &{$subrtRef}(@parm);
        }
    }
    close FASTA;
}

#--- write fa file, can extend some base at the end ---
## will not change the original sequence
sub write_fasta_file{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $seq_Sref = $parm{SeqSref};
    my $fa_file = $parm{FaFile};
    my $segname = $parm{SegName};
    my $linebase = $parm{LineBase} || 50;
    my $CL_Extend_Len = $parm{CircleExtLen} || 0;
    my $split_N_bool = $parm{split_N} || 0;

    # how to extend, or not
    my $new_seg = '';
    if($CL_Extend_Len > 0){
        $new_seg = $$seq_Sref . substr($$seq_Sref, 0, $CL_Extend_Len);
    }
    elsif($CL_Extend_Len < 0){
        $new_seg = substr($$seq_Sref, 0, length($$seq_Sref)+$CL_Extend_Len);
    }
    else{
        $new_seg = $$seq_Sref;
    }

    # create fasta file
    open (FA,Try_GZ_Write($fa_file)) || die "fail write $fa_file: $!\n";
    if( !$split_N_bool ){
        print FA ">$segname\n";
        print FA substr($new_seg, $_ * $linebase, $linebase)."\n" for ( 0 .. int( (length($new_seg)-1) / $linebase ) );
    }
    else{ # split N
        my @N_split_seg = split /N+/i, $new_seg;
        for (my $i = 0; $i <= $#N_split_seg; $i++){
            print FA ">$segname-part$i\n";
            print FA substr($N_split_seg[$i], $_ * $linebase, $linebase)."\n" for ( 0 .. int( (length($N_split_seg[$i])-1) / $linebase ) );
        }
    }
    close FA;
}

#--- construct BWA index ---
sub BWA_index_fasta{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $FaFile = $parm{FaFile};
    my $bwa = $parm{bwa};

    # check fasta file existence
    cluck_and_exit "<ERROR>\tno fasta file supplied. (read_fasta_file)\n" unless defined $FaFile;
    file_exist(filePath=>$FaFile, alert=>1);

    # at least BWA 0.7.13 automatically select index algriothm based on reference size.
    # check this link, https://www.biostars.org/p/53546/
    # generally, virus genome is smaller than 2GB, so '-a is' is preferable
    my $cmd = "$bwa index $FaFile 2>/dev/null";
    trible_run_for_success($cmd, 'IndexByBwa', {cmd_Nvb=>1, esdo_Nvb=>1});
}

#--- faidx and dict ---
sub Faidx_Dict_fasta{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $FaFile = $parm{FaFile};
    my $samtools = $parm{samtools};

    # check fasta file existence
    cluck_and_exit "<ERROR>\tno fasta file supplied. (read_fasta_file)\n" unless defined $FaFile;
    file_exist(filePath=>$FaFile, alert=>1);

    (my $dict_file = "$FaFile.dict") =~ s/\.fa\.dict/\.dict/;
    my $cmd = "($samtools faidx $FaFile 2>/dev/null) && ($samtools dict $FaFile 2>/dev/null 1>$dict_file)";
    trible_run_for_success($cmd, 'Faidx_Dict', {cmd_Nvb=>1, esdo_Nvb=>1});
}

1; ## tell the perl script the successful access of this module.
