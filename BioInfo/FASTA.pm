package BioFuse::BioInfo::FASTA;

use strict;
use warnings;
use BioFuse::Util::GZfile qw/ Try_GZ_Read Try_GZ_Write /;
use BioFuse::Util::Log qw/ stout_and_sterr warn_and_exit /;
use BioFuse::Util::Sys qw/ file_exist /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              read_fasta_file
              write_fasta_file
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::FASTA';
#----- version --------
$VERSION = "0.31";
$DATE = '2018-10-30';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @functoion_list = qw/
                        read_fasta_file
                        write_fasta_file
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
    if( !defined $FaFile){
        warn_and_exit "<ERROR>\tno fasta file supplied. (read_fasta_file)\n";
    }
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

    my $Option_Href = (@_ && $_[0] =~ /$MODULE_NAME/) ? $_[1] : $_[0];

    my $seq_Sref = $Option_Href->{SeqSref};
    my $fa_file = $Option_Href->{FaFile};
    my $segname = $Option_Href->{SegName};
    my $linebase = $Option_Href->{LineBase} || 50;
    my $CL_Extend_Len = $Option_Href->{CircleExtLen} || 0;
    my $split_N_bool = $Option_Href->{split_N} || 0;

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

1; ## tell the perl script the successful access of this module.
