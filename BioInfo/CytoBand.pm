package BioFuse::BioInfo::CytoBand;

use strict;
use warnings;
use BioFuse::Util::GZfile qw/ Try_GZ_Read /;
use BioFuse::Util::Log qw/ stout_and_sterr /;
use BioFuse::Util::Array qw/ Get_Two_Seg_Olen /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              load_cytoband
              get_cytoband
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::BioInfo::CytoBand';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-11-15';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @functoion_list = qw/
                        load_cytoband
                        get_cytoband
                     /;

#---- load cytoband info ---
sub load_cytoband{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $cytoBand_Href = $parm{cytoBand_Href};
    my $cytoBand_file = $parm{cytoBand_file};

    open (CBD, Try_GZ_Read($cytoBand_file)) || die"fail read cytoBand_file: $!\n";
    while(<CBD>){
        next if(/^#/);
        my ($refseg,$StP,$EdP,$band) = (split)[0,1,2,3];
        $cytoBand_Href->{$refseg}->{$.} = {band => $band, StP => $StP, EdP => $EdP};
    }
    close CBD;

    stout_and_sterr "[INFO]\tread cytoBand_file ok!\n"
                         ."\t$cytoBand_file\n";
}

#--- return the cytoband of given region on ref_seg ---
sub get_cytoband{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $cytoBand_Href = $parm{cytoBand_Href};
    my $refseg = $parm{refseg};
    my $regionStP = $parm{regionStP};
    my $regionEdP = $parm{regionEdP};

    my @band;
    for my $Idx (sort {$a<=>$b} keys %{$cytoBand_Href->{$refseg}}){
        my $bandHref = $cytoBand_Href->{$refseg}->{$Idx};
        push @band, $bandHref->{band} if Get_Two_Seg_Olen($bandHref->{StP}, $bandHref->{EdP}, $regionStP, $regionEdP);
    }

    return (scalar(@band) ? join(',',@band) : 'NA');
}

1; ## tell the perl script the successful access of this module.
