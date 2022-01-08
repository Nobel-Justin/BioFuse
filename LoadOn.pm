package BioFuse::LoadOn;

use strict;
use warnings;
use FindBin qw/ $RealBin /;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              $V_Href
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::LoadOn';
#----- version --------
$VERSION = "0.50";
$DATE = '2018-11-01';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#------ variants -------
our $V_Href = {

    #++++++++++++++#
    # Main Version #
    #++++++++++++++#

    MainName => 'BioFuse.pl',
    Version => '0.13',
    Date => '2022-01-09',
    AUTHOR => 'Wenlong Jia',
    EMAIL => 'wenlongkxm@gmail.com',

    # functions
    RealBin => $RealBin,
    func => {},
    command => undef,
    argv_Aref => undef,
    run_mode => 0,

    ## manual
    HELP => 0,
    HELP_INFO => {},

    # for debug to keep some key intermediate folders/files
    in_debug => 0,

    #++++++++++#
    # settings #
    #++++++++++#

    # genome reference
    whole_genome => undef,
    cytoBand_file => undef,

    # gene annotation
    gtf => undef,
    psl => undef,
    Start_codon => ['ATG'],

    #++++++++++#
    # variants #
    #++++++++++#

    # intermediated containers


    #++++++++++#
    # software #
    #++++++++++#
};

#--------- functions in this pm --------#
my @functoion_list = qw/
                        load_variants_dict
                     /;

#--- load variant dict Href ---
sub load_variants_dict{
    return $V_Href;
}

#--- 
1; ## tell the perl script the successful access of this module.
