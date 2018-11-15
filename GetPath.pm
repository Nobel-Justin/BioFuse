package BioFuse::GetPath;

use strict;
use warnings;
use File::Spec::Functions qw/ catfile /;
use BioFuse::Util::Log qw/ warn_and_exit /;
use BioFuse::LoadOn;

require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
my ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
	          GetPath
			/;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::GetPath';
#----- version --------
$VERSION = "0.01";
$DATE = '2018-11-15';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--- get path ---
sub GetPath{
	# options
	shift if (@_ && $_[0] =~ /$MODULE_NAME/);
	my %parm = @_;
	my $filekey = $parm{filekey};

    # BioFuse files
    if($filekey eq 'func_json'){
            return catfile($V_Href->{RealBin},'functions.json');
    }
    if($filekey eq 'BioFuseperlBin'){
            return catfile($V_Href->{RealBin},'BioFuse.pl');
    }
    # templete
	# if($filekey eq ''){
	# 	return catfile();
	# }

	# reach here, not found
	warn_and_exit "<ERROR>\tunrecognized filekey $filekey in $MODULE_NAME.\n";
}

#--- 
1; ## tell the perl script the successful access of this module.
