package BioFuse::Util::Web;

use strict;
use warnings;
use LWP::UserAgent;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              urlToHtmlText
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Util::Web';
#----- version --------
$VERSION = "0.02";
$DATE = '2021-01-10';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';


#--------- functions in this pm --------#
my @function_list = qw/
                        urlToHtmlText
                     /;

#--- query webpage text by url ---
sub urlToHtmlText{
    # options
    shift if (@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $url = $parm{url};
    my $timeout = $parm{timeout} || 30;

    my $uagent = LWP::UserAgent->new;
    $uagent->timeout($timeout);
    $uagent->env_proxy;
    my $response = $uagent->get($url);
    if ($response->is_success) {
        return "HTTP/1.0 200 OK\n" . $response->decoded_content;
    }
    else{
        # try curl
        my $html_text = `curl '$url' 2>/dev/null`;
        if( length($html_text) != 0 ){
            return "HTTP/1.0 200 OK\n" . $html_text;
        }
        else{
            return 'fail';
        }
    }
}

1; ## tell the perl script the successful access of this module.
