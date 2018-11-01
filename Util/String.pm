package BioFuse::Util::String;

use strict;
use warnings;
use BioFuse::Util::Log qw/ warn_and_exit /;
require Exporter;

#----- systemic variables -----
our (@ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
our ($VERSION, $DATE, $AUTHOR, $EMAIL, $MODULE_NAME);
@ISA = qw(Exporter);
@EXPORT = qw/
              getStrRepUnit
              getStrUnitRepeatTime
            /;
@EXPORT_OK = qw();
%EXPORT_TAGS = ( DEFAULT => [qw()],
                 OTHER   => [qw()]);

$MODULE_NAME = 'BioFuse::Util::String';
#----- version --------
$VERSION = "0.31";
$DATE = '2018-10-18';

#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

#--------- functions in this pm --------#
my @functoion_list = qw/
                        getStrRepUnit
                        getStrUnitRepeatTime
                     /;

#--- find given string's given-mode tandem-repeat unit and repeat-time ---
# search mode:
# w: whole-body exact-match, it's default
# s: start, from left
# e: end, from right
# g: global, no location required
sub getStrRepUnit{
    shift @_ if(@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $str = $parm{str};
    my $mode = $parm{mode} || 'w';

    $str = 'empty' unless( defined $str );
    my $sAnchor = ($mode =~ /^[ws]$/ ? '^' : '');
    my $eAnchor = ($mode =~ /^[we]$/ ? '$' : '');
    if( $str =~ /$sAnchor(.+)\1+$eAnchor/ ){
        my $repUnit = $1;
        my $repTime = length($&) / length($1);
        # 'g' mode needs to reverse and search again
        # as PERL regex is left-greedy
        if(    $mode eq 'g'
            && (reverse $str) =~ /$sAnchor(.+)\1+$eAnchor/
            && length($&) > length($repUnit) * $repTime
        ){
            $repUnit = reverse $1;
            $repTime = length($&) / length($1);
        }
        # find sub-repeats in repUnit
        my @unitAna = &getStrRepUnit( str => $repUnit, mode => 'w' );
        return ( $repTime * $unitAna[0], $unitAna[1] );
    }
    else{
        return ( 1, $str );
    }
}

#-------- calculate StrTest has how many StrUnit exact-match (allow partial) ---------#
## note that StrTest must wholely digested by StrUnit
sub getStrUnitRepeatTime{
    shift @_ if(@_ && $_[0] =~ /$MODULE_NAME/);
    my %parm = @_;
    my $StrTest = $parm{StrTest};
    my $StrUnit = $parm{StrUnit};
    my $TestEnd = $parm{TestEnd} || 'start'; # or 'end'
    my $CaseIns = $parm{CaseIns} || 0; # case in-sensitive
    my $MustInt = $parm{MustInt} || 0; # repeat time must be integer
    # my $CharAvo = $parm{CharAvo}; # Array-ref chars to avoid (not-implemented)

    # check
    if(    !defined $StrTest
        || !defined $StrUnit || length($StrUnit) == 0
        || ( $TestEnd ne 'start' && $TestEnd ne 'end' )
    ){
        warn_and_exit "<ERROR>\tWrong inputs in getStrUnitRepeatTime.\n"
                            ."\tStrTest: $StrTest\n"
                            ."\tStrUnit: $StrUnit\n"
                            ."\tTestEnd: $TestEnd\n"
                            ."\tCaseIns: $CaseIns\n";
    }

    # case in-sensitive
    if( $CaseIns ){
        $$_ = uc($$_) for ( \$StrTest, \$StrUnit );
    }

    # test
    my $repeat_time = 0;
    my $StrTestLen = length($StrTest);
    my $StrUnitLen = length($StrUnit);
    while( $StrTestLen >= $StrUnitLen ){
        my $subStrTest;
        if( $TestEnd eq 'start' ){
            $subStrTest = substr($StrTest, 0, $StrUnitLen);
            $StrTest = substr($StrTest, $StrUnitLen);
        }
        elsif( $TestEnd eq 'end' ){
            $subStrTest = substr($StrTest, -1 * $StrUnitLen);
            $StrTest = substr($StrTest, 0, -1 * $StrUnitLen);
        }
        # same as the StrUnit
        if( $subStrTest eq $StrUnit ){
            $repeat_time++;
            $StrTestLen -= $StrUnitLen;
        }
        else{
            return 0; # not match
        }
    }
    # whether StrTest has remains?
    if( $StrTestLen == 0 ){
        return $repeat_time;
    }
    elsif( $MustInt ){
        return 0; # do not allow paritial StrUnit
    }
    else{
        my $subStrUnit;
        if( $TestEnd eq 'start' ){
            $subStrUnit = substr($StrUnit, 0, $StrTestLen);
        }
        elsif( $TestEnd eq 'end' ){
            $subStrUnit = substr($StrUnit, -1 * $StrTestLen);
        }
        # same as the partial StrUnit
        if( $StrTest eq $subStrUnit ){
            $repeat_time += $StrTestLen / $StrUnitLen;
            return $repeat_time;
        }
        else{
            return 0;
        }
    }
}

1; ## tell the perl script the successful access of this module.
