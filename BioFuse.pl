#!/usr/bin/perl -w
use strict;
use BioFuse::BioFuse;
use BioFuse::RunFunc;

my ($VERSION, $DATE, $AUTHOR, $EMAIL);

#----- version --------
$VERSION = "0.01";
$DATE = '2018-11-15';
#----- author -----
$AUTHOR = 'Wenlong Jia';
$EMAIL = 'wenlongkxm@gmail.com';

BioFuse::RunFunc->options_alert_and_run( argv_Aref => \@ARGV );