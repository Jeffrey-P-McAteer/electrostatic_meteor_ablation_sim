#!/usr/local/bin/perl

use bgl::test;
use strict;
use warnings;

if ($#ARGV<1) {
    print "Usage: bgltest.pl bgltest_file task\n";
    exit;
}
my $task = "nothing";
my $filename = $ARGV[0];
$task = $ARGV[1];
my $bgltest1;
if ($task eq "start") {
    $bgltest1 = bgltest_read($filename);
    bgltest_submit($bgltest1);
    bgltest_write($bgltest1,$filename);
} elsif ($task eq "kill") {
    $bgltest1 = bgltest_read($filename);
    bgltest_kill($bgltest1);
} elsif ($task eq "status") {
    $bgltest1 = bgltest_read($filename);
    my @status = bgltest_status($bgltest1);
    my $count=1;
    foreach my $stat (@status) {
	print "Job ($count) in ",$bgltest1->{jobs}->{names}->[$count-1],
	" has status '$stat'\n";
	$count++;
    }
} elsif ($task eq "collect") {
    $bgltest1 = bgltest_read($filename);
    my @results = bgltest_collect($bgltest1,["wall","sys","efield"]);
    my @param_keys = sort (keys %{$bgltest1->{params}});
    if ($#results>=0) {
	print sprintf(" %30.30s","name");
	print sprintf(" %10s","errTime");
	print sprintf(" %10s","wall");
	print sprintf(" %10s","sys");
	print sprintf(" %10s","efield");
	foreach my $key (@param_keys) {
	    next if ($key eq "options");
	    print sprintf(" %10s","$key");
	}
	print "\n";
	my $count=0;
	foreach my $result (@results) {
	print sprintf(" %30.30s",
		      $bgltest1->{jobs}->{names}->[$result->{jobnum}]);
	print sprintf(" %10.5g",$result->{errTime});
	print sprintf(" %10.5g",$result->{wall});
	print sprintf(" %10.5g",$result->{sys});
	print sprintf(" %10.5g",$result->{efield});
	foreach my $key (@param_keys) {
	    next if ($key eq "options");
	    print sprintf(" %10s",
			  $bgltest1->{params}->{$key}->[$result->{jobnum}]);
	}
	print "\n";
	}
    }
} else {
    print "You asked me to do $task, what did you mean by that ...?\n";
    print "Here are some useful tasks:\nsetup\nkill\nstatus\ncollect\n";
}


