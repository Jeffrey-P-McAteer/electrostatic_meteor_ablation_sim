#!/usr/local/bin/perl

# defensive coding
use strict;
use warnings;

# use module that has routines for testing
use bgl::test;

my $srcdir = $ARGV[0];

print "\n\nExample 1\n\n";
# -------------------------------------------------------------------------- #
# example 1
# -------------------------------------------------------------------------- #
# Uses template1.in
# Tests efield_petsc_algorithm=0 (iterative)
# Looks at the effect of increasing domains and copies, but
# by keeping the problem's size constant, eg number of particles, mesh size
# These could be long tests...

# some useful parameters
my $nxmax = 512;
my $npdrequired = 58982400;
my $npdmax =  14745600;
my $nout = 5;
my $nt = 20;
my $eq_algorithm=0; #iterative
my $ncopy=16;
my $options = "-log_summary";
my @sub_pctypes = ("jacobi","ILU","ICC");
my @pctypes = ("JACOBI","BJACOBI","ASM");
my %params=();

my $count=0;
print sprintf ("%10s %10s %10s %10s %10s %10s %10s %10s \n",
	       "count",
	      "ndomains",
	      "nproc",
	      "npd",
	      "npd/max",
	      "nx",
	      "nxtotal",
	      "ny");
for(my $itype = 0; $itype<scalar(@pctypes);$itype++) {
    for(my $nproc=$ncopy;$nproc<=128;$nproc*=2 ) {
	$params{"ndomains"}->[$count] = $nproc/$ncopy;
	$params{"nproc"}->[$count] = $nproc;
	$params{"npd"}->[$count] = $npdrequired/($nproc);
	next if ($params{"npd"}->[$count]>$npdmax);
	$params{"nx"}->[$count] = $nxmax/$nproc*$ncopy;
	$params{"ny"}->[$count] = $nxmax;
	$params{"nt"}->[$count] = $nt;
	$params{"nout"}->[$count] = $nout;
	$params{"efield_petsc_algorithm"}->[$count] = $eq_algorithm;
	
	$params{"options"}->[$count] = 
	    "$options -pc_type $pctypes[$itype]";
	print sprintf ("%10d %10.5g %10.5g %10.5g %10.5g ".
		       "%10.5g %10.5g %10.5g \n",
		       $count+1,
		       $params{"ndomains"}->[$count],
		       $params{"nproc"}->[$count],
		       $params{"npd"}->[$count],
		       $params{"npd"}->[$count]/14745600,
		       $params{"nx"}->[$count],
		       $params{"nx"}->[$count]*$params{"ndomains"}->[$count],
		       $params{"ny"}->[$count]);
	$count++;
    }
}


my $bgltest1 = bgltest_init("$srcdir/template1.in","bgltest1",\%params);
bgltest_write($bgltest1,"bgltest1.yml");

# print "\n\nExample 2\n\n";
# # -------------------------------------------------------------------------- #
# # example 2
# # -------------------------------------------------------------------------- #
# # Uses template1.in
# # Tests efield_petsc_algorithm=0 (iterative)
# # Looks at the effect of increasing domains and copies, but
# # by keeping the problem's size constant, eg number of particles, mesh size
# # These could be long tests...

# # some useful parameters
# $nxmax = 512;
# $npdmax = 14745600;
# $nout = 5;
# $nt = 20;
# $eq_algorithm=0; #iterative
# %params=();

# $count=0;
# print sprintf ("%10s %10s %10s %10s %10s %10s %10s %10s \n",
# 	       "count",
# 	      "ndomains",
# 	      "nproc",
# 	      "npd",
# 	      "npd/max",
# 	      "nx",
# 	      "nxtotal",
# 	      "ny");
# for(my $nproc=8;$nproc<=64;$nproc*=2 ) {
#     for( my $ncopy=4;$ncopy<=4;$ncopy*=2) {
# 	$params{"ndomains"}->[$count] = $nproc;
# 	$params{"nproc"}->[$count] = $nproc*$ncopy;
# 	$params{"npd"}->[$count] = $npdmax*4/($nproc*$ncopy);
# 	$params{"nx"}->[$count] = $nxmax/($nproc);
# 	$params{"ny"}->[$count] = $nxmax;
# 	$params{"nt"}->[$count] = $nt;
# 	$params{"nout"}->[$count] = $nout;
# 	$params{"efield_petsc_algorithm"}->[$count] = $eq_algorithm;
# 	print sprintf ("%10d %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g \n",
# 		       $count+1,
# 		       $params{"ndomains"}->[$count],
# 		       $params{"nproc"}->[$count],
# 		       $params{"npd"}->[$count],
# 		       $params{"npd"}->[$count]/14745600,
# 		       $params{"nx"}->[$count],
# 		       $params{"nx"}->[$count]*$params{"ndomains"}->[$count],
# 		       $params{"ny"}->[$count]);
# 	$count++;
#     }
# }

# #my $bgltest2 = bgltest_init("$srcdir/template1.in","bgltest2",\%params);
# #bgltest_write($bgltest2,"bgltest2.yml");

