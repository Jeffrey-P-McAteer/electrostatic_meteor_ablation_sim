#!/bin/sh

eppic="../../src/eppic"
first_input="$srcdir/restart_test1_parta.i"
second_input="$srcdir/restart_test1_partb.i"


$eppic $first_input || exit 1
$eppic $second_input || exit 1
exit 0



