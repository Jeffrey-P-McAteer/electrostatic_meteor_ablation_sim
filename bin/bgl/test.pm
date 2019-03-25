


package bgl::test;

use strict;
use base 'Exporter';

our @EXPORT = qw(bgltest_init 
		 bgltest_read 
		 bgltest_write 
		 bgltest_submit
		 bgltest_kill
		 bgltest_status
		 bgltest_collect);
use YAML;
use Cwd;
use Date::Calc;

my $var = "0";
my $exepath = $ENV{'EPPIC_TEST_BIN_PATH'};

sub bgltest_init{
    my %new_bgltest=();
    my $templateName = shift;
    $new_bgltest{name} = shift;
    $new_bgltest{params} = shift;
    # test parameters, make sure arrays are the same length
    die "Need to have parameters to run tests!\n"
	if (scalar(keys %{$new_bgltest{params}})==0);
    
    my @param_keys = keys %{$new_bgltest{params}};
    my $nvals = $#{$new_bgltest{params}->{$param_keys[0]}};
    for (my $ikey=1; $ikey<scalar(@param_keys); $ikey++) {
	if ($nvals != $#{$new_bgltest{params}->{$param_keys[$ikey]}}) {
	    print "Number of tests for param ",$param_keys[$ikey],
	    " not equal to default param: ",$param_keys[0],"\n";
	    exit;
	}
    }

    # required params:
    die "'nproc' required parameter for tests!\n" 
	if (! defined $new_bgltest{params}->{nproc});
    if (! defined $new_bgltest{params}->{options}) {
	for (my $ival=0;$ival<$nvals;$ival++) {
	    $new_bgltest{params}->{options}->[$ival] = '';
	}
    }

    # put code here to build directories and files
    die "$templateName does not exist!\n"
	if ( ! -f "$templateName" );
    open TEMPLATE,"<","$templateName";
    my @templateLines = <TEMPLATE>;
    close TEMPLATE;
    my @params = keys %{$new_bgltest{params}};
    my $nruns = scalar(@{$new_bgltest{params}->{$params[0]}});

    # make the names after the parameter values
    my @names = ();
    for(my $irun=0;$irun<$nruns;$irun++){
	my $name = $new_bgltest{name};
	foreach my $param (@params){
	    next if ($param eq "options");
	    $name.="_".$new_bgltest{params}->{$param}->[$irun];
	}
	push @names,$name;
    }

    # make a jobs key, with only their names for now
    $new_bgltest{jobs}->{names}=\@names;

    my $cwd = getcwd;

    for(my $irun=0;$irun<$nruns;$irun++){
	my $name = $names[$irun];
	# make directory
	if (! -d $name) {
	    mkdir $name or die "Can't make directory $name";
	}
	
	# write input
	open TESTINPUT,">","$name/$name.i";
	my @values;
	foreach my $param (@params){
	    push @values,$new_bgltest{params}->{$param}->[$irun];
	}
	print TESTINPUT gen_input(\@templateLines,\@params,\@values);
	close TESTINPUT;
	
	# write sub script
	open TESTSUB,">","$name/$name.sub";

	print TESTSUB gen_subScript($name,"$cwd/$name",$exepath,
				    $new_bgltest{params}->{nproc}->[$irun],
				    $new_bgltest{params}->{options}->[$irun]);
	close TESTSUB;

	# delete old out err files if they exist
	unlink "$name/$name.out" if ( -f "$name/$name.out" );
	unlink "$name/$name.err" if ( -f "$name/$name.err" );
    }

    return \%new_bgltest;
}

sub bgltest_write{
    my $bgltest = shift;
    my $filename = shift;
    open WRITEFILE,">","$filename" or 
	die "Can't open $filename for bgltest_write";
    print WRITEFILE Dump($bgltest);
    close WRITEFILE;
}

sub bgltest_read{
    my $filename = shift;
    open READFILE,"<","$filename" or 
	die "Can't open $filename for bgltest_write";
    my $bgltest = Load(join '',<READFILE>);
    close CLOSEFILE;
    return $bgltest;
}

sub gen_input{
    my $templateLines = shift;
    my @templateCopy = @$templateLines;
    my $params = shift;
    my $values = shift;
    my $returnStrg = "";
    foreach my $line (@templateCopy) {
	for(my $iparam=0;$iparam<scalar(@$params);$iparam++){
	    my $paramCaps = "\U$$params[$iparam]";
	    #print "Substituting $paramCaps with $$values[$iparam]\n";
	    $line =~ s/\@$paramCaps\@/$$values[$iparam]/;
	}
	$returnStrg.=$line;
    }
    return $returnStrg;
}

sub gen_subScript{
    my $name = shift;
    my $jobpath = shift;
    my $exepath = shift;
    my $nproc = shift;
    my $options = shift;
    my $scriptStg = "# @ numprocs = $nproc\n".
	"# @ working_path = $jobpath\n".
	"# @ jobname = $name\n".
	"# @ exepath = $exepath\n".
	"# @ options = $options\n";
# big script, so using << 'EOF' notation to copy it in
    $scriptStg.=<< 'EOF'
#
### required, do not change
# @ job_type = parallel
#
### required, do not change
# @ executable = /bgl/BlueLight/ppcfloor/bglsys/bin/mpirun
#
### required, the arguments to mpirun, e.g.
# @ arguments = -verbose 1 -np $(numprocs) -cwd $(working_path) -exe $(exepath)/eppic -args "$(working_path)/$(jobname).i $(options)"
#
### required, you must specify a runtime limit for your job
### format: hours:minutes:seconds or minutes:seconds or seconds
# @ wall_clock_limit = 300:00
#
### a directory for mpirun io files (defaults to cwd at time of submission)
# @ initialdir = $(working_path)
#
### the next 3 are pathnames for mpirun's stdin, stdout, and stderr
### they can be fullpaths or relative to initialdir
### $(jobid) will be assigned by ll
# @ input = /dev/null
# @ output = $(jobname).out
# @ error = $(jobname).err
#
### the project that the job's time should be charged to 
### blank means to use your default project
### specify one of your other projects to override the default
### the command "groups" lists all of your projects
# @ group = eregion
#
### send email when job completes
# @ notification = complete
#
### required, do not change
# @ queue
EOF
;
	return $scriptStg;
}


sub bgltest_submit{
    my $bglTest=shift;
    my @jobnames = @{$bglTest->{jobs}->{names}};
    my @jobids = ();
    foreach my $name (@jobnames) {
	if (( -f "$name/$name.i" ) && (-f "$name/$name.sub")){
#	    llsubmit: Processed command file through Submit Filter: "/usr/local/etc/loadl/submit_filter".
#           llsubmit: The job "fe1.bgl.bu.edu.61332" has been submitted.
	    my $sysResult = `llsubmit $name/$name.sub 2>&1`;
	    my $jobid = 0;
	    print "Sysresult= $sysResult\n";
	    foreach my $line (split "\n",$sysResult) {
		if ($line =~ m/bu\.edu\.(.*)\"/) {
		    $jobid = $1;
		}
	    }
	    push @jobids,$jobid;
	}
    }
    $bglTest->{jobs}->{ids} = \@jobids;
    
}

my %statusToTxt=(
		 -2 => "unknown",
		 -1 => "failed",
		 0 => "never started",
		 1 => "succeeded",
		 2 => "running");


my %statusToNum=(
		 "unknown"   => -2,
		 "failed"    => -1,
		 "never started" =>  0,
		 "succeeded" =>  1,
		 "running"   =>  2);


sub bgltest_single_status{
    my $jobname=shift;

    
    # are the out,err files created
    my $outfile="";
    my $errfile="";
    if ( -d $jobname ){
	$outfile = "$jobname/$jobname.out"
	    if ( -f "$jobname/$jobname.out" );
	$errfile = "$jobname/$jobname.err"
	    if ( -f "$jobname/$jobname.err" );
    }
    if (($outfile eq "") && ($errfile eq "")) {    
	return $statusToNum{"never started"};
    }
    # collect times from err file, timers.dat
    open ERR,"<","$errfile" or die "Can't open $errfile\n";
    my @errlines = <ERR>;
    close ERR;

    open OUT,"<","$outfile" or die "Can't open $outfile\n";
    my @outlines = <OUT>;
    close OUT;

    # is the err file finished
    my $has_finished=0;
    foreach my $line (@errlines) {
	if ($line =~ m/Exit status/){
	    $has_finished++;
	}
    }
    if ($has_finished==0){
	return $statusToNum{"running"};
    }

    # did it finish with a nice message in the out file
    my $finished_happy=0;
    foreach my $line (@outlines) {
	$finished_happy++
	    if ($line =~ m/EPPIC ending normally at:/);
    }
    if ($finished_happy) {
	return $statusToNum{"succeeded"};
    } else {
	return $statusToNum{"failed"};
    }
    return  $statusToNum{"unknown"};
}

sub bgltest_status{
    my $bgltest=shift;
    my @status;
    foreach my $jobname (@{$bgltest->{jobs}->{names}}){
	push @status,$statusToTxt{bgltest_single_status($jobname)};
    }
    return @status;
}


sub bgltest_collect{
    # get and test parameters
    my $bgltest=shift;
    my $timersParams=shift;
    
    my @resultList;
    my $count=0;
    foreach my $jobname (@{$bgltest->{jobs}->{names}}){
	push @resultList, bgltest_collect_onejob($jobname,$count,$timersParams)
	    if (bgltest_single_status($jobname)==
		$statusToNum{"succeeded"});
	$count++;
    }
    return @resultList;
}

sub bgltest_collect_onejob{
    my $jobname=shift;
    my $jobnum=shift;
    my $timersParams=shift;

    # make defaults incase job is not done
    my %results=(errTime => 0);
    $results{jobnum} = $jobnum;

    foreach my $param (@$timersParams) {
	$results{$param} = 0;
    }

    # another safty check
    if (bgltest_single_status($jobname)!=
	$statusToNum{"succeeded"}){
	print "Returning with no resutls because of status\n";
	 return \%results;
     }

     # search for out, err and timers.dat files
     my $outfile="$jobname/$jobname.out";
     my $errfile="$jobname/$jobname.err";
     my $timersfile="";
     if ( -d "$jobname/domain000" ) {
	 $timersfile = "$jobname/domain000/timers.dat"
	     if ( -f "$jobname/domain000/timers.dat" );
     } else {
	 $timersfile = "$jobname/timers.dat"
		 if ( -f "$jobname/timers.dat" );
     }
     if ( $timersfile eq "") {
	 print "There's something wrong with timers.dat file ",
	 "for job $jobname\n";
	 return \%results;
     }

     # collect times from err file, timers.dat
     open ERR,"<","$errfile" or return %results;
     my @errlines = <ERR>;
     close ERR;

     my $start_unformated;
     my $end_unformated;
     foreach my $line (@errlines) {
	 if ($line =~ m/<(.*)>.*IO - Threads initialized/) {
	     $start_unformated = $1;
	 }
	 if ($line =~ m/<(.*)>.*Job successfully terminated/) {
	     $end_unformated = $1;
	 }
     }
    
    $results{errTime} = (time_from_unformated($end_unformated)- 
			 time_from_unformated($start_unformated));

    
# TIMERS:
    open TIMERS,"<","$timersfile";
    
    my @timerslines = <TIMERS>;
    close TIMERS;

# Drop headings line
    shift @timerslines;
# Drop 0th time step
    shift @timerslines;
# manually enter headings
    my %keysnames = (
		     "iter" => 0,
		     "wall" => 1,
		     "sys"  => 2,
		     "vadv" => 3,
		     "xadv" => 4,
		     "charge" => 5,
		     "collect" => 6,
		     "efield" => 7,
		     "output" => 8,
		     "fluid"  => 9);
    # values will hold averages of the above headings
    my @values=();

    # sum values
    my $count=0;
    foreach my $line (@timerslines) 
    {
	chomp(my @localvalues = split " ",$line);
	$count++;
	for (my $ival=0;$ival<1+scalar(keys %keysnames);$ival++)
	{
	    $values[$ival]+=$localvalues[$ival];
	}
    }
    
    return \%results if ($count == 0) ;
    foreach my $param (@$timersParams)
    {
	$results{$param} = $values[$keysnames{$param}]/$count
    }
    
    # return hash with start, end, average of each timers.dat col.
    return \%results;
}

sub bgltest_kill{
    my $bglTest = shift;
    foreach my $jobid (@{$bglTest->{jobs}->{ids}}) {
	my $sysResult = `llcancel $jobid`;
	foreach my $line (split "\n",$sysResult) {
	    print "Outline: $line \n";
	}
    }
}

sub time_from_unformated{
    my $unformated_string=$_[0];
    my $year = Date::Calc::This_Year();
    my @datevals = split ' ',$unformated_string;
    my $month = Date::Calc::Decode_Month($datevals[0]);
    my $day = $datevals[1];
    (my $hour, my $minute, my $seconds) = split ":",$datevals[2];
    (my $second, my $nanosecond) = split '\.',$seconds;
    my $time = Date::Calc::Date_to_Time($year,
					$month,
					$day,
					$hour,
					$minute,
					$second);
    return $time;
}



1;
