#!@SHELL@
# a small script to submit eppic jobs
# bsub -q your-queue twister.com your-eppic.i

#=====================  Configuration  =====================
# printing 
verbose=0

#===========================================================
# you should not need to edit below these lines
#===========================================================

# Variables set by configure
package_dir=@prefix@
package=`echo @PACKAGE@ | sed '@program_transform_name@'`

#================== Command line options ===================

while getopts 'vp:' option
  do
  case "$option" in 
      "p")      opt_p="$OPTARG"
	  ;;
      "v")      verbose=1
	  ;;
      ?)	logMesg "Error in reading arguents! Exiting ... " 
	  exit 1
	  ;;
  esac
done

# dump all the flags; shift off the beginning of the command line
shift `expr $OPTIND - 1`


#===========================================================
mk_domains() # a small function to create domain directories
#===========================================================
{
    i=0
    while [ $i -lt $1 ]; do
	name="domain00$i"
	if [ $i -le 10 ]; then; name="domain0$i"; fi
	if [ $i -le 100 ]; then; name="domain$i"; fi
	logMesg $name
	mkdir $name
	@ i +=1
    done
}




#===========================================================
logMesg() # a small function to record a message
#===========================================================
{
    if [ $verbose -ne 0 ]; then; echo "$1" >> $srcdir/g03.log; fi
}

# set the dir files will be copied from (to--afterwards)
srcdir=`pwd`

if [ "$HOSTNAME" = "" ]
then 
    HOSTNAME="twister"
fi

input="$1"
logMesg "#========Starting a eppic run======================#"
now=`date +"%r %a %d %h %Y"`
logMesg "Time: $now"
logMesg "Input: $input"

prefix=`basename $input '.i'`
if [ "$prefix" = "" ]
then
    echo "Wrong input!"
    exit 1
fi
logMesg "Prefix: $prefix"

scratchdir="/$HOSTNAME/scratch/$USER/$prefix" 
logMesg "Scratch Directory: $scratchdir"
if [ ! -d "$scratchdir" ] 
then 
	logMesg "-|Need to make scratch directory"
	mkdir -p "$scratchdir"
fi 

if [ -f "$scratchdir/$input" ]
then
    logMesg "<>Backing up old input "
    logMesg "  $scratchdir/$input "
    logMesg "  to $scratchdir/$input.bak"
    mv "$scratchdir/$input" $scratchdir/"$input.bak"
fi

cp $input $scratchdir/$input

nsubdomains=`sed -n -e s#nsubdomains *= *\(.*\) *$#\1#p`
logMesg "This run has $nsubdomains subdomians"
mk_domains $nsubdomains

logMesg "<>Changing directories to $scratchdir"
cd "$scratchdir"

logMesg "<>Executing eppic"
$package -procs $nprocs >! $package.log
logMesg "<>eppic done"
logMesg "<>Copying resulting files and .log file"
cp $package.log $srcdir
cp -r domain* $srcdir
logMesg "#==========Ending a eppic run======================#"