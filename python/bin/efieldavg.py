#!/usr/bin/python
'''
efieldavg
A script to read in the efield array and plot its magnitude vs. time. By
default it plots sqrt(E^2), but it can be made to plot the energy from
conserved.out

@todo Make the plots prettier
@todo Make printing to file an automatic option
@todo Setup option for skipping data so plotting is faster
'''

def usage():
    print '''
USAGE: efieldavg [-ch] [input]
    
Options:
    -h | --help      : print this message
    -c | --conserved : use the Efield energy from conserved.out
    input            : an optional input file, which defaults to eppic.i in
                       the current directory
    '''
    return

import eppic_io
import eppic_calc
import sys
import matplotlib.pyplot as plt
import numpy as np
import getopt

# handle options

try:                                
    opts, args = getopt.getopt(sys.argv[1:], "hc", ["help", "conserved"]) 
except getopt.GetoptError:           
    usage()                          
    sys.exit(2)

useconserved=0
for opt, arg in opts:                
    if opt in ("-h", "--help"):
        usage()                     
        sys.exit()                  
    elif opt in ("-c", "--conserved"): 
        useconserved=1
    else:
        print "opt is ",opt," arg is ",arg
# get data

eppicvars = eppic_io.read_input()
if eppicvars == None:
    sys.exit()

dt = float(eppicvars['dt'])
nout = int(eppicvars['nout'])

time = []
efield = []
if useconserved==1:
    conserveddat = eppic_io.read_conserved()

## this is because read_conserved does not return the correct time yet
    time = np.array(conserveddat[0])*dt*nout
    efield = conserveddat[2]
    plt.ylabel('E-field Energy')
    plt.plot(time,efield)

else: # get data from phi and calcualted E^2
    phiarray = eppic_io.read_phi(eppicvars)
    efieldarray = eppic_calc.efield_from_phi(phiarray)
    efieldsqr = np.sqrt(efieldarray[0,...]**2+efieldarray[1,...]**2)
    efieldndim = efieldsqr.shape
    timedim = len(efieldndim)-1
    efieldsqr_max=[]
    efieldsqr_avg=[]
    for it in range(0,efieldndim[timedim]):
        efieldsqr_max.append(np.max(efieldsqr[...,it]))
        efieldsqr_avg.append(np.sum(efieldsqr[...,it]))
    efieldsqr_avg=(np.array(efieldsqr_avg)
                   /float(efieldsqr.size/efieldndim[timedim]))
    
    efield = efieldsqr_avg
    time = np.arange(0,efieldndim[timedim])*dt*nout
    plt.plot(time,efieldsqr_avg,label='E$^2$avg')
    plt.plot(time,efieldsqr_max,label='E$^2$max')
    plt.legend()

# make plot
plt.xlabel('time (s)')
plt.show()

