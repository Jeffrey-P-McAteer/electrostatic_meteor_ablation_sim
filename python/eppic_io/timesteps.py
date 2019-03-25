

from domainfiles import domainfiles
import os
from stat import ST_SIZE

def timesteps(filename,nsubdomains=1,gridsize=None,params=None,bytes=4,
              nonp=[0,0],**keywords):
    '''timesteps
    Counts the number of timesteps in a file.
    '''
    nonp=[0,0]
    if (gridsize == None) and (params == None):
        print "Can't process timesteps without either grid size or params!"
        return -1
    
    area = 1
    if gridsize != None:
        for length in gridsize:
            area*=length
    else:
        if params.has_key('nx'):
            area*=int(params['nx'])
        if params.has_key('ny'):
            area*=int(params['ny'])
        if params.has_key('nz'):
            area*=int(params['nz'])

    if nsubdomains==0:
        filelist = [filename]
    else:
        filelist = domainfiles(nsubdomains,filename,**keywords)


    size=0

    if len(filelist) == 1:
        thisarea = area/gridsize[0]*(gridsize[0]+nonp[0]+nonp[1])
        size = int(os.path.getsize(file)/thisarea/bytes)

    else:
        # then the first one is special
        file = filelist[0]
        thisarea = area/gridsize[0]*(gridsize[0]+nonp[0])
        filesize = (os.stat(file))[ST_SIZE]
        #size = int(os.path.getsize(file)/thisarea/bytes)
        size = int(filesize/thisarea/bytes)

        # middle are the same
        for ifile in range(1,len(filelist)-1):
            file = filelist[ifile]
            filesize = (os.stat(file))[ST_SIZE]
            thissize = int(filesize/area/bytes)
            if thissize < size:
                size = thissize

        # last is also special

        file = filelist[-1]
        thisarea = area/gridsize[0]*(gridsize[0]+nonp[1])
        filesize = (os.stat(file))[ST_SIZE]
        thissize = int(filesize/thisarea/bytes)
        if thissize < size:
            size = thissize

    return size
