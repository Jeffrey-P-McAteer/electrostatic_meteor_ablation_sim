'''Functions inside:
readdomain - Read in multiple same named files in different domains and output 
an array.

'''

from readarray import readarray
from domainfiles import domainfiles
import numpy
import math
from timesteps import timesteps

def readdomain(filename, ndim, nsubdomains, order = None,
               precision = 'single',byte_order=None,
               nonp = [0,0],
               *extra_args,**extra_keywords):
    '''Read Domain
    Purpose:
    Read in multiple same named files in different domains and output an
     array. 

    Arguments:
    filename    - string, file to be read in
    
    ndim        - list or tuple, spatial dimensions of data in file (do not
                  add in the time dimension).
    
    nsubdomains - integer, number of domains from which files will be read in.
    
    order       - list or tuple, order in the dimensions of the array will be
                  read in, use the same order as what was used to write in.
                  Automatically sets itself to (1,0) if array is 2-d or
                  (2,1,0) if array is 3-d.
    
    precision   - string, set as either \'single\' or \'double\' depending on
                  what type the file is. Automatically set to \'single\'.

    byte_order  - a string telling how to handle the endian issue: \'S\', swap;
                  \'I\', ignore; \'L\', treat as little endian; \'B\', treat
                  as big endian.

    outdir      - directory that stores the domain*** directories, defaults to
                  \'.\'.

    nonp - The number of "ghost" cells to read on the first and last
                  domain. 
    Output:
    array       - An array of dimensions ndim with an extra time dimension and
                  the x dimension multiplied by nsubdomains.

    '''

    readarray_byte_order='I'

    domains = domainfiles(nsubdomains, filename,**extra_keywords)
    domain_num = numpy.arange(1,nsubdomains-1)


    ndimp = list(ndim)

    ndimp[0] = ndimp[0]+nonp[0]+nonp[1]
    
    smallestTime = timesteps(
        filename,
        nsubdomains=nsubdomains,
        gridsize=ndimp,
        nonp=nonp,#[0,0],
        **extra_keywords
        )
    
    startt=0
    endt=smallestTime
    skipt=1
    if extra_keywords.has_key('startt'):
        startt=extra_keywords['startt']
    if extra_keywords.has_key('endt'):
        endt=extra_keywords['endt']
        if endt < startt:
            endt=smallestTime
    if extra_keywords.has_key('skipt'):
         skipt=extra_keywords['skipt']


    outDtype = numpy.float32
    if precision == 'double':
        outDtype = numpy.float64

    totalTime = int(math.ceil((endt-startt)/skipt))
    
    try:
        ndimDomain = []
        for ix in range(len(ndim)):
            ndimDomain.append(ndim[ix])

        ndimDomain[0]*=nsubdomains
        ndimDomain[0]+=nonp[0]+nonp[1]

        output_array = numpy.zeros(ndimDomain+[totalTime],dtype=outDtype)
    except ValueError:
        raise ValueError


    if nsubdomains==1:
        new_array = readarray(domains[0], ndimp, order,
                              byte_order=readarray_byte_order,
                              *extra_args,**extra_keywords)
        output_array[0:ndimp[0],...,0:totalTime] = (
            new_array[...,0:totalTime])
    else:
        # first domain
        new_array = readarray(domains[0], ndimp, order,
                              byte_order=readarray_byte_order,
                              *extra_args,**extra_keywords)

        output_array[0:ndim[0]+nonp[0],...,0:totalTime] = (
            new_array[0:ndim[0]+nonp[0],...,0:totalTime])
        
        # middle domains
        i=0
        for i in domain_num:
            new_array =  readarray(domains[i], ndimp, order = order,
                                   precision = precision,
                                   byte_order=readarray_byte_order,
                                   *extra_args,**extra_keywords)

            xstart = (i)*ndim[0]+nonp[0]
            xend = (i+1)*ndim[0]+nonp[0]
            output_array[xstart:xend,...,0:totalTime] = (
                new_array[nonp[0]:nonp[0]+ndim[0],...,0:totalTime])
            
        #last domain
        i+=1
        new_array =  readarray(domains[i], ndimp, order = order,
                               precision = precision,
                               byte_order=readarray_byte_order,
                               *extra_args,**extra_keywords)
        
        
        xstart = (nsubdomains-1)*ndim[0]+nonp[0]
        xend = (nsubdomains)*ndim[0]+nonp[0]+nonp[1]
        output_array[xstart:xend,...,0:totalTime] = (
            new_array[nonp[0]:nonp[0]+nonp[1]+ndim[0],...,0:totalTime])
        


#                     # first domain
#     ndimThisdomain = list(ndim)
    
#     ndimThisdomain[0]+=nonp[0]
#     if nsubdomains==1:
#         ndimThisdomain[0]+=nonp[1]
#     new_array = readarray(domains[0], ndimThisdomain, order,
#                           byte_order=readarray_byte_order,
#                           *extra_args,**extra_keywords)


#     output_array[0:ndimThisdomain[0],...,0:totalTime] = (
#         new_array[...,0:totalTime])

#     # middle domains
#     for i in domain_num:

#         new_array =  readarray(domains[i], ndim, order = order,
#                                          precision = precision,
#                                          byte_order=readarray_byte_order,
#                                          *extra_args,**extra_keywords)

#         output_array[i*ndim[0]+nonp[0]:(i+1)*ndim[0]+nonp[0],...,0:totalTime] = (
#             new_array[...,0:totalTime])

#     #last domain
#     if nsubdomains>1:
#         ndimThisdomain = list(ndim)
#         ndimThisdomain[0]+=nonp[1]

#         new_array =  readarray(domains[i], ndimThisdomain, order = order,
#                                precision = precision,
#                                byte_order=readarray_byte_order,
#                                *extra_args,**extra_keywords)

#         xstart = (nsubdomains-1)*ndim[0]+nonp[0]
#         xend = (nsubdomains)*ndim[0]+nonp[0]+nonp[1]
#         output_array[xstart:xend,...,0:totalTime] = (
#             new_array[...,0:totalTime])


   ## check endian
    if byte_order == None:
        if numpy.isnan(numpy.max(output_array)):
            output_array = output_array.newbyteorder('S')
            if numpy.isnan(numpy.max(output_array)):
                # still nan so data is just bad, undo reorder
                output_array = output_array.newbyteorder()
    else:
        output_array = output_array.newbyteorder(byte_order)

    return output_array
        
