'''Functions inside:
combine_domains - take same named files in different domains and combine them 
into a smaller set of larger files

'''
import readarray
import shutil
import os
import domainfiles
import numpy
import write_array

def combine_domains(filename, output_dir, order, combine = 2, paramfile = 'eppic.i'):
    '''Combine Domains
    Purpose: Combining same named files in different domains into a smaller set
    of larger files.

    Arguments:
    filename - string, name of file to be combined.
    output_dir - string, output path to place domains in.
    order - tuple or list, order in which dimensions will be read in (generally
    (2,1,0)).
    combine - integer, how many consecutive files to combine into one, try to
    make the combine divisible by the number of domains, automatically set to 2.
    paramfile - string, the file to read in to get parameters from,
    automatically set to eppic.i .

    Output: None, but result is a smaller set of larger files in a pathname, 
    output_dir.

    '''
    #get paramaters from eppic.i and extra paramaters
    params = eppic_input.read(filename = paramfile)
    ndim = []
    ndim.append(int(params['nx']))
    if 'ny' in params:
        ndim.append(int(params['ny']))
    if 'nz' in params:
        ndim.append(int(params['nz']))
    nout_avg = int(params['nout_avg'])
    nsubdomains = int(params['nsubdomains'])
    ndim_combine = []
    for i in ndim:
        ndim_combine.append(i/nout_avg)

    #make directories in which the combined files will be put
    domains = domainfiles.domainfiles(nsubdomains, filename)
    output_dir = output_dir + '/'
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
        for num in range(nsubdomains/2):
            os.mkdir(output_dir+'domain'+str(num).zfill(3)+'/')
    
    #combine files and put them in directories
    for i in (range(len(domains)/combine)):
        arrays = []
        for x in range(combine):
            arr = readarray.readarray(domains[x+(i*combine)], ndim_combine,
                                       order)
            arrays.append(arr)

        combine_arr = numpy.concatenate(tuple(arrays), axis = 0)
        print combine_arr.shape
        write_array.write_array(combine_arr, output_dir+domains[i], order)

    #copy and add in new variables for eppic.i
    shutil.copy('eppic.i', output_dir)
    eppic = open(output_dir+ 'eppic.i', 'ab')
    eppic.write('\nnx = ' + str(ndim_combine[0]*combine*nout_avg))
    eppic.write('\nnsubdomains = ' + str(nsubdomains/combine))
