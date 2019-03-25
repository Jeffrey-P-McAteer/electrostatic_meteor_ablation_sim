import domainfiles
import numpy
import readarray_ftarray as read

## file found in ./src/output_FTarray.cc

## idl routines:
## ./idl/spectral/fill_k_array.pro -- takes structure and makes trad. array
## idl structure: fk_struct.ikx -> array of ikx values
##                fk_struct.iky -> array of iky values
##                fk_struct.fk  -> array of complexe values

def readdomain_ftarray(filename, ndim, nsubdomains, order = None,
                precision = 'single'):
    domains = domainfiles.domainfiles(nsubdomains, filename)
    domain_number = numpy.arange(1, nsubdomains)
    output_dictionary = {}
    output_array = numpy.zeros(ndim, 'complex')
    
    for i in domain_number:
        dictionary = read.readarray(domains[i], len(ndim.shape))
        output_dictionary.update(dictionary)
    
    for x in range(ndim[0]):
        for y in range(ndim[1]):
            if output_dictionary.has_key((x,y)):
                output_array[x,y]=output_dictionary[(x,y)]

    return output_array

### what gets written in output_FTarray.cc
    
#     unsigned short int ik[NDIM]={INDICIES(ikx_global,iky_global,ikz)};
#     fwrite(&ik, sizeof(ik[0]), NDIM, fname);
#     // Copy data to single precision
#     float data[2]={real(A_trans.cdata(INDICIES(iky,ikx,ikz))),
#                    imag(A_trans.cdata(INDICIES(iky,ikx,ikz)))};
#     fwrite(&data,sizeof(data[0]),2,fname);
#     ic++;


