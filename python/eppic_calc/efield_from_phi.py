
import numpy
from reform import reform

def efield_from_phi(phi,dx=1,dy=1,dz=1):
    '''Caculates the electric field using phi, which it is assumed has a final
    dimension representing time.'''
    phishape = phi.shape
    efieldshape = []
    for dim in range(0,phi.ndim-1):
        if dim != 1:
            efieldshape.append(phishape[dim])
    efieldshape.append(len(efieldshape)-1)
    nt = phishape[phi.ndim-1]

    if len(efieldshape) == 3:
        efield = numpy.array(numpy.gradient(reform(phi[...,0]),dx,dy))
    elif len(efieldshape) == 4:
        efield = numpy.array(numpy.gradient(reform(phi[...,0]),dx,dy,dz))

    dim_with_time = list(efield.shape)
    dim_with_time.append(1)
    efield.shape = (dim_with_time)
    for it in range(1,nt):
        if len(efieldshape) == 3:
            thisefield = numpy.array(
                numpy.gradient(reform(phi[...,it]),dx,dy))
        elif len(efieldshape) == 4:
            thisefield = numpy.array(
                numpy.gradient(reform(phi[...,it]),dx,dy,dz))
        thisefield.shape = (dim_with_time)
        efield = numpy.concatenate((efield,thisefield),efield.ndim-1)

    return efield
