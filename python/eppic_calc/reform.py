
import numpy

def reform(array):
    '''Like IDLs reform routine, it drops dimensions that have length 1'''

    newdim = []
    for dim in array.shape:
        if dim != 1:
            newdim.append(dim)

    newarray = array.copy()
    newarray.shape = (newdim)
    return newarray
