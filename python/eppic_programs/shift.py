'''Functions inside:
shift - adjust an array that is periodic

'''

import numpy

def shift(array):
    '''Shift:

    Puprose: Adjust a 2D array that is periodic such that the sides are moved
    toward the center, i.e. the 4 corngers in the orignial array are now in the
    center of the array.

    Arguments:
    array - 2D array that will be shifted.

    Output:
    A 2D array of the same size, but with the values shifted.

    '''

    x = array.shape[0]
    y = array.shape[1]
    copy_array = numpy.copy(array)

    half_x = numpy.copy(copy_array[:x/2,:])
    otherhalf_x = numpy.copy(copy_array[x/2:,:])

    copy_array[:x/2,:] = otherhalf_x
    copy_array[x/2:,:] = half_x

    half_y = numpy.copy(copy_array[:,:y/2])
    otherhalf_y = numpy.copy(copy_array[:,y/2:])

    copy_array[:,:y/2] = otherhalf_y
    copy_array[:,y/2:] = half_y

    return copy_array
