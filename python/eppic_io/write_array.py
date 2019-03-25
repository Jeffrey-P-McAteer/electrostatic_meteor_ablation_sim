"""Functions:
write_array - takes an array and converts it into a file

"""
import numpy

def write_array(array, filename, order = None):
    """Write Array
    Purpose: Take in an array, and create a raw binary file containing the data
    Arguments:
    array - arrary, array to be used
    filename - string, file and perhaps pathname in which the array will 
    converted to.
    order - tuple or list, order in which the dimensions of the array will be 
    read in, automatically set to (1,0) if 2-d, and (2,1,0) if 3-d. 
    Output: None, but result is a binary file containing array to be read in 
    with order, order, and in pathway, filename.
    
    """
    #If no order is inputted, order will be (1,0) if array is 2-d or (2,1,0) if
    #array is 3-d
    if order == None:
        array_dim = array.shape
        if len(array_dim) == 3:
            order = (1,0)

        if len(array_dim) == 4:
            order = (2,1,0)

    #Re-adjusting order so that it will read in properly
    new_order = list(order)
    new_order.append(len(order))
    new_order.reverse()
    
    #Transpose array's dimensions according to new order.
    new_array = numpy.transpose(array, new_order)

    #Write array to file.
    new_array.tofile(filename)
