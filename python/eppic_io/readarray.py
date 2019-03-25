"""Functions inside: 
readarray - Turns a data file into an array
"""
import math
import numpy as np
import os
from stat import ST_SIZE

def readarray(filename, ndim, order = None, precision = 'single',skipt=1,
              startt=0,endt=-1,byte_order=None,*extra_args,**extra_keywords): 
   """Read Array
   Purpose: Take a binary file of data and output an array

   Arguments:
   filename  - string, file and perhaps pathname to be read in.
   ndim      - tuple or list, dimensions of array (consider only the spatial 
               dimensions, the time dimension is determined inside the
               function).
   order     - tuple or list, order in which the dimensions will be read in.
               Use the same order as was used to write in, regardless of
               whether it was through idl or python. Automatically set to (1,0)
               if 2-d or (2,1,0) if 3-d
   precision - string, set as either 'single' or 'double' depending on what
               type of precision the file is. Automatically set to single.
   skipt     - the cadence of reading (default eq 1)
   startt    - number of time steps to skip at the beginning (default 0)
   endt      - number of time steps to read (default -1); endt<=0, read all
               possible; endt>0, read endt steps from startt
   

   Output:
               An array of size ndim with an additional time dimension.

   """


   #If no order is inputted, order will be (1,0) if array is 2-d or (2,1,0) if 
   #array is 3-d
   if order == None:
      order = tuple([i for i in range(len(ndim)-1,-1,-1)])

   #Dealing with precision issues between single and double
   if precision == 'single':
      bitnum = 4
      dtype = np.float32
   
   if precision == 'double':
      bitnum = 8
      dtype = np.float64

   #targetNdim are the dimensions for each timestep to be read in, real_ndim is 
   #output_dimensions with an extra dimension of one added to represent time.
   targetNdim = list(ndim)
   real_ndim = targetNdim[:]
   real_ndim.append(1)
   order = list(order)
   rOrder = list(order)
   rOrder.reverse()

   recordNdim = [targetNdim[i] for i in order]
   
   #The amount of units being read in per unit time.
   num_per_t = 1
   for x in real_ndim:
      num_per_t *= x


   #The axis along which the timesteps will be concatenated(the last axis:time)
   axis = len(real_ndim)-1
   
   #Determining time dimension
   stats = os.stat(filename)
   filesize = stats[ST_SIZE]
   totaltime = int(filesize/(bitnum*num_per_t))

   if endt > totaltime:
      raise ValueError("\nendt is too large (%d),"%endt+
                       "\nmust be no larger than %d"%totaltime)
      
   timedim = totaltime/skipt
   timedim -= int(math.ceil(startt/skipt))

   if (endt > 0):
      timedim -= int(math.ceil((totaltime-endt)/skipt))

   try:
      output_array = np.zeros((ndim)+[timedim],dtype=dtype)
   except ValueError:
      raise ValueError

      
   #Initial portion of unit time to be read in
   fileobj = open(filename, 'rb')
   fileobj.seek(startt*num_per_t*bitnum,0)


   readArray = np.fromfile(file = fileobj,
                             dtype = dtype,
                             count = num_per_t)
   output_array[...,0] = np.transpose(
      np.reshape(readArray,targetNdim)
      ,rOrder)

   #Rest of file read in unit by unit, until it hits the end of file.
   for i in range(1,timedim):
      # skip fixed ammount
      fileobj.seek((skipt-1)*num_per_t*bitnum,1)

      inarray = np.fromfile(file = fileobj,
                            dtype = dtype,
                            count = num_per_t)
      try:
         reshapeArray = np.reshape(inarray, targetNdim)
      except ValueError as err:
         msg = (str(err)
                + "\nin shape: "
                + str(inarray.shape)
                + "\ttarget shape: "
                + str(targetNdim)
                + "\tnum_per_t: "
                + str(num_per_t)
                + "\tat time: "
                + str(i))
         raise ValueError(msg)
      output_array[...,i] = (np.transpose(reshapeArray, rOrder))


   ## check endian
   if byte_order == None:
      if np.isnan(np.max(output_array)):
         output_array = output_array.newbyteorder('S')
         if np.isnan(np.max(output_array)):
            # still nan so data is just bad, undo reorder
            output_array = output_array.newbyteorder()


   return output_array

