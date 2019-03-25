import os
from stat import ST_SIZE
import scipy
import struct

def readarray(filename, dim = 2):
    stats = os.stat(filename)
    filesize = stats[ST_SIZE]
    f = open(filname, 'rb')
    output= {}
    if dim == 3:
        line_num = int(filesize/14)
        for i in range(line_num):
            x = f.read(14)
            line = struct.unpack('hhhff', x)
            output[(line[0],line[1],line[2])] = complex(line[3],line[4])
    if dim == 2:
        line_num = int(filesize/12)
        for i in range(line_num):
            x = f.read(12)
            line = struct.unpack('hhff', x)
            output[(line[0],line[1])] = complex(line[2],line[3])

    return output
