'''Functions inside:
read - Read\'s a file, specifically eppic.i and makes a dictionary out of it.

'''

def read_input(filename="eppic.i"):
    """read_input

    Purpose: Parse EPPIC input from the file, filename, and returns a
    dictionary

    Arguments:
    filename - string, inherent value is 'eppic.i'

    Output: dictionary of strings containing the values contained in eppic.i 

    """
    
    eppic_var = {}
    # Test the filename value
    try:
        file = open(filename,"r")
    except IOError,e:
        print e
        raise Exception("The file, '" + filename + "' cannot open!")
    
    # parse the input
    filelines = file.read().splitlines()
    line = 0
    dist = -1
    while line < len(filelines):

        # drop quotes
        if filelines[line].lstrip().startswith(";"):
            del filelines[line]
        else:
            # drop blank lines
            if filelines[line].strip() == '':
                del filelines[line]
                continue
            else:
                # strip comments from back end, get var value pair
                try:
                    (var,value) = (
                        filelines[line].split(";",1).pop(0)
                        ).lstrip().split("=",1)
                    var=var.strip()
                    value=value.strip()
                except ValueError:
                    # drop cases where only var exists, ie no '='
                    del filelines[line]
                    continue
                
                # check to see if starting dist specific info
                if var == 'dist':
                    dist = int(value)
                    if not(eppic_var.has_key("dist"+str(dist))):
                        eppic_var["dist"+str(dist)] = {}
                    eppic_var["dist"+str(dist)][var] = dist
                    line+=1
                    continue
                else:
                    if var == 'fndist':
                        
                        dist = -1
                    else: 
                        if dist >= 0:
                            # append dist number if not there
                            if not(var.endswith(str(dist))):
                                var+=str(dist)
                            if not(eppic_var.has_key("dist"+str(dist))):
                                eppic_var["dist"+str(dist)] = {}
                            eppic_var["dist"+str(dist)][var]=value
                            line+=1    
                            continue
                    if eppic_var.has_key(var):
                        print '''
                        Warning: \'%s\' is changing value from
                        \'%s\' to \'%s\'''' % (var,
                                               eppic_var[var],value)
                    eppic_var[var] = value.strip()
                    line+=1


                

    return eppic_var

if __name__ == "__main__":
    eppic_var = read()
    print "\n".join(["%s='%s'" % (k, v) for k, v in eppic_var.items()])

    if eppic_var.has_key("nx"):
        if eppic_var.has_key("ny"):
            nx = int(eppic_var["nx"])
            ny = int(eppic_var["ny"])
            print ("Need to read in data that is %d by %d, or %d total points" 
                   % (nx,ny,nx*ny))
    if eppic_var.has_key("eps"):
        print "The value of eps^2 is %e" % (float(eppic_var["eps"])**2)
    eppic_var = read("wrongfile.in")
    print "\n".join(["%s='%s'" % (k, v) for k, v in eppic_var.items()])
    


