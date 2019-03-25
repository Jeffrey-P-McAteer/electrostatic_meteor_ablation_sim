
def read_conserved(outdir="."):
    time = []
    particle_energy = []
    field_energy = []

    try:
        file = open(outdir+"/domain000/conserved.out","r")
    except IOError:
        print "The file, '" + outdir + "/domain000/conserved.out' cannot open!"
        return None
    filelines = file.read().splitlines()
    file.close()
    linenum=0
    for line in filelines:
        if (linenum == 0):
            linenum+=1
            continue
#        print "Fileline(%d): %s" % (linenum,line)
        vals = line.split()
        # in the future should be 
        #         time.append(float(vals[0])) # and change next two lines too
        time.append(float(linenum-1))
        particle_energy.append(float(vals[0]))
        field_energy.append(float(vals[1]))
        linenum+=1
    return(time,particle_energy,field_energy)


