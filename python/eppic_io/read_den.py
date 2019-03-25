

from readdomain import readdomain

def read_den(eppicvars,*args,**keywords):
    '''Reads the domain decomposed densities from the domain*/*bin files

    Inputs:
    eppicvars  -- is a dictionary that contains the EPPIC input file in the
                  form of "variable" => "value" pairs. It\'s need to figure
                  out what the actual grid dimensions are.

    Outputs:
    denarrays  -- A list of density arrays.

    Optional Keywords
    dist       -- Specifies a single distribution to collect, otherwise all
                  distributions are gathered. 

    '''
    ndim = [1,1,1]
    nout_avg = 1
    ndist = 1
    if eppicvars.has_key('nout_avg'):
        nout_avg = int(eppicvars['nout_avg'])
    if eppicvars.has_key('nx'):
        ndim[0] = int(eppicvars['nx'])/nout_avg
    if eppicvars.has_key('ny'):
        ndim[1] = int(eppicvars['ny'])/nout_avg
    if eppicvars.has_key('nz'):
        ndim[2] = int(eppicvars['nz'])/nout_avg
    if eppicvars.has_key('ndist'):
        ndist = int(eppicvars['ndist'])
    denarrays = []
    if keywords.has_key('dist'):
        filename = 'den%d.bin'%keywords['dist']
        denarrays = readdomain(filename,ndim,int(eppicvars['nsubdomains']),
                              *args,**keywords)
    else:
        for dist in range(0,ndist):
            filename = 'den%d.bin'%dist
            denarrays.append(
                readdomain(filename,ndim,int(eppicvars['nsubdomains']),
                           *args,**keywords))

    return denarrays
