

from readdomain import readdomain

def read_phi(eppicvars,*extra_args,**extra_keywords):
    ndim = [1,1,1]
    nout_avg = 1
    if eppicvars.has_key('nout_avg'):
        nout_avg = int(eppicvars['nout_avg'])
    if eppicvars.has_key('nx'):
        ndim[0] = int(eppicvars['nx'])/nout_avg
    if eppicvars.has_key('ny'):
        ndim[1] = int(eppicvars['ny'])/nout_avg
    if eppicvars.has_key('nz'):
        ndim[2] = int(eppicvars['nz'])/nout_avg
    nonperiodic=[0,0]


    if eppicvars.has_key('boundary_type'):
        if int(eppicvars['boundary_type']):
            nonperiodic=[1,2]

    phiarray = readdomain('phi.bin',ndim,int(eppicvars['nsubdomains']),
                          nonp=nonperiodic,
                          *extra_args,**extra_keywords
                          )

    return phiarray
