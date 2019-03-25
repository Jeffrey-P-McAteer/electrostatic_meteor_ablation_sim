

def domainfiles(nsubdomains,filename,*extra_args,**extra_keywords):

    outdir='.'
    if extra_keywords.has_key('outdir'):
        outdir = extra_keywords['outdir']

    filelist = []
    for domain in range(0,nsubdomains):
        filelist.append(outdir+'/domain%03d/%s'%(domain,filename))
    return filelist
