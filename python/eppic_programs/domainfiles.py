'''Functions inside:
domainfiles - Makes a list containing domain000/(specified filename) to a specified number.

'''

def domainfiles(nsubdomains,filename):
    '''Domainfiles
    Purpose: Make a list containing domain000/filename to a specified number, 
    e.g. [domain000/filename, domain001/filename, ..., domain035/filename]
   
    Arguments:
    nsubdomains - integer, The number of domains inside list.
    filename - string, The specified filename to be appended to each domain
   
    Output: List containing domain000/filename to a specified number, 
    e.g. [domain000/filename, domain001/filename, ..., domain035/filename]

    '''

    domains = []
    for idomain in range(0,nsubdomains):
        # str() converts numbers to strings
        # zfill() creates a string of a certain width, by padding 0's to the 
        # front
        domains.append("domain"+str(idomain).zfill(3)+"/"+filename)
                
    return domains


if __name__ == "__main__":
    domainfiles(16,'den0.bin')
