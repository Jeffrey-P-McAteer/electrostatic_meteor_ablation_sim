
import math
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import eppic_calc

def array_range(array,fit_type,take_abs = False):
    '''
    Find the range of relevant values in each dimension of array.

    array    - the numpy array

    fit_type - a vector the length of the array\'s dimensions,
               possible values:
                 0 - average other directions, fit to gaussian centered at
                     midpoint. Accept 2 standard deviations. 
                 1 - average other directions, fit to gaussian centered at
                     zero. Accept 2 standard deviations.
                     
    take_abs - take the absolute value of array before finding range. 
    '''

    lengths = array.shape
    ndim = len(lengths)
    if len(fit_type) != ndim:
        print '''
        Warning: wrong length for fit_type (%d should be %d), using zeros!
        '''%(len(fit_type),ndim)
        fit_type = np.zeros(ndim)
    for dim in range(ndim):
        print "for dim %d, array has length %d"%(dim,lengths[dim])
        if lengths[dim] < 2:
            continue
        sumArray = array.copy()
        if take_abs:
            sumArray = abs(sumArray)
        nsummed=0
        gridSummed = 1
        for sumDim in range(ndim):
            if sumDim == dim:
                continue
            sumArray = sumArray.sum(axis=sumDim-nsummed)
            gridSummed*=lengths[sumDim]
            nsummed+=1
        sumArray = eppic_calc.reform(sumArray)/gridSummed
        if fit_type[dim] == 0:
            fitFunc = lambda p,x : p[0]*p[0]*(np.exp(-1.*x*x/p[1]/p[1]))
            xarray = np.linspace(-lengths[dim]/2,lengths[dim]/2-1,lengths[dim])
            guessParams = [1.,1.]
        if fit_type[dim] == 1:
            fitFunc = lambda p,x : p[0]*p[0]*(np.exp(-1.*x*x/p[1]/p[1]))
            xarray = np.linspace(0,lengths[dim]-1,lengths[dim])
            guessParams = [math.sqrt(sumArray[0]),1.]

            
        errFunc = lambda p,x,y : fitFunc(p,x)-y
        fitParams, success = optimize.leastsq(errFunc,guessParams,
                                              args=(xarray,sumArray))
        print "Success? ",success
        print "Guess: ",guessParams
        print "Fit:   ",fitParams
        resultArray = fitFunc(fitParams[:],xarray)
#         subplot = plt.axes()
#         subplot.plot(xarray,resultArray,'r-')
#         subplot.plot(xarray,sumArray,'ro')
#         plt.show()
