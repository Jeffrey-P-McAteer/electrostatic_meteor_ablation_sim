# get version from subversion replacement variable
__version__ = "$Rev$"

#from matplotlib import rc
#rc('text',usetex=True)
import numpy as np
import matplotlib.pyplot as plt

def image_plot(in_array, subplot= None,**params):
    '''Plots the array using a color map. This routine mimicks that
    originally created for IDL.

    Arguments:
    in array      -- 2D array or tuple/list of arrays to be plotted

    subplot       -- If defined, it is used to make the plot, otherwise a
                     new figure and subplot are created.

    Keywords
    x             -- a two element list describing the extent of the x
                     direction. x[0] is the smallest x value and x[1] is the
                     largest. 

    y             -- same as x, but for the y direction. 
    
    aspect        -- integer or string, ratio of y-scale to x-scale, useful 
                     for true to units scales; for example if the the y range 
                     is 2x as large then the aspect should be 2 so that the 
                     grid spacing in the image remains true to the unit size,
                     or use 'equal' to automatically make units equal. 
                     
    legend        -- string, say 'yes', or ignore the option, creates a color 
                     bar on the left hand side

    title         -- string, Plot title

    [xyz]title    -- string, [XYZ]-axis title

    cbtitle       -- string, legend (colorbar, ie cb) title

    nlevels       -- integer, number of contour lines

    zlog          -- plot the z range in log base 10

    interp        -- intperolation is used.

    noshow        -- don\'t show the plot,but return it. 

    xrange        -- range of x values to plot; plots closest possible value
                     if the range does not correspond directly to an array
                     cell.

    yrange        -- range of y values to plot; some constraints as xrange. 

    zrange        -- range of z values to plot

    colormap      -- Force a specific color map. 

    Options to be added?

    contour       -- string, 'yes' if you want a contour plot, ignore otherwise

    decades      -- an alternative means to set the z range; number of
                     decades smaller than the max value to plot.
                     Especially useful for zlog plots.

    palette       -- set the color palette, acceptable values are: gray, color
                     ... more to come

    filename      -- output file name; if it does not exist, the plot is to
                     the screen


    format        -- output format. This is ignored if filename paramter is
                     not set. Possible values: ps (postscript),
                     png (portable network graphics), and svg (scalable
                     vector graphics). The default is \'ps\'; there are
                     actually many other formats tha gnuplot can handle, but
                     these are the only formats setup in this routine.

    target        -- output target; this has different meanings depending on
                     the format.

                     For postscript, there are three choices:
                     fullpage, which produces a landscape plot that takes the
                     full page; pagewidth, which produces a 6 in wide figure;
                     and colwidth, which produces a 3 in wide figure for two
                     colomn publications. If fullpage is used, it will be
                     plotted in postscript format, otherwise, it will be
                     plotted in encapsulated format.

                     For png and svg, the output is intended to be for the web,
                     so small(640x480), medium(800x600) and large(1024x768)
                     options are available, though if one of the postscript
                     options are given (fullpage, pagewidth, or colwidth), then
                     a 300dpi figure of that size is produced. 
    
    
    Obsolete Options:
    window_scale  -- scales image to the size of window by default, otherwise
                     scales window to image size.

    nlabels       -- number of labels on the legend/color bar.

    zoom          -- zoom to the "most interesting" part of the image.


    '''


    if subplot == None:
        subplot = plt.axes()

    if params == None:
        params = {}
        
    array_size = in_array.shape
    array_copy = in_array[:,:]
        
    if params.has_key('zlog'):
        array_copy = np.log10(in_array)

    plot_array = np.transpose(array_copy)
    
    #----------------------------------------------
    # some notes about variables, same apply for y
    #----------------------------------------------
    # xidx    -> the range of the array that is plotted.
    #            if xsize is 10, and all values are plotted, xidx is [0,11]
    # x       -> the x values corresponding to the 1st and last array elements
    # xrange  -> optional values in x units used to plot portion of the xrange
    # dx      -> the conversion scale factor from xidx to x
    # xextent -> the value of x limits, just smaller/larger than x or xrange by
    #            half a data point (dx)
        
    xsize= array_size[0]
    ysize = array_size[1]


    if params.has_key('x'):
        x = [float(params['x'][0]),float(params['x'][1])]
    else:
        x = [0.,float(xsize)-1]

    if params.has_key('y'):
        y = [float(params['y'][0]),float(params['y'][1])]
    else:
        y = [0.,float(ysize)-1]

    xRange = list(x)
    yRange = list(y)
    if params.has_key('xrange'):
        xRange[0] = float(params['xrange'][0])
        xRange[1] = float(params['xrange'][1])
    dx = (x[1]-x[0])/float(xsize-1)

    if params.has_key('yrange'):
        for idx,rval in enumerate(params['yrange']):
            yRange[idx] = float(rval)
    dy = (y[1]-y[0])/float(ysize-1)

    xidx = [int(max(0,xRange[0]-x[0])/dx),
            int(min(x[1]-x[0],xRange[1]-x[0])/dx)+1]
    if xidx[1]<xidx[0]:
        raise "error in xrange! "
    xextent = [x[0]+(xidx[0]-.5)*dx,x[0]+(xidx[1]-.5)*dx]

    yidx = [int(max(0,yRange[0]-y[0])/dy),
            int(min(y[1]-y[0],yRange[1]-y[0])/dy)+1]
    if yidx[1]<yidx[0]:
        raise "error in yrange! "
    yextent = [y[0]+(yidx[0]-1.0)*dy,y[0]+(yidx[1]-1.0)*dy]
    
    if params.has_key('nlevels'):
        nlevels = params['nlevels']

    if params.has_key('contour'):
        raise "'contour' option not coded yet"

    map = plt.cm.jet
    if params.has_key('colormap'):
        map = params['colormap']
    
    ## remember plot_array is transposed (for some god awful reason)
    im = subplot.imshow(plot_array[yidx[0]:yidx[1],xidx[0]:xidx[1]],
                        cmap=map,
                        origin='lower'
                        )
    im.set_interpolation('nearest')

    if params.has_key('interp'):
        if params['interp'] == 1:
            im.set_interpolation('bilinear')
        else:
            im.set_interpolation('bicubic')



    if params.has_key('title'):
        plt.title(params['title'])

    if params.has_key('xtitle'):
        subplot.set_xlabel(params['xtitle'])

    if params.has_key('ytitle'):
        subplot.set_ylabel(params['ytitle'])

    if params.has_key('aspect'):
        subplot.set_aspect(params['aspect'])

    # set extent after aspect is defined because the option 'auto' does not
    # behave well if the extent is defined.
    im.set_extent((xextent[0],xextent[1],yextent[0],yextent[1]))

    if params.has_key('yrange'):
        print "yrange ",params['yrange']
    
    if params.has_key('legend'):
        cb = plt.colorbar(im)
        if params.has_key('cbtitle'):
            cb.set_label(params['cbtitle'])

    if params.has_key('zrange'):
        zrange = params['zrange']
        im.set_clim(zrange[0],zrange[1])

    if params.has_key('noshow'):
        return subplot
    else:
        plt.show()
        
