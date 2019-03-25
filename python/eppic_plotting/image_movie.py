

import os
import tempfile
import shutil
from image_plot import image_plot
import matplotlib.pyplot as plt
try:
    from matplotlib.animation import FuncAnimation
except:
    print '''
    image_movie requires matplotlib 1.1.0 or greater. It expects the animation
    module.
    '''
    raise 


def image_movie(movie_frames,selfnormalize=False,delay=25,movie_prefix='movie',
                show_movie=True,save=False,print_time=True,timescale=1,
                timetitle="Time Step %4d",format='mpg',delTmpFiles=True,
                *extra_args,**extra_keywords):
    '''image_movie
    Calls image_plot for each 2D slice of the 3D array movie_frames.

    Takes all arguments image_plot takes.
    @see image_plot

    Arguments:
    ----------
    movie_frames  - An (exactly) 3D array, where the first two dimensions are
                    the images and the last dimension is the time. 

    Keywords:
    ---------
    selfnormalize - Don\'t set the zrange to the max and min for all time. 
                    image_plot\'s zrange keyword takes precedence. 

    delay         - The number of frames convert should delay, the default is
                    25. 

    movie_prefix  - A string used to name the movie, the default is 'movie'

    show_movie    - Use the animate program to show the movie once created.
                    The default is True.

    save          - Copy the movie to the current working directory. The
                    default is False.

    print_time    - Print a subtitle (bottom of the figure) with the
                    time step info; the default is True.

    timescale     - Scale the integer time step; the default is 1.

    timetitle     - Format the time title, use string modifiers such as \%f
                    to incorporate the time step.
                    
    format        - Set movie output format. Defaults to "mpg". "gif", "avi"
                    and "flv" are other options.

    delTmpFiles   - Delete temporary files; by default this is True, if False,
                    then the routine prints the path to the tmp files.

    \todo theres no test for an invalid format
    '''

    # test movie_frames
    dims = movie_frames.shape
    if (len(dims) != 3):
        raise "movie_frames array does not have the correct dimensions."
    nt=dims[2]
    if (nt < 2):
        raise "not enough movie frames to make a movie."

    # setup directory
    tmp_outdir = tempfile.mkdtemp()

    # set range
    if not selfnormalize:
        if not extra_keywords.has_key('zrange'):
            extra_keywords['zrange'] = [
                movie_frames.min(),
                movie_frames.max()
                ]


# in development: 

#     fig = plt.figure()
#     subplot = fig.add_subplot(1,1,1)
    
#     it = 0
    
#     def frameUpdater()
#         '''A routine to plot one frame.'''
#         subplot.cla()
#         if print_time:
#             fig.suptitle(timetitle%(it*timescale),x=.5,y=.01,va='baseline')
#         image_plot(frame,subplot=subplot,noshow=1,**extra_keywords)
#         it=(it+1)%nt
#         #pass

#     def frameGetter():
#         '''A routine to return the current frame'''
#         pass

#     ani = animation()

    # loop over frame plotting things, writing file
    for it in range(0,nt):
        plt.clf()
        subplot = fig.add_subplot(1,1,1)
        if print_time:
            fig.suptitle(timetitle%(it*timescale),x=.5,y=.01,va='baseline')
        frame = movie_frames[:,:,it]
        image_plot(frame,subplot=subplot,noshow=1,**extra_keywords)
        plt.savefig(tmp_outdir+"/movie%05d.png"%it,
                    dpi=72,papertype='letter',format='png')

    # use convert to make a movie
    ## convert -delay 25 *png movie.gif
    if format== 'gif':
        command=[
            'convert',
            '-delay',
            "%d"%delay,
            tmp_outdir+'/*png',
            tmp_outdir+'/movie.gif'
            ]
    if format=='flv':
        command=[
            'mencoder',
            '\'mf://'+tmp_outdir+'/*.png\'',
            '-mf type=png:fps=%d'%(10),
            '-ovc lavc',
            '-oac copy',
            '-of lavf',
            '-really-quiet',
            '-lavcopts vcodec=flv:vbitrate=500:trell:v4mv:last_pred=3',
            '-lavfopts i_certify_that_my_video_stream_does_not_use_b_frames',
            '-o '+tmp_outdir+'/movie.flv'
            ]
    if format=='avi':
        command=[
            'mencoder',
            '\'mf://'+tmp_outdir+'/*.png\'',
            '-mf type=png:fps=%d'%(60./delay),
            '-ovc lavc',
            '-oac copy',
            '-really-quiet',
            #'-of lavf',
            '-lavcopts vcodec=wmv2',#:vbitrate=500:trell:v4mv:last_pred=3',
            '-lavfopts i_certify_that_my_video_stream_does_not_use_b_frames',
            '-o '+tmp_outdir+'/movie.avi'
            ]
    if format=='mpg':
        command=[
            'mencoder',
            '\'mf://'+tmp_outdir+'/*.png\'',
            '-mf type=png:fps=%d'%(60./delay),
            '-ovc lavc',
            '-really-quiet',
            '-oac copy',
            #'-of lavf',
            '-lavcopts vcodec=mpeg4',#:vbitrate=500:trell:v4mv:last_pred=3',
            '-lavfopts i_certify_that_my_video_stream_does_not_use_b_frames',
            '-o '+tmp_outdir+'/movie.mpg'
            ]

    print "movie command:\n"
    print " ".join(command)
    os.system(" ".join(command))

    if show_movie:
        command = ['mplayer',tmp_outdir+'/movie.'+format]
        if format=='gif':
            command=['animate',tmp_outdir+'/movie.gif']

        os.system(" ".join(command))
        
    # copy movie to current directory
    if save:
        try:
            shutil.move(tmp_outdir+'/movie.'+format,"./"+movie_prefix+'.'+format)
#             print "moved %s to %s"%(tmp_outdir+'/movie.'+format,
#                                     "./"+movie_prefix+'.'+format)
        except:
            print "Failed to move %s to %s"%(tmp_outdir+'/movie.'+format,
                                             "./"+movie_prefix+'.'+format)
                    
    # Clean up the directory
    if delTmpFiles:
        shutil.rmtree(tmp_outdir,True)
    else:
        print tmp_outdir," is the tmp directory..."



    
