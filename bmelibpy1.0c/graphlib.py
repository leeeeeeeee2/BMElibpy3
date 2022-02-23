import matplotlib
matplotlib.use('GTK')
# from matplotlib.matlab import *
from math import *
from scipy import *
from scipy.stats import *
from scipy.interpolate import *
from string import join

from iolib import *
from genlib import *
from graphlib import *
from modellib import *
from statlib import *


def graphlibkwargs():
    '''
    Help on the use of kwargs in BMElibpy graphlib functions

    Dimitri D''Or - 12 November 2004

    kwargs is a KeyWord arguments dictionary. It is frequently used in python to configure optional
    arguments to a function using keywords.

    There are two ways to set kwargs :

    1. Simultaneously with the function call:
       ex. markerplot(c,z,sizelim,zrange=[5,50],markeredgecolor='r',markerfacecolor=[0,1,0])
           where zrange is an optional argument member of the args list
                 markeredgecolor is a key in the kwargs dictionary, with value 'r'
                 markerfacecolor is a key in the kwargs dictionary, with value [0,1,0]

    2. Before the function call:
       ex. kwargs={'markeredgecolor':'r','markerfacecolor':[0,1,0]}
           markerplot(c,z,sizelim,zrange=[5,50],kwargs)
           where zrange is an optional argument member of the args list
                 markeredgecolor is a key in the kwargs dictionary, with value 'r'
                 markerfacecolor is a key in the kwargs dictionary, with value [0,1,0]

    Note the difference in syntax for the kwargs arguments.

    This way of doing is useful because functions such as LEN can be applied to kwargs.
    Keyword arguments can also be added within the function itself. For an example, see in MARKERPLOT (in graphlib) 

    '''


def colorplot(c,z,*args,**kwargs):
    '''
    Translated by Dimitri D''Or - November 2004

    Status: not completed: colormap, ishold

    % colorplot                 - plot of values at 2-D coordinates using colormaps (Jan 1,2001)
    %
    % Plot colored symbols for the values of a vector at a set of two dimensional
    % coordinates. The function uses a colormap such that the color displayed at
    % these coordinates is a function of the corresponding values.
    %
    % SYNTAX :
    %
    % colorplot(c,z,args,kwargs)
    %
    % INPUT :
    %
    % c           n by 2 matrix of coordinates for the locations to be displayed.
    %                    Each line corresponds to the vector of coordinates at a location,
    %                    so the number of columns is equal to two.
    % z           n by 1 column vector of values to be coded as colors.
    % *args              list containing optional arguments consisting here in:
    %                    - map     string  that contains the name of the color map to be used.
    %                                      See the help about graph3d.m for more information on
    %                                      color maps. E.g., map='hot' yields a black-red-yellow-white
    %                                      gradation of colors.
    %                    - zrange  1 by 2  optional vector specifying the minimum and maximum value of
    %                                      the range of z values scaled with the symbol sizes.  
    %                                      The default is zrange=[min(z) max(z)]
    %                                      If zrange is not mentioned then the default values are used.
    % **kwargs           dictionary containing optional arguments for defining plot object properties
    %                    Execute get(H), where H is a plot handle, to see a list of plot object
    %                    properties and their current values. Execute set(H) to see a list of
    %                    plot object properties and legal property values. See also the help
    %                    for plot.m.
    %                    Type help(graphlibkwargs) for detailed explanation on how to use kwargs
    %
    % NOTE :
    % 
    % For example,
    %
    % colorplot(c,z,'hot',[5,50],kwargs)
    %
    % where kwargs={'marker':'^', 'markeredgecolor':'b', 'markersize':10}
    %
    % will plot red triangles with a black border that have a MarkerSize value
    % equal to 10. By default, colorplot(c,z) will use circles with a MarkerSize
    % equal to 10 and with a MarkerEdgeColor equal to the default color.
    '''

    if not args:
        map='hot'

    noptions=len(kwargs)

    if len(args)<=1:
      zrange=[min(z[:]), max(z[:])]

    [n,d]=c.shape
    if d != 2: return 'c must be a n by 2 matrix'
    if prod(z.shape) != z.shape[0]: return 'z must be a column vector'
    if z.shape[0] != n: return 'c and z must have the same number of rows'

    c=take(c,(find(logical_not(isnan(z)))))
    z=take(z,(find(logical_not(isnan(z)))))
    if len(z) == 0: return

#    test = (hold==True)
    minz=zrange[0]
    maxz=zrange[1]

    if maxz==minz:
        return 'At least two values of z must be different'

    colormap(map)
    map=colormap
    nc=colormap.shape[0]

    n=len(z)
    for i in xrange(0,n):
        index=((z[i]-minz)/(maxz-minz))*(nc-1)+1

        if index<1:
            Color=map[0,:]
        elif index>map.shape[0]:
            Color=map[-1,:]
        else:
            indexl=floor(index)
            indexu=ceil(index)
            if indexl==indexu:
                Color=map[index,:]
            else:
                finterp1d=interp1d([indexl, indexu],[map[indexl,0], map[indexu,0]],'linear')
                Color[0]=finterp1d(index)
                finterp1d=interp1d([indexl, indexu],[map[indexl,1], map[indexu,1]],'linear')
                Color[1]=finterp1d(index)
                finterp1d=interp1d([indexl, indexu],[map[indexl,2], map[indexu,2]],'linear')
                Color[2]=finterp1d(index)

    a=plot(c[i,0],c[i,1],'o')
    if kwargs:
        kwargs['MarkerFaceColor']=Color
        set(a,**kwargs)
    else:
        set(a,markersize=10,markerfacecolor=Color)
    hold(True)

    #%
    #% Set the color axis for the colorbar
    #%
    ax=axis
    patch(0,0,0)    #% Somehow this helps so the next line has an effect
    caxis(zrange)
    axis(ax)

##     if test==0:
##         hold(False)

    show()


def histscatterplot(Z,*args):
    '''
    Translated by Dimitri D''Or - November 2004

    Status: approved

    % histscatterplot           - histograms and scatter plots (Jan 1,2001)
    %
    % Plot the histograms for a set of variables and the corresponding
    % scatterplots between variables.
    %
    % SYNTAX :
    %
    % histscatterplot(Z,args)
    %
    % INPUT : 
    %
    % Z          n by nv matrix of values for the different variables, where
    %                    each column corresponds to the values of one variable.
    %                    There must be at least two variables.
    % *args              list containing optional arguments consisting here in:
    %                    - nbins      scalar  number of bins to be used for drawing
    %                                 the histograms. If args is omitted from the
    %                                 input list of variables, the function will use
    %                                 the default number of bins as specified in hist.m.
    '''

    ##%%% Initialize the parameters

    (n,nv)=Z.shape
    if nv==1:
        return 'The Z matrix should at least have two columns'

    ##%%% display the histograms and scatter plots

#    test= (hold==True)
    for i in xrange(0,nv):
        for j in xrange(i,nv):
            subplot(nv,nv,i*nv+j+1)
            if i==j:
                if not args:
                    hist(Z[:,i])
                else:
                    hist(Z[:,i],args[0])
                labels = get(gca(),'xticklabels')
                set(labels,fontsize = 6)
                labels = get(gca(),'yticklabels')
                set(labels,fontsize = 6)
                title(join(['Variable',str(i)]),fontsize=8)
            else:
                plot(Z[:,i],Z[:,j],'.')
                labels = get(gca(),'xticklabels')
                set(labels,fontsize = 6)
                labels = get(gca(),'yticklabels')
                set(labels,fontsize = 6)
                title(join(['Couple',str(i),'-',str(j)]),fontsize=8)

#    if test==0:
#        hold(False)
    show()


def markerplot(c,z,sizelim,*args,**kwargs):
    '''
    Translated by Dimitri D''Or - November 2004

    Status: approved, just the test on hold  status to solve

    % markerplot                - plot of values at 2-D coordinates using markers (Jan 1,2001)
    %
    % Plot the values of a vector at a set of two dimensional coordinates
    % using symbols of varying sizes such that the size of the displayed
    % symbols at these coordinates is a function of the corresponding values. 
    %
    % SYNTAX :
    %
    % markerplot(c,z,sizelim,args,kwargs) 
    %
    % INPUT :
    %
    % c           n by 2 matrix of coordinates for the locations to be displayed.
    %                    Each line corresponds to the vector of coordinates at a
    %                    location, so the number of columns is equal to two.
    % z           n by 1 column vector of values to be coded as markers.
    % sizelim     1 by 2 vector that contains the minimum and maximum values in pixels
    %                    for the size of the symbols to be displayed. The minimum and
    %                    maximum size values are associated with the minimum  and maximum
    %                    values in z, respectively. The size of the symbols for values in
    %                    between these minimum and maximum are obtained by linear interpolation.
    % *args              list containing optional arguments consisting here in:
    %                    - zrange  1 by 2  optional vector specifying the minimum and maximum value of
    %                                      the range of z values scaled with the symbol sizes.  
    %                                      The default is zrange=[min(z) max(z)]
    %                                      If zrange is not mentioned then the default values are used.
    % **kwargs           dictionary containing optional arguments for defining plot object properties
    %                    Execute get(H), where H is a plot handle, to see a list of plot object
    %                    properties and their current values. Execute set(H) to see a list of
    %                    plot object properties and legal property values. See also the help
    %                    for plot.m.
    %                    Type help(graphlibkwargs) for detailed explanation on how to use kwargs
    %
    % NOTE :
    %
    % For example,
    %
    % markerplot(c,z,sizelim,[5,50],kwargs)
    % 
    % where sizelim=[5, 20]
    %       kwargs={marker='^', markeredgecolor='b', markerfacecolor=[1 0 0]}
    %
    % will plot red triangles with a blue border that have a 
    % MarkerSize value between 5 and 20 pixels. By default, 
    % markerplot(c,z,sizelim) will use disks with the default 
    % properties for plot.m.
    '''

    noptions=len(kwargs)

    if not args:
        zrange = [min(z[:]), max(z[:])]
    else:
        zrange = [args[0][0], args[0][1]]

    (n,d)=c.shape
    if d != 2:  return 'c must be a n by 2 matrix'
    if prod(z.shape) != z.shape[0]: return 'z must be a vector'
    if z.shape[0] != n: return 'c and z must have the same number of rows'

    c=take(c,(find(logical_not(isnan(z)))))
    z=take(z,(find(logical_not(isnan(z)))))
    if len(z) == 0: return

    test = (hold==True)
    minz=zrange[0]
    maxz=zrange[1]

    if maxz==minz:
        return 'At least two values of z must be different'

    n=len(z)
    finterp1d=interp1d([minz, maxz],sizelim,'linear')
    for i in xrange(0,n):
        if z[i]<minz:
            s=[sizelim[0]]
        elif z[i]>maxz:
            s=[sizelim[1]]
        else:
            s=finterp1d(z[i])
        a=plot(c[i,0],c[i,1],'o')
        if kwargs:
            kwargs['markersize']=s[0]
            set(a,**kwargs)
        else:
            set(a,markersize=s[0])
        hold(True)

    if test==0:
        hold(False)

    show()


def poleplot(c,z,*args,**kwargs):
    '''
    Translated by Dimitri D''Or - November 2004

    Status: not completed (ishold,plot3), not tested

    % poleplot                  - 3-D perspective plot at 2-D coordinates using poles (Jan 1,2001)
    %
    % Plot the values of a vector at a set of two dimensional coordinates 
    % in a three dimensional perspective using poles of various lengths.
    % The lengths of the poles are proportional to the corresponding values.
    % Positive values have upward orientated poles, negative values have
    % downward orientated poles.
    %
    % SYNTAX :
    %
    % poleplot(c,z,args,kwargs) 
    %
    % INPUT :
    %
    % c           n by 2 matrix of coordinates for the locations to be displayed.
    %                    Each line corresponds to the vector of coordinates at a location,
    %                    so the number of columns is equal to two.
    % z           n by 1 column vector of values to be coded as poles.
    % *args              list containing optional arguments consisting here in:
    %                    - z0  scalar  optional value specifying the basic level for the poles
    %                                  If z0 is not mentioned, then z0 is set to 0.
    % **kwargs           dictionary containing optional arguments for defining plot object properties
    %                    Execute get(H), where H is a plot handle, to see a list of plot object
    %                    properties and their current values. Execute set(H) to see a list of
    %                    plot object properties and legal property values. See also the help
    %                    for plot.m.
    %                    Type help(graphlibkwargs) for detailed explanation on how to use kwargs
    %
    % NOTE :
    %
    % For example,
    %
    % poleplot(c,z,args,kwargs)
    %
    % where args=-2
    %       kwargs={marker='^', markeredgecolor='b', markerfacecolor=[1 0 0]}
    %
    % will plot red triangles with a blue border. The basis level is set at -2. By default,
    % poleplot(c,z) will use points with the default properties
    % for plot3.m.
    '''

    noptions=len(kwargs)

    if not args:
        z0=0

#    test=(ishold==True)
    x=c[:,0]
    y=c[:,1]
    xmin=min(x)
    xmax=max(x)
    ymin=min(y)
    ymax=max(y)

    for i in xrange(0,max(z.shape)):
        a=plot3(x[i],y[i],z[i],'o')
        if kwargs:
            set(a,**kwargs)
        LineColor=get(a,'MarkerEdgeColor')
        b=plot3(array([x[i],x[i]]),array([y[i],y[i]]),array([z[i],z0]))
        hold(True)
        if not isinstance(LineColor,str):
            set(b,color=LineColor)
    grid(True)

#    if test==0:
#        hold(False)

    show()

def probaplot(softpdftype,nl,limi,probdens,*args):
    '''
    Translated by Dimitri D''Or - November 2004

    Status: completed, not tested

    % probaplot                 - Plots the probabilistic data (Jan 1, 2001) 
    %
    % [h]=probaplot(softpdftype,nl,limi,probdens) plots the probabilistic
    %     soft data and returns the handle h. See help probasyntax for 
    %     explanation of the soft data defined by softpdftype, nl, limi
    %     and probdens.
    % 
    % [h]=probaplot(softpdftype,nl,limi,probdens,args) plots the probabilistic
    %     soft data using optional arguments args:
    %     *args    list containing optional arguments consisting here in:
    %              - S     linetype indication.
    %                      character string made from one element
    %                      from any or all the following 3 colunms:
    % 
    %                      y     yellow        .     point              -     solid
    %                      m     magenta       o     circle             :     dotted
    %                      c     cyan          x     x-mark             -.    dashdot 
    %                      r     red           +     plus               --    dashed   
    %                      g     green         *     star
    %                      b     blue          s     square
    %                      w     white         d     diamond
    %                      k     black         v     triangle (down)
    %                                          ^     triangle (up)
    %                                          <     triangle (left)
    %                                          >     triangle (right)
    %                                          p     pentagram
    %                                          h     hexagram
    %
    %              - idx   vector of index for plotting only the associated soft pdf
    % 
    % OUTPUT :
    % h        vector of length ns of handles to the lines representing the 
    %          soft data at each soft data points, where ns=length(nl)
    %
    '''

    if not args:
        S='b-'
        idx=arrayrange(0,len(nl))
        
    if len(args) == 1:
        idx=arrayrange(0,len(nl))

    ns=len(nl)

    for i in xrange(0,len(idx)):
        ist=idx[i]
        if softpdftype == 1:
            x=kron(limi[ist,0:nl[ist]],array([1, 1]))
            y=concatenate((array([0]), kron(probdens[ist,0:nl[ist]-1],array([1, 1])),array([0])),1)
        elif softpdftype == 2:
            x=concatenate((limi[ist,0], limi[ist,0:nl[ist]], limi[ist,nl[ist]]),1)
            y=concatenate((array([0]), probdens[ist,0:nl[ist]], array([0])),1)
        elif softpdftype == 3:
            x=kron(arrayrange(limi[ist,0],limi[ist,2],limi[ist,1]),array([1, 1]))
            y=concatenate((array([0]), kron(probdens[ist,0:nl[ist]-1],array([1, 1])), array([0])),1)
        elif softpdftype == 4:
            x=concatenate((array([limi[ist,0]]), arrayrange(limi[ist,0],limi[ist,2],limi[ist,1]), array([limi[ist,2]])),1)
            y=concatenate((array([0]), probdens[ist,0:nl[ist]], array([0])),1)
        hold(True)
        h[i]=plot(x,y,S)

    return h


def stringplot(c,s,**kwargs):
    '''
    Translated by Dimitri D''Or - November 2004

    Status: approved

    % stringplot                - plot of text strings at 2-D coordinates (Jan 1,2001)
    %
    % Plot arbitrary text strings at a set of two dimensional coordinates.
    %
    % SYNTAX :
    %
    % stringplot(c,str,Property,Value); 
    %
    % INPUT :
    %
    % c           n by 2 matrix of coordinates for the locations to be displayed.
    %                    Each line corresponds to the vector of coordinates at a
    %                    location, so the number of columns is equal to two.
    % s           n by 1 list of strings having as strings as there are values in z.
    %                    Each  string contains the text that must be
    %                    displayed at the corresponding coordinates specified in c.
    % **kwargs           dictionary containing optional arguments for defining plot object properties
    %                    Execute get(H), where H is a plot handle, to see a list of plot object
    %                    properties and their current values. Execute set(H) to see a list of
    %                    plot object properties and legal property values. See also the help
    %                    for plot.m.
    %                    Type help(graphlibkwargs) for detailed explanation on how to use kwargs
    %
    % NOTE :
    %
    % For example,
    %
    % stringplot(c,str,kwargs)
    %
    % where kwargs={'fontfize':30,'Color':[1 0 0]}
    %
    % will display the text with FontSize=30 and Color=[1 0 0].
    % By default, stringplot(c,str) will use the default properties
    % used by text.m.
    '''

    noptions=len(kwargs)

    n=len(s)
    for i in xrange(0,n):
#        plot(c[i,0],c[i,1],'.w')
        a=text(c[i,0],c[i,1],s[i])
        hold(True)
        if kwargs:
            set(a,**kwargs)

    minc1=floor(min(c[:,0]))
    minc2=floor(min(c[:,1]))
    maxc1=ceil(max(c[:,0]))
    maxc2=ceil(max(c[:,1]))

    axis([minc1, maxc1, minc2, maxc2])

    show()


def valplot(c,z,**kwargs):
    '''
    Translated by Dimitri D''Or - November 2004

    Status: approved

    % valplot                   - plot of values at 2-D coordinates using texts (Jan 1,2001)
    %
    % Plot the values of a vector at a set of two dimensional coordinates
    % by writing these values as text strings at the corresponding
    % coordinates. 
    %
    % SYNTAX :
    %
    % valplot(c,z,Property,Value); 
    %
    % INPUT :
    %
    % c           n by 2 matrix of coordinates for the locations to be displayed.
    %                    Each line corresponds to the vector of coordinates at a
    %                    location, so the number of columns is equal to two.
    % z           n by 1 column vector of values to be coded as strings.
    % **kwargs           dictionary containing optional arguments for defining plot object properties
    %                    Execute get(H), where H is a plot handle, to see a list of plot object
    %                    properties and their current values. Execute set(H) to see a list of
    %                    plot object properties and legal property values. See also the help
    %                    for plot.m.
    %                    Type help(graphlibkwargs) for detailed explanation on how to use kwargs
    %
    % NOTE :
    %
    % For example,
    %
    % valplot(c,z,kwargs)
    %
    % where kwargs={'fontfize':30,'Color':[1 0 0]}
    %
    % will display the text with FontSize=30 and Color=[1 0 0].
    % By default, valplot(c,z) will use the default properties
    % for text.m.
    '''
    
    noptions=len(kwargs)

    n=len(z)
    for i in xrange(0,n):
        a=text(c[i,0],c[i,1],str(z[i]))
        if kwargs:
            set(a,**kwargs)

    minc1=floor(min(c[:,0]))
    minc2=floor(min(c[:,1]))
    maxc1=ceil(max(c[:,0]))
    maxc2=ceil(max(c[:,1]))

    axis([minc1, maxc1, minc2, maxc2])

    show()
