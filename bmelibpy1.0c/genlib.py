#coding=utf-8
## import iolib
## import modellib
## import statlib
## import BMEmatlab

## import matplotlib
## matplotlib.use('GTK')
## import matplotlib.matlab
## import scipy

import matplotlib
matplotlib.use('GTK')
# from matplotlib.matlab import *
from scipy import *
from scipy.stats import *
from scipy.interpolate import *

from iolib import *
from genlib import *
from graphlib import *
from modellib import *
from statlib import *
import BMEmatlab


def aniso2iso(c,angle,ratio):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: debugged, not tested
    
    % aniso2iso                 - convert anisotropic to isotropic coordinates (Jan 1,2001)
    %
    % Transform a set of two dimensional or three dimensional coordinates
    % using rotations and dilatations of the axes, in order to map an
    % anisotropic space into an isotropic one. The geometric anisotropy
    % is characterized by the angle(s) of the principal axis of the ellipse
    % (in 2D) or ellipsoid (in 3D), and by the ratio(s) of the principal axis
    % length by the other axes lengths. Using this function, an ellipse
    % (ellipsoid) is thus mapped into a circle (sphere) having as radius the
    % length of the principal axis. The transformation consist in a rotation
    % of the axes followed by a dilatation. 
    % 
    % SYNTAX :
    %
    % [ciso]=aniso2iso(c,angle,ratio); 
    %
    % INPUT :
    %
    % c       n by d   matrix of coordinates for the locations in the anisotropic
    %                  space. A line corresponds to the vector of coordinates at a
    %                  location, so that the number of columns corresponds to the
    %                  dimension of the space. Only two dimensional or three dimensional
    %                  space coordinates can be processed by this function.
    % angle   1 by d-1 vector of angle values that characterize the anisotropy. 
    %                  In a two dimensional space, angle is the trigonometric angle
    %                  between the horizontal axis and the principal axis of the
    %                  ellipse. In a three dimensional space, spherical coordinates
    %                  are used, such that angle(1) is the horizontal trigonometric
    %                  angle and angle(2) is the vertical trigonometric angle for the
    %                  principal axis of the ellipsoid. All the angles are measured
    %                  counterclockwise in degrees and are between -90� and 90�.
    % ratio   1 by d-1 vector that characterize the ratio for the length of the axes
    %                  for the ellipse (in 2D) or ellipsoid (in 3D). In a two dimensional
    %                  space, ratio is the length of the principal axis of the ellipse
    %                  divided by the length of the secondary axis, so that ratio>1. In a
    %                  three dimensional space, ratio(1) is the length of the principal
    %                  axis of the ellipsoid divided by the length of the second axis, 
    %                  whereas ratio(2) is length of the principal axis of the ellipsoid
    %                  divided by the length of the third axis, so that ratio(1)>1 and
    %                  ratio(2)>1.
    %
    % OUTPUT :
    %
    % ciso    n by d   matrix having the same size as c, that gives the new coordinates
    %                  in the isotropic space.
    %
    % NOTE :
    %
    % It is possible to specify an additional index vector, taking integer values from 1
    % to nv. The values in index specify which one of the nv variable is known at each one
    % of the corresponding coordinates. The c matrix of coordinates and the index vector
    % are then grouped together using the MATLAB cell array notation, so that c={c,index}.
    % This allows to perform the same coordinate transformation at once on a set of possibly
    % different variables. The output variable ciso is then also a cell array that contains
    % both the new matrix of coordinates and the index vector.
    '''

    ##%%%%%% Determine the dimension of the space and set angle

    islist=isinstance(c,list)
    if islist==True:
      d=c[0].shape[1]
    else:
      d=c.shape[1]
    angle=angle*2*pi/360
    angle=-angle

    ##%%%%%% When d<2 or d>3, error

    if (d<2)or(d>3):
      return 'aniso2iso.m requires coordinates in a 2D or 3D space'

    ##%%%%%% Case for d=2

    if d==2:

        R=array([ cos(angle),  sin(angle),
                 -sin(angle),  cos(angle)])

        if islist==True:
            ciso=c
            ciso[0]=c[0]*R
            ciso[0][:,1]=ciso[0][:,1]*ratio
        else:
            ciso=c*R
            ciso[:,1]=ciso[:,1]*ratio

    ##%%%%%% Case for d=3

    if d==3:
        phi=angle[0]
        teta=angle[1]
        ratioy=ratio[0]
        ratioz=ratio[1]

        R1=array([ cos(phi),   sin(phi),  0,
                  -sin(phi),   cos(phi),  0,
                   0,          0,         1])
        R2=array([ cos(teta),  0,  sin(teta),
                   0,          1,  0,
                  -sin(teta),  0,  cos(teta)])
        R=R1*R2

        if islist==True:
            ciso=c
            ciso[0]=c[0]*R
            ciso[0][:,1]=ciso[0][:,1]*ratioy
            ciso[0][:,2]=ciso[0][:,2]*ratioz
        else:
            ciso=c*R
            ciso[:,1]=ciso[:,1]*ratioy
            ciso[:,2]=ciso[:,2]*ratioz

    return ciso

def avedupli2D(x,y,z):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: Approved

    % avedupli2D               - averages duplicate values for 2D coordinate x-y
    %
    %  all the duplicate values (with same x-y coordinates) are averaged.
    %  The returned coordinates xu and yu are unique, while za is the corresponding
    %  averaged value.  The index j is such that xu=x(j) and yu=y(j)
    %
    %  SYNTAX :
    % 
    %  [xu,yu,za,j] = avedupli2D(x,y,z);
    % 
    %  INPUT :
    %
    %  x   vector      vector of x coordinates
    %  y   vector      vector of y coordinates
    %  z   vector      vector of z values
    %
    %  OUTPUT :
    %
    %  xu and yu  vectors    vector of x-y coordinates such that (xu,yu) are a set of unique points
    %  za         vector     vector of corresponding z-values obtained by averaging the duplicate points
    %  j          vector     index such that xu=x(j) and yu=y(j)
    '''

    [msg,x,y,z,dumb1,dumb2] = BMEmatlab.xyzchk(x,y,z)
    if len(msg)!=0:
        return msg
    
    #% Need x,y and z to be column vectors

    if len(x.shape)==1:
        x = x[:,NewAxis]
        x=x.astype(Float32)
    
    if len(y.shape)==1:
        y = y[:,NewAxis]
        y=y.astype(Float32)
    
    if len(z.shape)==1:
        z = z[:,NewAxis]
        z = z.astype(Float32)

    sz = prod(x.shape)
    j = arange(0,sz)[:,NewAxis]
    j=j.astype(Float32)

    (sxyzj,dumb) = BMEmatlab.sortrows(concatenate((x,y,z,j),1),array([1,0]))
    xs = sxyzj[:,0]
    ys = sxyzj[:,1]
    zs = sxyzj[:,2]
    js = sxyzj[:,3]
    ind = concatenate((array([0]), (ys[1:] == ys[:-1] and xs[1:] == xs[:-1]), array([0])),0)
    if sum(ind) > 0:
        #%disp('Duplicate x-y data points detected: using average of the z values');
        fs = find(logical_and((ind[:-1] == 0),(ind[1:] == 1)))
        fe = find(logical_and((ind[:-1] == 1), (ind[1:] == 0))) 
        for i in xrange(0,len(fs)):
            #% averaging z values
            zs[fe[i]] = mean(zs[fs[i]:fe[i]+1])
        xs = take(xs,find(ind[1:]==0))
        ys = take(ys,find(ind[1:]==0))
        zs = take(zs,find(ind[1:]==0))
        js = take(js,find(ind[1:]==0))
        szs = len(xs)
        xs = reshape(xs,(szs,1))
        ys = reshape(ys,(szs,1))
        zs = reshape(zs,(szs,1))
        js = reshape(js,(szs,1))
    (xyzu,dumb) = BMEmatlab.sortrows(concatenate((xs,ys,zs,js),1),array([3]))
    xu=xyzu[:,0]
    yu=xyzu[:,1]
    za=xyzu[:,2]
    j=sort(js)

    return xu,yu,za,j

def avedupli(p,z):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: Approved
    
    % avedupli               - averages duplicate values
    %
    %  all the duplicate values (with same s/t p coordinates) are averaged.
    %  The returned coordinates pu are unique, while za is the corresponding
    %  averaged value.  The index i is such that pu=p(i,:)
    %
    %  SYNTAX :
    % 
    %  [pu,za,i] = avedupli(p,z);
    % 
    %  INPUT :
    %
    %  p          n by d     matrix of the coordinates of n points in a s/t space of dimension d
    %  z          n by 1     vector of z values
    %
    %  OUTPUT :
    %
    %  pu         nu by d    matrix of the coordinates of nu unique points in a space of dimension d, 
    %                           note that nu<=n
    %  za         nu by 1    vector of corresponding z-values obtained by averaging the duplicate points
    %  i          nu by 1    index such that pu=p(i,:)
    '''

    [n,d]=p.shape
    if z.shape[0]!=n:
        return 'z must have the same number of rows as p'
    if n==0:
        pu=p
        za=z
        i=[]
        return pu,za,i
    if len(z.shape)!=1:
        if z.shape[1]!=1:
            return 'z must have one column'

    z=ravel(z)

    if isdupli(p)!=True:
      pu=p
      za=z
      i=range(0,len(z[:]))
    else:
      [iu,ir]=finddupli(p)
      pu=take(p,iu)
      za=take(z,iu)
      i=iu
      nu=len(iu)
      for k in xrange(0,len(ir)):
          i=concatenate((i,array([ir[k][0]])))
          za=concatenate((za,array([mean(take(z,ir[k]))])))
      pu=take(p,i)

    return pu,za,i

def creategrid(minc, dc, nc):
    """
    
    Translated by Dimitri D''Or - October 2004

    Status: approved
    
    % creategrid                - create a regular grid in 1D, 2D or 3D (Jan 1,2001)
    %
    % Create a regular grid of Cartesian coordinates up to a three dimensional space.
    %
    % SYNTAX :
    %
    % [c]=creategrid(minc,dc,nc);
    %
    % INPUT :
    %
    % minc  1 by d   vector of coordinates for the origin of the grid. The number
    %                of elements is equal to the dimension of the space.
    % dc    1 by d   vector of values for the distances between two grid points
    %                along each axis. The number of elements is the same as for minc.
    % nc    1 by d   vector of values for the number of grid points along each axis.
    %                The number of elements is the same as for minc.
    %
    % OUTPUT :
    %
    % c     n by d   matrix of coordinates for the points of the grid. Each line
    %                corresponds to the vector of coordinates at one of the grid point.
    %                The coordinates in c are sorted in ascending order according to the
    %                successive space axes.
    """

    #print "Create grid started"
    d = minc.shape[0] 
    if d==1 :
        xMat=arrayrange(minc[0],minc[0]+nc[0]*dc[0],dc[0])
        return xMat
    elif d==2:        
        x = arrayrange(minc[0], minc[0]+ nc[0]*dc[0], dc[0])
        y = arrayrange(minc[1], minc[1]+ nc[1]*dc[1], dc[1])
        [c1,c2]=meshgrid(x,y)
        return transpose(array([ravel(c1),ravel(c2)]))
    elif d==3 :
        x = arrayrange(minc[0], minc[0]+ nc[0]*dc[0], dc[0])
        y = arrayrange(minc[1], minc[1]+ nc[1]*dc[1], dc[1])
        z = arrayrange(minc[2], minc[2]+ nc[2]*dc[2], dc[2])
        [c1,c2,c3]=meshgrid(x,y,z)
        return transpose(array([ravel(c1),ravel(c2),ravel(c3)]))
        
def testCreateGrid():
    print "1D Grid"
    basicCoord = array([0])
    step = array([1])
    count = array([10])
    c = creategrid(basicCoord, step, count)
    #print c
    print c.shape
    print "2D Grid"
    basicCoord = array([0,0])
    step = array([1,1])
    count = array([10,10])
    c = creategrid(basicCoord, step, count)
    #print c
    print c.shape
    print "3D Grid"
    basicCoord = array([0,0,0])
    step = array([1,1,1])
    count = array([2,2,2])
    c = creategrid(basicCoord, step, count)
    print c
    print c.shape
    print "Finished"

    

def col2mat(c,z):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: debugged, not tested
    
    % col2mat                   - columnwise defined gridded variable toward matrix (Jan 1,2001)
    %
    % Convert, e.g. the format of the data output from invdist.m and
    % kernelsmoothing.m, so that the converted data can be directly
    % displayed using graphical functions like mesh.m or pcolor.m. The
    % col2mat.m function converts a columnwise defined gridded variable
    % at a set of planar coordinates into a matrix of values and two
    % separate vectors of coordinates along each axis.
    %
    % SYNTAX :
    %
    % [c1,c2,Z]=col2mat(c,z);
    %
    % INPUT :
    %
    % c     n by 2    matrix of coordinates for the nodes of a rectangular grid
    %                 in a two dimensional space. Each line in c corresponds to
    %                 the coordinates of a grid node, so the number of columns is
    %                 equal to 2, and the number of lines n is equal to n1*n2,
    %                 where n1 and n2 are the number of columns and lines for the
    %                 rectangular grid, respectively.
    % z     n by 1    vector of values at the coordinates specified in c.
    %
    % OUTPUT :
    %
    % c1    n1 by 1   vector of coordinates for the grid nodes along the first axis.
    % c2    n2 by 1   vector of coordinates for the grid nodes along the second axis.
    % Z     n1 by n2  matrix of values that correspond to the values in z.
    %
    % NOTE :
    %
    % It is not needed that the values in c are sorted in ascending order, as
    % the function is taking care of it, but the coordinates for all the points
    % of the rectangular grid should appear once and only once in the c matrix.
    % There is no constraint about the regularity of the spacing for the nodes
    % of the grid along the axes, as long as the grid is rectangular.
    '''
    
    c=c.astype(Float32)
    z=z.astype(Float32)
    cz=concatenate((c,z),1)
    S=sortrows(cz,array([0,1]))    #%%% sort jointly the coordinates and values
    n1=sum(S[0,0]==S[:,0])         #%%% compute the number of columns in the grid
    n2=sum(S[0,1]==S[:,1])         #%%% compute the number of lines in the grid
    Z=reshape(S[:,2],(n1,n2))      #%%% reshape the z vector accordingly
    c1=reshape(S[:,0],(n2,n1))     #%%% extract the coordinates along the lines
    c1=c1[:,0]
    c2=reshape(S[:,1],(n2,n1))     #%%% extract the coordinates along the columns
    c2=c2[0,:]

    return c1,c2,Z



def combinedupli(ch,zh,method='ave',cs=array([]),softpdftype=2,nl=array([]),limi=array([]),probdens=array([])):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: not completed: dependencied to complete : softpdftypeCheckArgs, proba2probdens, probacat,
                                                      probaStudentT, probasplit, probacombinedupli
    
    % combinedupli               - combine duplicate hard and soft data
    %
    %  Each set of duplicated hard and soft data (with same space/time
    %  coordinates) are combined into a hard or soft datum.
    %  First the hard data is combined by replacing each set of duplicate hard 
    %  data by either a unique hard datum (such as the average, minimum, maximum 
    %  or quantile of the set of duplicated hard values), or by a soft pdf (such
    %  as a student T or the histogram of the duplicated hard data).
    %  Then if there are soft data, each duplicated set containing one hard datum 
    %  and one or more soft data is replaced by the hard datum, and each set of 
    %  duplicated soft data is combined using probacombinedupli.m
    %  The returned coordinates chC and csC do not have any duplicated data point.
    %
    %  SYNTAX :
    % 
    %  [chC,zhC,csC,softpdftypeC,nlC,limiC,probdensC] = combinedupli(ch,zh,method,cs,softpdftype,nl,limi,probdens)
    % 
    %  INPUT :
    %
    %  ch          nh by d    matrix of the coordinates of nh hard data points 
    %                         in a space of dimension d
    %  zh          n by 1     vector of zh hard data values
    %  method      char array or numeric value, defining the method to combine
    %                         each set of duplicated hard data as follow
    %              'ave'      - average value
    %              'min'      - minimum value 
    %              'max'      - maximum value 
    %              q          - where 0<q<1, q-quantile of each set of duplicated hard data
    %              'studentT' - creates a student T soft pdf for each set of 
    %                           duplicated hard data using probaStudentT.m
    %              'hist'     - creates a histogram soft pdf for each set of 
    %                           duplicated hard data using histscaled.m
    %              default is 'ave', which is the same as using avedupli.m
    %  cs          ns by d    matrix of the coordinates of ns soft data points 
    %                         in a space of dimension d
    %                         default is []
    %  softpdftype scalar     indicates the type of soft pdf representing the  
    %                         probabilitic soft data at the coordinates specified in cs.  
    %                         softpdftype may take value 1, 2, 3 or 4, as follow:
    %                         1 for Histogram, 2 for Linear, 3 for Grid histogram, 
    %                         and 4 for Grid Linear. (see probasyntax for more explanations)
    %                         default is 2
    %   nl          ns by 1   vector of the number of interval limits. nl(i) is the number  
    %                         of interval limits used to define the soft pdf for soft data 
    %                         point i. (see probasyntax for more explanations)
    %                         default is []
    %   limi        ns by l   matrix of interval limits, where l is equal to
    %                         either max(nl) or 3 depending of the softpdftype.
    %                         limi(i,:) are the limits of intervals for the i-th 
    %                         soft data. (see probasyntax for more explanations)
    %                         default is []
    %   probdens    ns by p   matrix of probability density values, where p is 
    %                         equal to either max(nl)-1 or max(nl), depending on the 
    %                         softpdftype. probdens(i,:) are the values of the probability 
    %                         density corresponding to the intervals for the i-th soft data 
    %                         defined in limi(i,:). (see probasyntax for more explanations)
    %                         default is []
    %
    %  OUTPUT :
    %
    %   chC          nhC by d  matrix of the coordinates of nhC unique hard data points 
    %                          in a space of dimension d. Note that nhC <= nh
    %   zhC          nhC by 1  vector of hard data values at chC
    %   csC          nsC by d  matrix of the coordinates of nsC unique soft data points 
    %                          in a space of dimension d. Note that (nhC+nsC) <= (nC+nC) 
    %   softpdftypeC scalar    indicates the type of soft pdf of combined soft data
    %   nlC          nsC by 1  vector of the number of interval limits.
    %   limiC        nsC by l  matrix of interval limits
    %   probdensC    nsC by p  matrix of probability density values
    %
    % NOTE:
    %
    % softpdftype must be equal to 2 when method is equal to 'studentT' or 'hist'
    %
    % [chC,zhC]=combinedupli(ch,zh,'ave') returns the same output variables as
    % [pu,za] = avedupli(ch,zh)
    '''

    ## if nargin<3, method='ave'; end
    ## if nargin<4, cs=[]; end
    ## if nargin<5, softpdftype=2; end
    ## if nargin<6, nl=[]; end
    ## if nargin<7, limi=[]; end
    ## if nargin<8, probdens=[]; end

    ## nbins=20;  % number of bins used when method is equal to 'hist'

    ## if ~isnumeric(ch) | ~isnumeric(cs), error('ch and cs must be numeric');  end;
    ## [nh,d]=size(ch);
    ## if size(zh,1)~=nh, error('zh must have the same number of rows as ch'); end
    ## if nh~=0 & size(zh,2)~=1, error('zh must have one column'); end
    ## if ~isnumeric(method) & ~ischar(method), error('method must be numeric or character'); end
    ## if isnumeric(method)
    ##   if length(method)~=1, error('when method is numeric, it must have a length of 1');
    ##   elseif ~(0<method & method<1), error('when method is numeric, it must be such that 0<method<1');
    ##   else 
    ##     P=method; 
    ##     method='quantile';
    ##   end
    ## end
    ## [ns]=size(cs,1);
    ## if ns~=0 & nh~=0 & size(cs,2)~=d, error('cs must have the same number of columns as ch'); end
    ## if isempty(softpdftype), softpdftype==2; end;
    ## softpdftypeCheckArgs(softpdftype,nl,limi,probdens);
    ## if size(limi,1)~=ns, error('limi must have the same number of rows as cs'); end
    ## if size(probdens,1)~=ns, error('probdens must have the same number of rows as cs'); end
    ## switch method
    ##   case {'ave','min','max','quantile'}
    ##     methodtype='hard';
    ##   case 'studentT'
    ##     methodtype='soft';
    ##     if softpdftype~=2, 
    ##       error('When method is equal to ''studentT'' then softpdftype must be equal to 2'); 
    ##     end
    ##   case 'hist'
    ##     methodtype='soft';
    ##     if softpdftype~=2, 
    ##       error('When method is equal to ''hist'' then softpdftype must be equal to 2'); 
    ##     end
    ##   otherwise
    ##     error('Bad value for method');
    ## end

    ## csC=[];
    ## softpdftypeC=softpdftype;
    ## nlC=[];
    ## limiC=[];
    ## probdensC=[];

    ## %%%%%% Combine the hard data
    ## if nh==0 | ~isdupli(ch)
    ##   chC=ch;
    ##   zhC=zh;
    ## else
    ##   [iu,ir]=finddupli(ch);
    ##   chC=ch(iu,:);
    ##   zhC=zh(iu);
    ##   nhC=length(iu);
    ##   nsC=0;  
    ##   for k=1:length(ir)
    ##     switch methodtype
    ##       case 'hard'
    ##         nhC=nhC+1;
    ##         chC(nhC,:)=ch(ir{k}(1),:);
    ##         switch method
    ##           case 'ave', zhC(nhC,1)=mean(zh(ir{k}));
    ##           case 'min', zhC(nhC,1)=min(zh(ir{k}));
    ##           case 'max', zhC(nhC,1)=max(zh(ir{k}));
    ##           case 'quantile', zhC(nhC,1)=quantest(zh(ir{k}),P);
    ##         end        
    ##       case 'soft'
    ##         nsC=nsC+1;
    ##         csC(nsC,:)=ch(ir{k}(1),:);
    ##         switch method
    ##           case 'studentT', 
    ##             zave(nsC,1)=mean(zh(ir{k}));
    ##             zvar(nsC,1)=var(zh(ir{k}));
    ##             nobvs(nsC,1)=length(zh(ir{k}));
    ##           case 'hist',
    ##             [n,x]=histscaled(zh(ir{k}),nbins);
    ##             dx=diff(x(1:2));
    ##             limi1=[x(1)-dx x' x(end)+dx];
    ##             probdens1=[0 n' 0];
    ##             nl1=length(limi1);
    ##             [probdens1]=proba2probdens(softpdftype,nl1,limi1,probdens1);
    ##             if nsC==1
    ##               nlC=nl1;
    ##               limiC=limi1;
    ##               probdensC=probdens1;
    ##             else
    ##               [softpdftype,nlC,limiC,probdensC]=probacat(softpdftype,nlC,limiC,probdensC,...
    ##                 softpdftype,nl1,limi1,probdens1);
    ##             end
    ##         end
    ##     end
    ##   end
    ##   switch method
    ##     case 'studentT'
    ##       [softpdftype,nlC,limiC,probdensC]=probaStudentT(zave,zvar,nobvs);
    ##   end
    ## end

    ## if isempty(nl) 
    ##   return; 
    ## else
    ##   csC=[csC;cs];
    ##   [softpdftype,nlC,limiC,probdensC]=probacat(softpdftype,nlC,limiC,probdensC,...
    ##     softpdftype,nl,limi,probdens);
    ##   if ~isdupli([chC;csC]), return; end;
    ## end

    ## %%% Remove each set of duplicate soft data collocated with a hard datum
    ## nhC=size(chC,1);
    ## [iu,ir]=finddupli([chC;csC]);
    ## idxRemove=[];
    ## for k=1:length(ir)
    ##   if min(ir{k})<=nhC
    ##     idxRemove=[idxRemove;ir{k}(ir{k}>nhC)-nhC];
    ##   end
    ## end

    ## [csRemove,csC,nlRemove,limiRemove,probdensRemove,nlC,limiC,probdensC]=probasplit(csC,softpdftype,nlC,...
    ##   limiC,probdensC,idxRemove);

    ## [csC,nlC,limiC,probdensC,idxFailed] = probacombinedupli(csC,softpdftype,nlC,limiC,probdensC,'auto');

    ## if length(idxFailed)>0
    ##   str1='Some duplicated soft data points resulted in zero soft pdf when using multiplicative combination.';
    ##   str2='Use probacombinedupli on the soft data alone to see if there is something wrong with the soft data.';
    ##   str3='If using probacombinedupli does not reveal any failed multiplicative combination then the problem';
    ##   str4='might be with combining hard AND soft duplicates.  For these problematic duplicate data points we';
    ##   str5=' had to use an additive combination of soft data (see probacombinedupli for more explanations).';
    ##   disp(sprintf('warning in combinedupli: \n%s\n%s\n%s\n%s\n%s',str1,str2,str3,str4,str5));
    ## end

    nbins=20           # % number of bins used when method is equal to 'hist'

    if (isinstance(ch,array)!=true) or (isinstance(cs,array)!=true):
        return 'ch and cs must be numeric'
    [nh,d]=ch.shape
    if zh.shape[0]!=nh:
        return 'zh must have the same number of rows as ch'
    if (nh!=0) and (zh.shape[1]!=1):
        return 'zh must have one column'
    if (isinstance(method,Float)!=True) and (isinstance(method,str)!=True):
        return 'method must be numeric or character'
    if isinstance(method,Float):
      if len(method)!=1:
          return 'when method is numeric, it must have a length of 1'
      elif (0<method & method<1)!=True:
          return 'when method is numeric, it must be such that 0<method<1'
      else: 
        P=method 
        method='quantile'
    [ns]=cs.shape[0]
    if ns!=0 and nh!=0 and cs.shape[1]!=d:
        return 'cs must have the same number of columns as ch'
    if softpdftype==[]:
        softpdftype==2
    softpdftypeCheckArgs(softpdftype,nl,limi,probdens)
    if limi.shape[0]!=ns:
        return 'limi must have the same number of rows as cs'
    if probdens.shape[0]!=ns:
        return 'probdens must have the same number of rows as cs'
    if method=='ave' or method=='min' or method=='max' or method=='quantile':
        methodtype='hard'
    elif method=='studentT':
        methodtype='soft'
        if softpdftype!=2:
            return 'When method is equal to ''studentT'' then softpdftype must be equal to 2'
    elif method== 'hist':
        methodtype='soft'
        if softpdftype!=2: 
          return 'When method is equal to ''hist'' then softpdftype must be equal to 2'
    else:
        return 'Bad value for method'

    csC=array([])
    softpdftypeC=softpdftype
    nlC=array([])
    limiC=array([])
    probdensC=array([])

    ##%%%%%% Combine the hard data
    if nh==0 or isdupli(ch)==True:
        chC=ch
        zhC=zh
    else:
      [iu,ir]=finddupli(ch)
      chC=ch[iu,:]
      zhC=zh[iu]
      nhC=len(iu)
      nsC=0  
      for k in xrange(0,len(ir)):
          if methodtype=='hard':
              nhC=nhC+1
              chC[nhC,:]=ch[ir[k][0],:]
              if method=='ave':
                  zhC[nhC,0]=mean(zh[ir[k]])
              elif method=='min':
                  zhC[nhC,0]=min(zh[ir[k]])
              elif method=='max':
                  zhC[nhC,0]=max(zh[ir[k]])
              elif method=='quantile':
                  zhC[nhC,0]=quantest(zh[ir[k]],P)
          elif methodtype=='soft':
              nsC=nsC+1
              csC[nsC,:]=ch[ir[k][0],:]
              if method=='studentT':
                  zave[nsC,0]=mean(zh[ir[k]])
                  zvar[nsC,0]=var(zh[ir[k]])
                  nobvs[nsC,0]=len(zh[ir[k]])
              elif method=='hist':
                  [n,x]=histscaled(zh[ir[k]],nbins)
                  dx=diff(x[0:1])
                  limi1=concatenate((reshape(array([x[0]-dx]),(1,1)),reshape(x,(len(x),1)),reshape(array([x[-1]+dx]),(1,1))),0)
                  probdens1=concatenate((reshape(array([0]),(1,1)),reshape(n,(len(n),1)),reshape(array([0]),(1,1))),0)
                  nl1=len(limi1)
                  [probdens1]=proba2probdens(softpdftype,nl1,limi1,probdens1)
                  if nsC==1:
                      nlC=nl1
                      limiC=limi1
                      probdensC=probdens1
                  else:
                      [softpdftype,nlC,limiC,probdensC]=probacat(softpdftype,nlC,limiC,probdensC,softpdftype,nl1,limi1,probdens1)
 
      if method=='studentT':
          [softpdftype,nlC,limiC,probdensC]=probaStudentT(zave,zvar,nobvs)

    if len(nl)==0:
        return chC,zhC,csC,softpdftypeC,nlC,limiC,probdensC
    else:
        csC=concatenate((csC,cs))
        [softpdftype,nlC,limiC,probdensC]=probacat(softpdftype,nlC,limiC,probdensC,softpdftype,nl,limi,probdens)
        if isdupli(concatenate((chC,csC),0))!=True:
            return chC,zhC,csC,softpdftypeC,nlC,limiC,probdensC

    ##%%% Remove each set of duplicate soft data collocated with a hard datum
    nhC=chC.shape[0]
    [iu,ir]=finddupli(concatenate((chC,csC),0))
    idxRemove=[]
    for k in xrange(0,len(ir)):
        if min(ir[k])<=nhC:
            idxRemove.append(ir[k][ir[k]>nhC]-nhC)
    asarray(idxRemove)

    [csRemove,csC,nlRemove,limiRemove,probdensRemove,nlC,limiC,probdensC]=probasplit(csC,softpdftype,nlC,limiC,probdensC,idxRemove)

    [csC,nlC,limiC,probdensC,idxFailed] = probacombinedupli(csC,softpdftype,nlC,limiC,probdensC,'auto')

    if len(idxFailed)>0:
                   print 'warning in combinedupli:'
                   print 'Some duplicated soft data points resulted in zero soft pdf when using multiplicative combination.'
                   print 'Use probacombinedupli on the soft data alone to see if there is something wrong with the soft data.'
                   print 'If using probacombinedupli does not reveal any failed multiplicative combination then the problem'
                   print 'might be with combining hard AND soft duplicates.  For these problematic duplicate data points we'
                   print ' had to use an additive combination of soft data (see probacombinedupli for more explanations).'

    return chC,zhC,csC,softpdftypeC,nlC,limiC,probdensC


def coord2dist(c1,c2):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: debugged, not tested
    
    % coord2dist                - distance matrix from coordinates (Jan 1,2001)
    %
    % Computes the Euclidean distance matrix between two sets of coordinates.
    %
    % SYNTAX :
    %
    % [D]=coord2dist(c1,c2);
    %
    % INPUT :
    %
    % c1      n1 by d   matrix of coordinates for the locations in the first set.
    %                   A line corresponds to the vector of coordinates at a location,
    %                   so the number of columns is equal to the dimension of the space.
    %                   There is no restriction on the dimension of the space.
    % c2      n2 by d   matrix of coordinates for the locations in the second set, using
    %                   the same conventions as for c1. 
    %
    % OUTPUT :
    %
    % D       n1 by n2  Euclidean distance matrix between the coordinates in c1 and c2.
    %                   The number of lines and columns for D corresponds to the number
    %                   of lines in c1 and the number of lines in c2, respectively.
    %
    % NOTE :
    %
    % It is possible to specify additional n1 by 1 and n2 by 1 index vectors, taking
    % integer values from 1 to nv. The values in the index vectors specify which one
    % of the nv variable is known at each one of the corresponding coordinates. The c1
    % and c2 matrices of coordinates and the index vectors are then grouped together
    % using the MATLAB cell array notation, so that c1={c1,index1} and c2={c2,index2}.
    '''
        
    islist=isinstance(c1,list)
    if islist==True :
        c1=c1[0]
        c2=c2[0]

    if len(c1.shape)==1 :
        c1=reshape(c1,(1,len(c1)))
    if len(c2.shape)==1 :
        c2=reshape(c2,(1,len(c2)))
        
    n1=c1.shape[0]
    n2=c2.shape[0]
    unit1=ones((n1,1))
    unit2=ones((n2,1))
    D=BMEmatlab.kron(unit2,c1)-BMEmatlab.kron(c2,unit1)
    if D.shape[1]==1 :
        D=abs(D)
    else :
        D=sqrt(transpose(sum(transpose(D**2))))

    D=reshape(D,(n1,n2))

    return D



def coord2K(c1,c2,model,param,filtmodel=None):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: not competed: dependencies to complete
    
    % coord2K                   - covariance/variogram matrix from coordinates (Jan 1,2001)
    %
    % Compute the covariance or variogram matrix between two sets of coordinates,
    % based on the Euclidean distances between these sets of coordinates.
    %
    % SYNTAX :
    %
    % [K]=coord2K(c1,c2,model,param,filtmodel);
    %
    % INPUT :
    %
    % c1        n1 by d  matrix of coordinates for the locations in the first set.
    %                    A line corresponds to the vector of coordinates at a location,
    %                    so the number of columns is equal to the dimension of the space.
    %                    There is no restriction on the dimension of the space.
    % c2        n2 by d  matrix of coordinates for the locations in the second set, using
    %                    the same conventions as for c1. 
    % model     string   that contains the name of the variogram or covariance model that
    %                    is used for the computation.
    % param     k by 1   vector of values for the parameters of model, according to the
    %                    convention for the corresponding variogram or covariance model
    %                    (see the MODELS directory).
    % filtmodel scalar   optional value specifying if the stochastic model component is to
    %                    be included (filtmodel=1) or is to be filtered out (filtmodel=0).
    %                    See the krigingfilter.m function in the KRIGING directory.
    %                    
    %
    % OUTPUT :
    %
    % K         n1 by n2 covariance or variogram matrix between the coordinates in c1 and c2,
    %                    depending on the fact that model is a covariance or a variogram model.
    %                    The number of lines and columns for K corresponds to the number of
    %                    lines in c1 and the number of lines in c2, respectively.
    %
    % NOTE :
    %
    % For a detailed discussion about the coding of the model and param variables in
    % various situations (e.g., nested models, multivariate case, space/time case),
    % the reader is referred to the detailed description given for the kriging.m
    % function in the KRIGING directory.
    '''

    ## %%%%%% Build the general data structure for singular cases

    ## noindex=(iscell(c1)==0);
    ## nocell=(iscell(model)==0);

    ## if nocell==1 & noindex==1,          % case for 1 model and 1 variable
    ##   n1=size(c1,1);                    % compute the number of coordinates
    ##   n2=size(c2,1);
    ##   c1={c1,ones(n1,1)};               % create an index with values 1
    ##   c2={c2,ones(n2,1)};
    ##   nm=1;                             % set number of models equal to 1
    ##   nv=1;                             % set number of variables equal to 1
    ##   model={model};                    % create a cell array with model
    ##   np=length(param);                 % create a cell array with param
    ##   param={{param(1),param(2:np)}};
    ## end;

    ## if nocell==1 & noindex==0,          % case for 1 model and nv variables
    ##   n1=size(c1{1},1);                 % compute the number of coordinates
    ##   n2=size(c2{1},1);
    ##   nm=1;                             % set number of models equal to 1
    ##   nv=max([c1{2};c2{2}]);            % compute the number of variables 
    ##   model={model};                    % create a cell array with model
    ##   param={param};
    ## end;

    ## if nocell==0 & noindex==1,          % case for nm models and 1 variable
    ##   n1=size(c1,1);                    % compute the number of coordinates
    ##   n2=size(c2,1);
    ##   c1={c1,ones(n1,1)};               % create an index with values 1
    ##   c2={c2,ones(n2,1)};
    ##   nm=length(model);                 % compute the number of models 
    ##   nv=1;                             % set the number of variables to 1 
    ##   for i=1:nm,                       % create cell arrays with param vectors 
    ##     np=length(param{i});               
    ##     param{i}={param{i}(1),param{i}(2:np)};
    ##   end;
    ## end;

    ## if nocell==0 & noindex==0;          % case for nm models and nv variables
    ##   n1=size(c1{1},1);                 % compute the number of coordinates
    ##   n2=size(c2{1},1);
    ##   nm=length(model);                 % compute the number of models
    ##   nv=max([c1{2};c2{2}]);            % compute the number of variables 
    ## end;

    ## %%%%%% Initialize the parameters

    ## if nargin<5,
    ##   filtmodel=ones(nm,1);
    ## end;

    ## %%%%%% compute the distance matrix and initialize the size of K

    ## [isST,isSTsep,modelS,modelT]=isspacetime(model);

    ## if isSTsep,
    ##   [nparam,models]=nparammodels;
    ## end;

    ## if (n1==0)|(n2==0),
    ##   K=zeros(n1,n2);
    ##   return;
    ## else
    ##   if ~isST,
    ##     D=coord2dist(c1{1},c2{1});
    ##     K=zeros(size(D));
    ##   else
    ##     nd=size(c1{1},2)-1;
    ##     Ds=coord2dist(c1{1}(:,1:nd),c2{1}(:,1:nd));
    ##     Dt=coord2dist(c1{1}(:,nd+1),c2{1}(:,nd+1));
    ##     K=zeros(size(Ds));
    ##   end;
    ## end;

    ## %%%%%% Compute the covariance matrix or vector

    ## if n1==n2,                            % test if K is symmetric and set
    ##   issymm=prod(prod(double(c1{1}==c2{1})));    % issymm to 0 or 1 
    ## else
    ##   issymm=0;
    ## end;

    ## if issymm==1,                 % if K is symmetric 
    ##   for i=1:nv,                 % loop twice over the number of variables
    ##     indexi=find(c1{2}==i);    % and the number of models by selecting the 
    ##     for j=i:nv,               % appropriate subset of locations and model
    ##       indexj=find(c2{2}==j);  % and use the symmetry property of the matrix
    ##       for k=1:nm,
    ##         if filtmodel(k)==1,
    ## 	  C=param{k}{1};
    ##           if ~isST,
    ## 	    K(indexi,indexj)=K(indexi,indexj)+eval...
    ##                               ([model{k},'(D(indexi,indexj),[C(i,j),param{k}{2}])']);
    ##           end;
    ##           if (~isSTsep)&(isST),
    ## 	    K(indexi,indexj)=K(indexi,indexj)+eval...
    ##             ([model{k},'(Ds(indexi,indexj),Dt(indexi,indexj),[C(i,j),param{k}{2}])']);
    ##           end;
    ##           if isSTsep,
    ##             nps=nparam{strmatch(modelS{k},models,'exact')};
    ##             npt=nparam{strmatch(modelT{k},models,'exact')};
    ## 	    Ks=eval([modelS{k},'(Ds(indexi,indexj),[1,param{k}{2}(1:nps-1)])']);
    ## 	    Kt=eval([modelT{k},'(Dt(indexi,indexj),[1,param{k}{2}(nps:nps+npt-2)])']);
    ##             K(indexi,indexj)=K(indexi,indexj)+C(i,j)*Ks.*Kt;
    ##           end;
    ##         end;
    ##       end;
    ##       if i~=j,
    ##         K(indexj,indexi)=K(indexi,indexj)';
    ##       end;
    ##     end;
    ##   end;
    ## else                          % else K is not symmetric
    ##   for i=1:nv,                 % loop twice over the number of variables
    ##     indexi=find(c1{2}==i);    % and the number of models by selecting the
    ##     for j=1:nv,               % appropriate subset of locations and model
    ##       indexj=find(c2{2}==j);
    ##       for k=1:nm,
    ##         if filtmodel(k)==1,
    ## 	  C=param{k}{1};
    ##           if ~isST,
    ##             K(indexi,indexj)=K(indexi,indexj)+eval...
    ##                               ([model{k},'(D(indexi,indexj),[C(i,j),param{k}{2}])']);
    ##           end;
    ##           if (~isSTsep)&(isST),
    ## 	    K(indexi,indexj)=K(indexi,indexj)+eval...
    ##             ([model{k},'(Ds(indexi,indexj),Dt(indexi,indexj),[C(i,j),param{k}{2}])']);
    ##           end;
    ##           if isSTsep,
    ##             nps=nparam{strmatch(modelS{k},models,'exact')};
    ##             npt=nparam{strmatch(modelT{k},models,'exact')};
    ## 	    Ks=eval([modelS{k},'(Ds(indexi,indexj),[1,param{k}{2}(1:nps-1)])']);
    ## 	    Kt=eval([modelT{k},'(Dt(indexi,indexj),[1,param{k}{2}(nps:nps+npt-2)])']);
    ##             K(indexi,indexj)=K(indexi,indexj)+C(i,j)*Ks.*Kt;
    ##           end;
    ##         end;
    ##       end;
    ##     end;
    ##   end;
    ## end;

    ##%%%%%% Build the general data structure for singular cases
    
    noindex=(isinstance(c1,list)==False)
    nocell=(isinstance(model,list)==False)

    if nocell==True and noindex==True:          #% case for 1 model and 1 variable
        n1=c1.shape[0]                          #% compute the number of coordinates
        n2=c2.shape[0]
        c1=[c1,ones((n1,1))]                    #% create an index with values 1
        c2=[c2,ones((n2,1))]
        nm=1                                    #% set number of models equal to 1
        nv=1                                    #% set number of variables equal to 1
        model=[model]                           #% create a cell array with model
        np=len(param)                           #% create a cell array with param
        param=[[param[0],param[1:np]]]

    if nocell==True and noindex==False:         #% case for 1 model and nv variables
        n1=c1[0].shape[0]                       #% compute the number of coordinates
        n2=c2[0].shape[0]
        nm=1                                    #% set number of models equal to 1
        nv=max(concatenate((c1[1],c2[1])))      #% compute the number of variables
        model=[model]                           #% create a cell array with model
        param=[param]


    if nocell==False and noindex==True:          #% case for nm models and 1 variable
        n1=c1.shape[0]                           #% compute the number of coordinates
        n2=c2.shape[0]
        c1=[c1,ones((n1,1))]                     #% create an index with values 1
        c2=[c2,ones((n2,1))]
        nm=len(model)                            #% compute the number of models
        nv=1                                     #% set the number of variables to 1
        for i in xrange(0,nm):                   #% create cell arrays with param vectors
            np=len(param[i])
            param[i]=[param[i][0],param[i][1:np]]
 
    if nocell==False and noindex==False:         #% case for nm models and nv variables
        n1=c1[0].shape[0]                        #% compute the number of coordinates
        n2=c2[0].shape[0]
        nm=len(model)                            #% compute the number of models
        nv=max(concatenate((c1[1],c2[1])))       #% compute the number of variables

    ##%%%%%% Initialize the parameters

    if filtmodel==None:
        filtmodel=ones((nm,1))

    ##%%%%%% compute the distance matrix and initialize the size of K

    [isST,isSTsep,modelS,modelT]=isspacetime(model)

    if isSTsep:
      [nparam,models]=nparammodels

    if (n1==0) or (n2==0):
        K=zeros((n1,n2))
        return
    else:
      if isST==False:
          D=coord2dist(c1[0],c2[0])
          K=zeros(D.shape)
      else:
          nd=c1[0].shape[1]-1
          Ds=coord2dist(c1[0][:,0:nd],c2[0][:,0:nd])
          Dt=coord2dist(c1[0][:,nd+1],c2[0][:,nd+1])
          K=zeros(Ds.shape)

    ##%%%%%% Compute the covariance matrix or vector

    if n1==n2:                                      #% test if K is symmetric and set
        issymm=prod(prod(double(c1[0]==c2[0])))     #% issymm to 0 or 1 
    else:
      issymm=0

    if issymm==1:                      #% if K is symmetric
        for i in xrange(0,nv):         #% loop twice over the number of variables
            indexi=find(c1[1]==i)      #% and the number of models by selecting the
            for j in xrange(i,nv):     #% appropriate subset of locations and model
                indexj=find(c2[1]==j)  #% and use the symmetry property of the matrix
                for k in xrange(0,nm):
                    if filtmodel[k]==1:
                        C=param[k][0]
                        if isST!=True:
                            K[indexi,indexj]=K[indexi,indexj]+model[k](D[indexi,indexj],[C[i,j],param[k][1]])
                        if (isSTsep!=True) and (isST):
                            K[indexi,indexj]=K[indexi,indexj]+model[k](Ds[indexi,indexj],Dt[indexi,indexj],[C[i,j],param[k][1]])
                        if isSTsep:
                            nps=nparam[strmatch(modelS[k],models,'exact')]
                            npt=nparam[strmatch(modelT[k],models,'exact')]
                            Ks=modelS[k](Ds[indexi,indexj],[1,param[k][1][0:nps-1]])
                            Kt=modelT[k](Dt[indexi,indexj],[1,param[k][1][nps:nps+npt-2]])
                            K[indexi,indexj]=K[indexi,indexj]+C[i,j]*Ks*Kt
                if i!=j:
                    K[indexj,indexi]=transpose(K[indexi,indexj])
    else:                                #% else K is not symmetric
        for i in xrange(0,nv):           #% loop twice over the number of variables
            indexi=find(c1[1]==i)        #% and the number of models by selecting the
            for j in xrange(0,nv):       #% appropriate subset of locations and model
                indexj=find(c2[1]==j)
                for k in xrange(0,nm):
                    if filtmodel[k]==1:
                        C=param[k][0]
                        if isST!=True:
                            K[indexi,indexj]=K[indexi,indexj]+model[k](D[indexi,indexj],[C[i,j],param[k][1]])
                        if (isSTsep!=True) and (isST):
                            K[indexi,indexj]=K[indexi,indexj]+model[k](Ds[indexi,indexj],Dt[indexi,indexj],[C[i,j],param[k][1]])
                        if isSTsep:
                            nps=nparam[strmatch(modelS[k],models,'exact')]
                            npt=nparam[strmatch(modelT[k],models,'exact')]
                            Ks=modelS[k](Ds[indexi,indexj],[1,param[k][1][:nps-1]])
                            Kt=modelT[k](Dt[indexi,indexj],[1,param[k][1][nps:nps+npt-2]])
                            K[indexi,indexj]=K[indexi,indexj]+C[i,j]*Ks*Kt

    return K


def finddupli(p):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: Approved
    
    %  finddupli             - find the duplicate coordinates             
    %
    %  Find in a set of points those that are in duplicates (i.e. having the same coordinates)
    %
    %  SYNTAX :
    % 
    %  [iu,ir]=finddupli(p);
    % 
    %  INPUT :
    %
    %  p          n by d     matrix of the coordinates of n points in a space of dimension d
    % 
    %  OUTPUT :
    %
    %  iu        nu by 1    vector with the indices of the points that are not in duplicates
    %                       When there are no duplicates, then p[iu,:] is equal to p
    %                       otherwise nu<n and p[iu,:] is the subset of p that do not have duplicate
    %                       coordinates
    %  ir        1 by nr    list of the duplicate points.  
    %                       When there are no duplicates, then ir is empty
    %                       Otherwise ir[k] is a k=th cluster of duplicate points, so that
    %                       p(ir[k],:) is the k-th subset of p having all the same coordinates
    '''

    iu=transpose(arrayrange(0,p.shape[0]))
    if isdupli(p)!=True:
        ir=[]
    else:
        [n,d]=p.shape
        i = arrayrange(0,n)[:,NewAxis]
        [spi,dumb] = BMEmatlab.sortrows(concatenate((p,i),1))
        ps = spi[:,:d]
        iss = spi[:,d].astype(int)
        ind = ps[1:,0] == ps[0:-1,0]
        for k in xrange(1,d):
            ind = logical_and(ind,(ps[1:,k] == ps[:-1,k]))
        ind = concatenate((zeros(1,),ind,zeros(1,)))
        fs = find(logical_and((ind[:-1] == 0),(ind[1:] == 1)))
        fe = find(logical_and((ind[:-1] == 1), (ind[1:] == 0)))
        irv=array([])
        ir=[]
        for k in xrange(0, len(fs)):
            ir.append(iss[fs[k]:fe[k]+1])
            irv=concatenate((irv,ir[k]))
        irv = asarray(irv)
        [iu,dumb]=BMEmatlab.setdiff(iu,irv)
        if len(iu)==0:
            iu=[]

    return iu,ir


def findpairs(c1,c2):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: debugged, not tested
    
    % findpairs                 - identify identical coordinates (Jan 1,2001)
    %
    % Find pairs of coordinates that are identical for two diferent
    % matrices of coordinates.
    %
    % SYNTAX :
    %
    % [index]=findpairs(c1,c2);
    %
    % INPUT :
    %
    % c1        n1 by d  matrix of coordinates, where d is the dimension
    %                    of the space.
    % c2        n2 by d  matrix of coordinates.
    %
    % OUTPUT :
    %
    % index     n by 2  vector of indices for identical coordinates
    %                   in c1 and c2. The first column of index 
    %                   refers to the c1 matrix, whereas the second
    %                   column of index refers to the c2 matrix.
    %
    % NOTE :
    %
    % It is possible to specify additional n1 by 1 and n2 by 1 index vectors,
    % taking integer values from 1 to nv. The values in the index vectors
    % specify which of the nv variable is known at each one of the corresponding
    % coordinates. The c1 and c2 matrices of coordinates and the index vectors
    % are then grouped together using the MATLAB cell array notation, so that
    % c1={c1,index1} and c2={c2,index2}. Two points are then considered as
    % identical if both their coordinates and their indexes are identical.
    % The c1 and c2 matrices cannot contain duplicated coordinates. 
    '''

    if isinstance(c1,list)==True:
        c1=concatenate((c1[0],c1[1]),1)
        c2=concatenate((c2[0],c2[1]),1)

    n1=c1.shape[0]
    n2=c2.shape[0]
    index=zeros((min((n1,n2)),2))
    compt=-1

    if n1<n2:
        for i in xrange(0,n1):
            distance=coord2dist(c2,c1[i,:])
            isnulldist=find(ravel(distance)==0)
            if len(isnulldist)!=0:
                compt=compt+1
                npairs=len(isnulldist)
                index[compt,:]=[i,isnulldist[0]]
    else:
        for i in xrange(0,n2):
            distance=coord2dist(c1,c2[i,:])
            isnulldist=find(ravel(distance)==0)
            if len(isnulldist)!=0:
                compt=compt+1
                npairs=len(isnulldist)
                index[compt,:]=[isnulldist[0],i]
                
    if compt==-1:
        index=[]
    else:
        index=index[:compt,:]

    return index

def index2noindex(c,z):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: debugged, not tested
    
    % index2noindex             - split grouped coordinates and values for several variables (Jan 1,2001)
    %
    % Do the reciprocal of the conversion made by noindex2index.m with reversed
    % input and output variables.
    %
    % SYNTAX :
    %
    % [c,z]=index2noindex(c,z);
    %
    % INPUT :
    %
    % c    1 by 2    cell array, where the first cell is a n by d matrix of
    %                coordinates and the second cell is a n by 1 vector of
    %                index values ranging from 1 to nv. There is no restriction
    %                on the dimension of the space.
    % z    n by 1    vector of values at the coordinates specified in c.
    %
    % OUTPUT :
    %
    % c    1 by nv   cell array, where each cell contains the ni by d matrix
    %                of coordinates for each of the nv variables, with ni the
    %                number of coordinates for the ith variable and d the dimension
    %                of the space. 
    % z    1 by nv   cell array, where each cell contains the column vector of
    %                values for one of the nv variables at the coordinates
    %                specified in c.
    '''

    index=c[1]
    nv=max(index)
    ctemp=[]
    ztemp=[]
 
    for i in xrange(0,nv[0]+1):
        isvari=find(ravel(index)==i)
        ctemp.append(take(c[0],isvari))
        ztemp.append(take(z,isvari))
    c=ctemp
    z=ztemp

    return c,z


def invdist(ck,ch,zh,power,nhmax,dmax,options=None):
    """
    
    Translated by Dimitri D''Or - October 2004

    Status: debugged, not tested
    
    % invdist                   - prediction using an inverse distance weighting (Jan 1,2001)
    %
    % Estimate at a set of locations the values of a variable, based on
    % the knowledge of the hard data values at another set of locations.
    % The method uses an inverse distance weighting, so that the
    % estimated value at a location is a weighted linear combination of
    % the hard data values at surrounding locations. The weights are
    % positive and proportional to the inverse of a power function of the
    % Euclidian distances between the estimation locations and the locations
    % where hard data values are known.
    %
    % SYNTAX :
    %
    % [zk]=invdist(ck,ch,zh,power,nhmax,dmax,options);
    %
    % INPUT :
    %
    % ck       nk by d   matrix of coordinates for the estimation locations.
    %                    A line corresponds to the vector of coordinates at
    %                    an estimation location, so the number of columns in
    %                    ck corresponds to the dimension of the space. There
    %                    is no restriction on the dimension of the space.
    % ch       nh by d   matrix of coordinates for the hard data locations,
    %                    with the same convention as for ck.
    % zh       nh by 1   vector of values for the hard data at the coordinates
    %                    specified in ch.
    % power    scalar    value of the power to be used in the computation of the weights.
    %                    Special cases are power=Inf for the nearest neighbour estimate
    %                    and power=0 for an equal weighting. Values between 0 and Inf
    %                    correspond to intermediate situations.
    % nhmax    scalar    maximum number of hard data values that are considered for the
    %                    computations at each estimation location.
    % dmax     scalar    maximum distance between an estimation location and existing hard
    %                    data locations. All hard data locations separated by a distance
    %                    smaller than dmax from an estimation location will be included in
    %                    the estimation process for that location, whereas other hard data
    %                    locations are neglected.
    % options  scalar    optional parameter that can be used if default value is not
    %                    satisfactory (otherwise this vector can simply be omitted from the
    %                    input list of variables). options is equal to 1 or 0, depending if
    %                    the user wants or does not want to display the order number of the
    %                    location which is currently processed by the function. Default
    %                    value is 0.
    %
    % OUTPUT :
    %
    % zk       nk by 1   vector of estimated values at the estimation locations. Values coded
    %                    as NaN mean that no estimation has been performed at that location,
    %                    due to the lack of available data in the neighbourhood.
    """
    
    ##%%%%%% Initialize the parameters
    
    if options==None :
      options=array([0])

    nk=ck.shape[0]           #% nk is the number of estimation points
    nh=ch.shape[0]           #% nh is the number of hard data points
    zk=zeros((nk,),Float32)
    zk[:]=NAN

    ##%%%%%% Main loop starts here

    for i in xrange(0,nk) :
      ck0=ck[i,:]
      [chlocal,zhlocal,dh,sumnhlocal,index]=neighbours(ck0,ch,zh,nhmax,dmax)
      if sumnhlocal>0 :
          if power==Inf :
              findmind=find(dh==min(dh))
              zk[i]=zhlocal[findmind[0]]
          else :
              lam=(1./dh)**power
              lam=lam/sum(lam)
              zk[i]=transpose(lam)*zhlocal
    
      if options[0]==1 :
          print 'str(i) + '/' + str(nk)'

    zk=reshape(zk,(len(zk),1))

    return zk
    
def isdupli(p):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: Approved
    
    %  isdupli                - True if there are duplicate coordinates             
    %
    %  SYNTAX :
    % 
    %  [c]=isdupli(p);
    % 
    %  INPUT :
    %
    %  p          n by d     matrix of the coordinates of n points in a space of dimension d
    % 
    %  OUTPUT :
    %
    %  c         logical     True if there are duplicates, False otherwise
    '''

    pr=BMEmatlab.unique(p,'rows')
    c = p.shape[0] > pr[0].shape[0]

    return c


def iso2aniso(ciso,angle,ratio):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: debugged, not tested
    
    % iso2aniso                 - convert isotropic to anisotropic coordinates (Jan 1,2001)
    %
    % Do the transformation which is reciprocal to the transformation made by
    % the aniso2iso.m function, so that applying successively both transformation
    % has no effect on the coordinates. The function maps a set of coordinates
    % in an isotropic space into a set of coordinates in an anisotropic one. 
    %
    % SYNTAX :
    %
    % [c]=iso2aniso(ciso,angle,ratio); 
    %
    % INPUT :
    %
    % ciso    n by d   matrix of coordinates for the locations in the isotropic
    %                  space. A line corresponds to the vector of coordinates at a
    %                  location, so that the number of columns corresponds to the
    %                  dimension of the space. Only two dimensional or three dimensional
    %                  space coordinates can be processed by this function.
    % angle   1 by d-1 vector of angle values that characterize the anisotropy. 
    %                  In a two dimensional space, angle is the trigonometric angle
    %                  between the horizontal axis and the principal axis of the
    %                  ellipse. In a three dimensional space, spherical coordinates
    %                  are used, such that angle(1) is the horizontal trigonometric
    %                  angle and angle(2) is the vertical trigonometric angle for the
    %                  principal axis of the ellipsoid. All the angles are measured
    %                  counterclockwise in degrees and are between -90� and 90�.
    % ratio   1 by d-1 vector that characterize the ratio for the length of the axes
    %                  for the ellipse (in 2D) or ellipsoid (in 3D). In a two dimensional
    %                  space, ratio is the length of the principal axis of the ellipse
    %                  divided by the length of the secondary axis, so that ratio>1. In a
    %                  three dimensional space, ratio(1) is the length of the principal
    %                  axis of the ellipsoid divided by the length of the second axis, 
    %                  whereas ratio(2) is length of the principal axis of the ellipsoid
    %                  divided by the length of the third axis, so that ratio(1)>1 and
    %                  ratio(2)>1.
    %
    % OUTPUT :
    %
    % c       n by d   matrix having the same size as ciso, that gives the new coordinates
    %                  in the anisotropic space.
    %
    % NOTE :
    %
    % It is possible to specify an additional index vector, taking integer values from 1
    % to nv. The values in index specifies which of the nv variable is known at each one
    % of the corresponding coordinates. The ciso matrix of coordinates and the index vector
    % are then grouped together using the MATLAB cell array notation, so that ciso={ciso,index}.
    % This allows to perform the same coordinate transformation at once on a set of possibly
    % different variables. The output variable c is then also a cell array that contains
    % both the new matrix of coordinates and the index vector.
    '''

    ##%%%%%% Determine the dimension of the space and set ratio 

    isindex=isinstance(ciso,list)
    if isindex==True:
        d=ciso[0].shape[1]
    else:
        d=ciso.shape[1]
    angle=angle*2*pi/360
    ratio=1/ratio

    ##%%%%%% When d<2 or d>3, error

    if (d<2) or (d>3):
        return 'iso2aniso.m requires coordinates in a 2D or 3D space'

    ##%%%%%% Case for d=2

    if d==2:
        R=[ cos(angle),  sin(angle),
           -sin(angle),  cos(angle)]

        if isindex==True:
            c=ciso
            c[0][:,1]=c[0][:,1]*ratio
            c[0]=c[0]*R
        else:
            c=ciso
            c[:,1]=c[:,1]*ratio
            c=c*R

    ##%%%%%% Case for d=3

    if d==3:
        phi=angle[0]
        teta=angle[1]
        ratioy=ratio[0]
        ratioz=ratio[1]

        R1=[ cos(phi),   sin(phi),  0,
            -sin(phi),   cos(phi),  0,
             0,          0,         1]
        R2=[ cos(teta),  0,  sin(teta),
             0,          1,  0,
            -sin(teta),  0,  cos(teta)]
        R=R2*R1

        if isindex==1:
            c=ciso
            c[0][:,1]=c[0][:,1]*ratioy
            c[0][:,2]=c[0][:,2]*ratioz
            c[0]=c[0]*R
        else:
            c=ciso
            c[:,1]=c[:,1]*ratioy
            c[:,2]=c[:,2]*ratioz
            c=c*R

    return c


def isspacetime(model):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: debugged, not tested
    
    % isspacetime               - syntaxical analysis of the string 'model' (Jan 1,2001)
    %
    % Syntaxical detection of a space-time covariance/variogram
    % model inside the character string named model.
    %
    % A separable covariance model is coded as 'covmodelS/covmodelT',
    % where covmodelS refers to the spatial covariance function and
    % covmodelT refers to the temporal covariance function.
    % A non-separable covariance model is coded as 'covmodelST'.
    %
    % SYNTAX :
    %
    % [isST,isSTsep,modelS,modelT]=isspacetime(model);
    %
    % INPUT :
    %
    % model     string   that contains the name of the covariance/
    %                    variogram model.
    %
    % OUTPUT :
    %
    % isST      scalar   equal to 1 for a space-time model and equal
    %                    to 0 otherwise.
    % isSTsep   scalar   equal to 1 for a space-time separable model
    %                    and equal to 0 otherwise.
    % modelS    string   that contains the name of the spatial covariance
    %                    model if the space-time model is separable, or
    %                    is empty otherwise
    % modelT    string   that contains the name of the temporal covariance
    %                    model if the space-time model is separable, or
    %                    is empty otherwise.
    %
    % NOTE :
    %
    % If model is a cell array and the space-time model is separable,
    % modelS and modelT are cell arrays having the same dimension as model.
    '''

    modelS=[]
    modelT=[]
    onemodel=0

    if isinstance(model,list)!=True:
        model=[model]
        onemodel=1

    nc=len(model)
    isST=zeros((nc,1))
    isSTsep=zeros((nc,1))
    for i in xrange(0,nc):
        findslash=model[i].find('/')
        if len(findslash)>0:
            isST[i]=1
            isSTsep[i]=1
            modelS[i]=model[i][0:findslash]
            modelT[i]=model[i][findslash+1:len(model[i])]
        findST='ST' in model[i]
        if findST==True:
            isST[i]=1

    if isST[0]:
        if (alltrue(logical_not(isST))==True):
            return 'All models must be jointly either space, time or space-time'
    else:
        if (alltrue(isST)!=1):
            return 'All models must be jointly either space, time or space-time'
    if isSTsep[0]==0:
        if (alltrue(logical_not(isSTsep))!=1):
            return 'All space-time models must be jointly either separable or non-separable'
    else:
        if (alltrue(isSTsep)!=1):
            return 'All space-time models must be jointly either separable or non-separable'


    isST=isST[0]
    isSTsep=isSTsep[0]

    if onemodel==1:
        if isSTsep==1:
            modelS=modelS[0]
            modelT=modelT[0]

    return isST,isSTsep,modelS,modelT


def kernelsmoothing(ck,ch,zh,v,nhmax,dmax,options=None):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: debugged, not tested
    
    % kernelsmoothing           - prediction using a Gaussian kernel smoothing (Jan 1,2001)
    %
    % Estimate at a set of locations the values of a variable, based on
    % the knowledge of the hard data values at another set of locations.
    % The method used here is a Gaussian kernel smoothing. The estimated 
    % value at a location is a weighted linear combination of hard data
    % values at surrounding locations. The weights are positive and 
    % proportional to the values of a Gaussian distribution evaluated at
    % the Euclidean distances between the estimation locations and the
    % locations where hard data values are known.  
    %
    % SYNTAX :
    %
    % [zk]=kernelsmoothing(ck,ch,zh,v,nhmax,dmax,options); 
    %
    % INPUT :
    %
    % ck       nk by d   matrix of coordinates for the estimation locations.
    %                    A line corresponds to the vector of coordinates at
    %                    an estimation location, so the number of columns in
    %                    ck corresponds to the dimension of the space. There
    %                    is no restriction on the dimension of the space.
    % ch       nh by d   matrix of coordinates for the hard data locations,
    %                    with the same convention as for ck.
    % zh       nh by 1   vector of values for the hard data at the coordinates
    %                    specified in ch.
    % v        scalar    variance of the isotropic Gaussian kernel distribution.
    %                    A higher values for v provides a higher smoothing for
    %                    the zk estimates.
    % nhmax    scalar    maximum number of hard data values that are considered for the
    %                    computations at each estimation location.
    % dmax     scalar    maximum distance between an estimation location and existing hard
    %                    data locations. All hard data locations separated by a distance
    %                    smaller than dmax from an estimation location will be included in
    %                    the estimation process for that location, whereas other hard data
    %                    locations are neglected.
    % options  scalar    optional parameter that can be used if default value is not
    %                    satisfactory (otherwise this vector can simply be omitted from the
    %                    input list of variables). options is equal to 1 or 0, depending if
    %                    the user wants or does not want to display the order number of the
    %                    location which is currently processed by the function. 
    %
    % OUTPUT :
    %
    % zk       nk by 1   vector of estimated values at the estimation locations. Values coded
    %                    as NaN mean that no estimation has been performed at that location,
    %                    due to the lack of available data in the neighbourhood.
    '''
    
    ##%%%%%% Initialize the parameters

    if options==None:
      options=array([0])

    nk=ck.shape[0]           #% nk is the number of estimation points
    nh=ch.shape[0]           #% nh is the number of hard data
    zk=zeros((nk,),Float32)*NAN

    ##%%%%%% Main loop starts here

    for i in xrange(0,nk):
      ck0=ck[i,:]
      [chlocal,zhlocal,dh,sumnhlocal,index]=neighbours(ck0,ch,zh,nhmax,dmax)
      if sumnhlocal>0 :
        lam=gausspdf(dh,array([0,v]))
        lam=lam/sum(lam)
        zk[i]=matrixmultiply(transpose(lam),zhlocal)

      if options[0]==1 :
        print 'str(i) + '/' + str(nk)'

    zk=reshape(zk,(len(zk),1))

    return zk



def neighbours(c0,c,Z,nmax,dmax):
    """
    
    Translated by Dimitri D''Or - October 2004

    Status: debugged, not tested
    
    % neighbours                - radial neighbourhood selection (Jan 1,2001)
    %
    % Select a subset of coordinates and variables based on
    % their distances from the coordinate c0.
    %
    % SYNTAX :
    %
    % [csub,Zsub,dsub,nsub,index]=neighbours(c0,c,Z,nmax,dmax);
    %
    % INPUT :
    %
    % c0      1 by d   vector of coordinates, where d is the dimension
    %                  of the space
    % c       n by d   matrix of coordinates
    % Z       n by k   matrix of values, where each column corresponds to
    %                  the values of a same variable and each line corresponds
    %                  to the values of the diferent variables at the corresponding
    %                  c coordinates
    % nmax    scalar     maximum number of lines of Z that must be kept.
    % dmax    scalar     maximum Euclidian distance between c0 and c.
    %
    % OUTPUT :
    %
    % csub    m by d   matrix which is a subset of lines of the c matrix (m<=n)
    % Zsub    m by k   matrix which is a subset of lines of the Z matrix.
    % dsub    m by 1   vector of distances between csub and c0.
    % nsub    scalar   length of the dsub vector.
    % index   m by 1   vector giving the ordering of the lines in csub with
    %                  respect to the initial matrix c.
    %
    % NOTE :
    %
    % 1-In the case of space/time coordinates, dmax is a vector of length 3,
    % and the last column of c0 and c is the temporal coordinate. In this case,
    % dmax(1) is the maximum spatial distance between the coordinate in c and c0,
    % dmax(2) is the maximum temporal distance between coordinate in c and c0,
    % and dmax(3) refers to the space/time metric, such that :
    %
    % space/time distance=spatial distance+dmax(3)*temporal distance.
    % 
    % The space/time distance is used to select the nmax closest coordinates c
    % from the estimation coordinates c0.
    %
    % 2- It is possible to specify additional 1 by 1 and n by 1 index vectors,
    % taking integer values from 1 to nv. The values in the index vectors specify
    % which of the nv variable is known at each one of the corresponding coordinates.
    % The c0 and c matrices of coordinates and the index vectors are then grouped
    % together using the MATLAB cell array notation, so that c0={c0,index0} and
    % c={c,index}.
    """

    ##%%%%%% Test if c is non empty and is cell array

    isST=isinstance(dmax,list)

    if type(nmax)!=array:
        nmax=array([nmax])
    
    if isST :
        if len(dmax)==2 :
            return 'Error: dmax must have 1 element for spatial cases and three elements for space/time cases'
        dmax=array([dmax])

    if isinstance(c,list) :
        noelements=(c[0].shape[0]==0)
        noindex=0
    else :
        noelements=(c.shape[0]==0)
        noindex=1

    if noelements==1 :
        Zsub=array([],Float32)
        dsub=array([],Float32)
        nsub=0
        index=array([])                 # array of integers
        if noindex==1 :
            csub=array([],Float32)
        else :
            csub=[[array([],Float32),array([])]]
    else :

    ##%%%%%% two steps data selection inside the neighbourhood 

      if noindex==1 :                                              #% if there is no index
          n=c.shape[0]                                             #% create one with all values=1
          c0=[c0,array([1])]
          c=[c,zeros((n,1))]
      else :
          n=c[0].shape[0]

      if isST!=True :                                              #% for the spatial or temporal case,
          d=coord2dist(c[0],c0[0])                                 #% compute the distances in space or time
          d=reshape(d,(d.shape[0],))
          index=find(d<=dmax)                                      #% find distances<=dmax
      else :                                                       #% for the space-time case
          nd=c[0].shape[1]-1                                       #% compute the dimension of the space
          ds=coord2dist(take(c[0],xrange(0,nd),1),c0[0][0:nd])     #% compute the distances in space
          dt=coord2dist(take(c[0],(nd,),1),array([c0[0][nd]]))     #% compute the distances in time
          ds=reshape(ds,(ds.shape[0],))
          dt=reshape(dt,(dt.shape[0],))
          index=find((ds<=dmax[0]) and (dt<=dmax[1]))              #% find distances in space & time<=dmax

    ##%%%%%% checking for the number of data

      nv=max(take(c[1],index))                                     #% determine maximum value of this index
      nv=nv[0]+1                                                   #% To convert the array in an Int
      nsub=NAN*ones((nv,),Int)                                     #% initialize the value of nsub to 0
      cindex=take(c[1],index)
      cindex=reshape(cindex,(cindex.shape[0],))
      indexi=[]
      for i in xrange(0,nv) :                                      #% loop over the variables
          select=find(cindex==i)                                   #% select coordinates for variable i
          indexi.append(take(index,select))                        #% select corresponding values of index
          nsub[i]=len(indexi[i])                                   #% compute the number of data for each variable
          if nsub[i]>nmax[i] :                                     #% when number of data exceeds nmax(i)
              if isST!=True :
                  di=take(d,indexi[i])
              else :
                  di=take(ds,indexi[i])+dmax[2]*take(dt,indexi[i]) #% combine distance in space and time
              di=reshape(di,(len(di),))
              S=argsort(di)                                        #% sort the distance in increasing order
              indexi[i]=take(S,xrange(0,nmax[i]))                  #% select the nmax(i) closest data
              nsub[i]=nmax[i]
      index=array([])
      for i in xrange(0,nv) :                                      #% concatenate the indexi{} as the variable index
          index=concatenate((index,indexi[i]))                     #% in increasing order of the index
      if isST!=True :                                              #% for the spatial or temporal case,
          dsub=take(d,index)                                       #% extract the subset of distances in space or time
      else :                                                       #% for the space-time case,
          dsub=concatenate((take(ds,(index,)),take(dt,(index,))))  #% extract the subset of distances in space and time
      nsub=sum(nsub)                                               #% compute the length of the dsub vector
      if len(index)!=0 :
          Zsub=take(Z,index,0)                                     #% extract the subset of lines in V
          if noindex==1 :                                          #% if there was no index
              csub=take(c[0],index,0)                              #% extract the subset of coordinates
          else :                                                   #% else extract the subset of coordinates and index
              csub[0]=take(c[0],index,0)
              csub[1]=take(c[1],index,0)
      else :
          Zsub=array([],Float32)
          dsub=array([],Float32)
          nsub=0
          if noindex==1 :
              csub=array([],Float32)
          else :
              csub=[[array([],Float32),array([])]]

    return (csub,Zsub,dsub,nsub,index)


def noindex2index(c,z):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: debugged, not tested
    
    % noindex2index             - separated coordinates and variables toward grouped (Jan 1,2001)
    %
    % Convert a set of matrix of coordinates and vector of values
    % into a single cell array of coordinates and a single array
    % of values. The cell array of coordinates has as first cell 
    % the grouped set of coordinates, and as second cell an index
    % vector with values ranging from 1 the number of input vector
    % of values.
    %
    % SYNTAX :
    %
    % [c,z]=noindex2index(c,z);
    %
    % INPUT :
    %
    % c    1 by nv   cell array, where each cell contains the ni by d matrix
    %                of coordinates for each of the nv variables, with ni the
    %                number of coordinates for the ith variable and d the dimension
    %                of the space. There is no restriction on the dimension of
    %                the space.
    % z    1 by nv   cell array, where each cell contains the column vector of
    %                values for one of the nv variables at the coordinates
    %                specified in c.
    %
    % OUTPUT :
    %
    % c    1 by 2    cell array, where the first cell is the vertical stacking
    %                of the matrices of coordinates in the input c variable, and
    %                the second cell is a vector of index values ranging from 1
    %                to nv. The n1 first index values are equal to 1, the following
    %                n2 index values are equal to 2, etc.
    % z    n by 1    vector that corresponds to a vertical stacking of the nv vector
    %                of values.
    %
    % NOTE :
    %
    % For example, if c1, c2 and c3 are the coordinates for the vector
    % of values z1, z2 and z3 having n1, n2 and n3 elements, respectively,
    % then
    %
    % [c,z]=noindex2index({c1,c2,c3},{z1,z2,z3});
    %
    % yield the output variables c and z such that
    %
    % c{1}=[c1;c2;c3];
    % c{2}=[1*ones(n1,1);2*ones(n2,1);3*ones(n3,1)];
    % z=[z1;z2;z3];
    '''

    ctemp=[]
    c1=NAN*ones((1,c[0].shape[1]))
    c2=NAN*ones((1,1))
    ctemp.append(c1)
    ctemp.append(c2)
    ztemp=reshape(array([NAN]),(1,1))
    nv=len(c)
    for i in xrange(0,nv):
        nic=c[i].shape[0]
        niz=z[i].shape[0]
        if nic!=niz:
            return 'Number of lines must correspond for the elements in the c and z cell arrays'
        ctemp[0]=concatenate((ctemp[0],c[i]))
        ctemp[1]=concatenate((ctemp[1],(i*ones((len(z[i]),1)))))
        ztemp=concatenate((ztemp,z[i]))

    c[0]=ctemp[0][1:]
    c[1]=ctemp[1][1:]
    z=ztemp[1:]

    return c,z


def split(c,z,index):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: debugged, not tested
    
    % split                     - split coordinates in two sets (Jan 1,2001)
    %
    % Split the matrix of coordinates and the  associated matrix of
    % values into two sets according to the index vector that refers
    % to line numbers in c.
    %
    % SYNTAX :
    %
    % [c1,c2,z1,z2]=split(c,z,index);
    %
    % INPUT :
    %
    % c    n by d       matrix of coordinates, where d is the dimension
    %                   of the space.
    % z    n by k       matrix of values, where each line refers to a set
    %                   of values at the corresponding c coordinates.
    %
    % OUPUT :
    %
    % c1   ni by d      matrix of coordinates, where ni is the length of
    %                   the index vector.
    % c2   (n-ni) by d  matrix of coordinates.
    % z1   ni by k      matrix of values.
    % z2   (n-ni) by k  matrix of values.
    %
    % NOTE :
    %
    % It is possible to specify an additional index vector for c, taking
    % integer values from 1 to nv. The values in this index specifies which
    % one of the nv variable is known at each one of the corresponding
    % coordinates. The c matrix of coordinates and the index vector are then
    % grouped together using the MATLAB cell array notation, so that c={c,index}.
    % This allows to perform the same coordinate transformation at once on a set
    % of possibly different variables. The output matrices c1 and c2 are then
    % cell arrays too, with the index splitted accordingly.
    '''

    if isinstance(c,list)!=True:
        c1=take(c,index)
        z1=take(z,index)
        c2=c
        c2[index,:]=[]
        z2=z
        z2[index,:]=[]
    else:
        c1=[take(c[0],index),take(c[1],index)]
        z1=take(z,index)
        c2=c
        c2[0][index,:]=[]
        c2[1][index]=[]
        z2=z
        z2[index,:]=[]

    return c1,c2,z1,z2



def trapezint(x,y,xmin,xmax):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: debugged, not tested
    
    % trapezint                 - integration using trapezoidal formula (Jan 1,2001)
    %
    % Compute the integral of a discretized definition of a function of
    % a single variable using a trapezoidal integration formula.
    %
    % SYNTAX :
    %
    % [I]=trapezint(x,y,xmin,xmax);
    %
    % INPUT :
    %
    % x      n by 1   vector of values for the variable.
    % y      n by 1   vector for the values of the function at the x values.
    % xmin   scalar   lower bound for the integration.
    % xmax   scalar   upper bound for the integration.
    %
    % OUTPUT :
    %
    % I      scalar   integral of the function y=f(x) between xmin and xmax.
    %
    % NOTE :
    %
    % Not-a-Number (NaN) values for x and y are automatically stripped out
    % when computing the integral. As the functions is using the interp1.m
    % MATLAB function, trapezint.m requires x be monotonic.
    '''

    x=sort(x)
    index=argsort(x)
    y=take(y,index)
    
    index=find(logical_and((isnan(x)!=True),(isnan(y)!=True)))
    x=take(x,index)
    y=take(y,index)

    if type(xmin)!=array:
        xmin=array([xmin])
    if type(xmax)!=array:
        xmax=array([xmax])        
    
    n=len(x)
    if xmin>x[0]:
        finterp1d=interp1d(x,y,'linear')
        ymin=finterp1d(xmin)
        cond=find(x>xmin)
        x=take(x,cond)
        y=take(y,cond)
        x=concatenate((xmin,x))
        y=concatenate((ymin,y))
    
    n=len(x)
    if xmax<x[n-1]:
        finterp1d=interp1d(x,y,'linear')
        ymax=finterp1d(xmax)
        cond=find(x<xmax)
        x=take(x,cond)
        y=take(y,cond)
        x=concatenate((x,xmax))
        y=concatenate((y,ymax))
    
    n=len(x)
    y=0.5*y[:n-1]+0.5*y[1:n]
    I=sum(y*diff(x))

    return I
