from math import *
import numpy
from scipy import *
from scipy.interpolate import *
from scipy.sparse import find

from modellib import *
import BMEmatlab
import mvnlib


def cdf2pdf(z,cdf):
    '''
    
    Translated by Didrik Pinte - July 2004

    Status: debugged, not tested
    
    % cdf2pdf                   - compute the pdf from the cdf (Jan 1,2001)
    %
    % Compute the values of the probability distribution function
    % based on a discrete definition of the cumulative distribution
    % function.
    %
    % SYNTAX :
    %
    % [pdf]=cdf2pdf(z,cdf);
    %
    % INPUT :
    %
    % z      n by 1   vector of values.
    % cdf    n by 1   vector of the cumulative distribution function values
    %                 at the z values.
    %
    % OUPUT :
    %          
    % pdf    n by 1   vector of the probability distribution function values
    %                 at the z values.
    %
    % NOTE :
    %
    % As the differentiation of the cumulative distribution function is
    % obtained using a finite difference scheme, it is recommended to
    % have a finely discretized definition of this distribution. 
    '''
    index1=argsort(z)                ## %%% Keep track of the rank for the sorted values
    z = sort(index1)
    cdf=take(cdf,index1)             ## %%% Sort the pdf values accordingly

    forward=diff(cdf)/diff(z)        ## %%% Compute the pdf using finite differences
    n=len(z)
    pdf=zeros(n,1);
    pdf[1]=forward[1];
    pdf[n]=forward[n-1];
    pdf[2:n-1]=0.5*forward[2:n-1]+0.5*forward[1:n-2];

    index2 = argsort(index1);        ## %%% Resort the values in the original order
    index1 = sort(index1)
    pdf=pdf(index2);

    return pdf


def cdfest(z,w=None):
    '''
    
    Translated by Didrik Pinte - July 2004

    Status: debugged, not tested
    
    % cdfest                    - experimental cumulative distribution (Jan 1,2001)
    %
    % Compute an estimate of the univariate cumulative distribution
    % function from a set of values for a variable. The theoretical
    % definition of this distribution is F(z)=P[Z<=z]. As the
    % distribution is computed from a discrete set of values, a
    % continuity correction is done such that the estimate Fest(z)=
    % (Pest[Z<=z]+Pest[Z<z])/2, where Pest(.) are probabilities
    % estimated from the frequencies.
    %
    % SYNTAX :
    %
    % [zfile,cdfzfile]=cdfest(z,w);
    %
    % INPUT :
    %
    % z         n by 1   vector of values.
    % w         n by 1   optional column vector of weights for the z values.
    %                    These weights must be positive and must sum to one
    %                    (see, e.g., decluster.m above). If w is not specified,
    %                    all the weights are taken as equal.
    %
    % OUPUT :
    %
    % zfile     m by 1   vector of z values sorted in ascending order, where
    %                    the duplicate z values have been removed. If there
    %                    are no duplicate values, zfile is simply the sorted
    %                    z vector (m<=n).
    % cdfzfile  m by 1   vector of estimated values of the cumulative
    %                    distribution function computed at the zfile values.
    '''
    
    n=len(z)
    if w==None:
        w=ones((n),Float32)/n

    index= argsort(z)
    z = sort(z)
    w = take(w,index)

    i=0                       ##%%% delete the duplicate values
    while i<(n-1) :   
        if z[i+1]==z[i]:
            w[i+1]=w[i+1]+w[i]
            mask = [1 for x in xrange(0,n)]
            mask[i]=0;
            z = compress(mask,z)
            w = compress(mask,w)
            n = n-1
        else: 
            i=i+1
    
    temp = w.tolist()
    temp.insert(0,0.0)
    w = array((temp))
    cdfzfile = cumsum( (0.5*w[0:n]) + (0.5*w[1:n+1]))

    return(z,cdfzfile)


def coregfit(d,V,o,model,param0,options=None):
    '''
    
    Translated by Didrik Pinte - July 2004

    Status: not debugged, not tested
    
    % BUG : verifier le fonctionnement global de la fonction :  probleme de stockage des resultats (array, list, etc.)
    
    % coregfit                  - fitting of the coregionalization model (Jan 1,2001)
    %
    % Use an iterated least squares based algorithm for fitting a
    % multivariate model of variograms or covariance functions that
    % have been estimated using the vario.m or covario.m function.
    % This joint model is also called the Linear Model of
    % Coregionalization. The various modeled (cross)variograms or
    % (cross)covariance functions are defined as a sum of the same
    % few basic models having identical parameters for all variables
    % except for the sill parameters (or the equivalent of the sill
    % parameters for variogram models).
    %
    % SYNTAX :
    %
    % [param]=coregfit(d,V,o,model,param0,options);
    %
    % INPUT :
    %
    % d         nc by 1     vector giving the sorted values of the mean distance separating
    %                       the pairs of points that belong to the same distance class. 
    % V         nv by nv    symmetric array of cells that contains the variograms and cross
    %                       variograms estimates, or the covariance and cross covariance
    %                       estimates, for the distance classes specified in d.
    %                       Diagonal cells contain the variograms or covariance, whereas
    %                       off-diagonal cells contain the cross variograms or covariances.
    %                       If Z is a column vector (only one variable), then V is simply a
    %                       column vector having same size as d.
    % o         nc by 1     vector giving the number of pairs of points that belong to the 
    %                       corresponding distance classes.
    % model     string      that contains the name of the variogram or covariance model for
    %                       which parameters are sought (see the MODELS directory for possible
    %                       names of variograms and covariance functions). 
    % param0    1 by k      initial guess values for the parameters of model, according to the
    %                       convention for the corresponding variogram or covariance model. As
    %                       the coregfit.m function only seeks for estimates of the sill parameter
    %                       of model for the different variables, the first value in param0 is
    %                       arbitrary. The other parameters are kept unchanged by the function.
    % options   1 by 1 or 5 vector of optional parameters that can be used if default values are
    %                       not satisfactory (otherwise this vector can simply be omitted from the
    %                       input list of variables), where :
    %                       options(1) displays the estimated and fitted multivariate models if
    %                       the value is set to one (default value is 0),
    %                       options(2) is the relative tolerance for the fitting (default value
    %                       is 1e-4),
    %                       options(3) is the maximum number of iterations (default value is 1000),
    %                       options(4) is 0 for ordinary least squares and 1 for weighted least
    %                       squares, where the weights are proportional to the number of pairs of
    %                       points in each distance class (default value is 1),
    %                       options(5) is 0 for least squares based solely on the variograms
    %                       (covariance functions), and is 1 for least squares fitting based on
    %                       cross variograms (cross covariance functions) too (default value is 1).
    % OUTPUT :
    %
    % param     1 by 2      cell array that contains the complete set of parameters (both estimated
    %                       and unchanged) associated with model (see the detailed description
    %                       given for kriging.m about the coding of these parameters for the various
    %                       possible cases).
    %
    % NOTE :
    %
    % For a detailed discussion about the coding of the model and param0
    % variables for nested models, the reader is referred to the detailed
    % description given for the kriging.m function. The model and param
    % output variables are formatted in a way that the user can use them
    % as input model and param variables for the kriging.m and associated
    % functions.
    '''
    ##%%%%%% Check if there are no NaN

    if len(nonzero(isnan(d))):
        error('Some distance classes do not contain pairs of points')

    ##%%%%%% Initialize the parameters 

    if options == None:        #% setup the default options if not provided
        options = [0,1e-4,1000,1,1]
    else:
        noptions=len(options)
        if noptions==1 :
            options.append(1e-4)
            options.append(1000)
            options.append(1)
            options.append(1)
            
    if V.shape[1]>1:        #% test if there is one variable
        nv=V.shape[1]       #% nv is the number of variables
        onevariable=0
    else:
        nv=1
        #V={V};              #% create a cell with the content of vector V
        onevariable=1
    
    nd=len(d);           #% nd is the number of distance classes

    if len(model)==1:                   # % if there is one model, create a second
        np=len(param0)                      # % identical one, evaluate the theoretical
        g=model[0](d,[1,param0[1:np-1]])      # % value at distances d and set onemodel=1 
        g=[g,g]
        nm=2
        onemodel=1
    else:                                    #% else evaluate the theoretical values at
        nm=len(model)                        #% distance d and set onemodel=0
        g= []
        for i in xrange(0,nm-1):
            np=len(param0[i]);
            print model[i]
            print str(i) + " " + str(np)
            print param0[i]
            g.append(model[i](d,param0[i]))
        onemodel=0

    if options[3]==1:       #% if weighted least squares is required
        w=o/sum(o)          # % set the weights proportional to the
    else:                    #% mean number of pairs of points
        w=ones(nd,1)/nd      #% else ordinary least squares is required
                             #% and the weigths are set equal to 1/nd

    Cinit=zeros((nv-1,nv-1),Float32)
    for i in xrange(0,nv-1):             #% define a diagonal initialization matrix
        for j in xrange(0,nv-1):         #% of coefficients
            if i==j:
                Cinit[i,i]=mean(V[i,i])
    
                                 #% set the initial matrix of C equal to the
                                 #% initialization matrix divided by the
                                 #% number of variogram models
    C= []
    [C.append(Cinit/nm) for x in xrange(nm) ]

    ##%%%%%% Iterations over the different structures 

    iter=0
    test=1
    WSSold=inf
    while ((test==1)&(iter<options[2])):
        iter+= 1
        for i in xrange(0,nm-1):
            index=range(0,nm)

            #index[i]=0    # A eclaircir !!!
            print index
            gt=copy.copy(g)
            gt[i]=None
            gi=g[i]
            Ct=copy.copy(C)
            Ct[i]=None
            dVt=[]
            
            for j in xrange(0,nv-1):
                if (len(dVt)<=j+1):
                    dVt.append([])
                for k in xrange(j,nv-1):
                    dVt[j].append(V[j,k])
                    for l in xrange(0,nm-1):
                        print dVt[j][k]
                        print gt[l]
                        print Ct[l]
                        # Erreur car Ct[l] n'existe pas --> voir ligne 281 !!
                        dVt[j,k]=dVt[j][k]-gt[l]*Ct[l][j,k]
            for j in xrange(0,nv-1):
                for k in xrange(j,nv-1):
                    Ctemp[j,k]=transpose((w*gi))*dVt[j][k]
                    Ctemp[k,j]=Ctemp[j,k]
            [Q,S] = MLab.eig(Ctemp)
            S = MLab.diag(S)
            condi = S<=0
            S=S*(~condi)
            S=S+condi*0.001*amin(S(~condi))
            S=MLab.diag(S)
            Ctemp=Q*S*transpose(Q)
            C[i]=Ctemp/(transpose(w)*(g[:,i]**2))
 

        ##%%%%%% Compute the fitting criterion

        WSS=0
        if options[4]==1:
            for i in xrange(0,nm-1):
                for j in xrange(0,nv-1):
                    for k in xrange(j,nv-1):
                        WSS=WSS+sum((w*(V[j,k]-g[:,i]*C[i][j,k]))**2)
                        if onemodel==1:
                            WSS=WSS/2

        else:
            for i in xrange(0,nm-1):
                for j in xrange(0,nv-1):
                    WSS=WSS+sum((w*(V[j,j]-g[:,i]*C[i][j,j]))**2)

        if abs((WSSold-WSS)/WSSold)< options[1]:
            test=0

        WSSold=WSS

    if onemodel==1:
        C=C[1]+C[2]
        C=[C]
        nm=1

    ##%%%%%% Display the fitting results 

    if options[0]==1:
        ##test=(ishold==1);
        for i in xrange(0,nv-1):
            for j in xrange(i,nv-1):
                minVij=amin(V[i,j])
                maxVij=amax(V[i,j])
                vfit=zeros(nd);
                for k in xrange(0,nm-1):
                    vfit=vfit+C[k][i,j]*g[:,k]
                digit = str(nv)+str(nv)+str((i)*nv+j)
                subplot(digit)
                plot(d,V[i,j],'.')
                plot(d,vfit)
                #set(gca,'FontSize',6)
                axis([0,amax(d),amin([0,-1.1*sign(minVij)*minVij]),amax([01.1*sign(maxVij)*maxVij])])
                plot([0,amax(d)],[0,0],':')
                xlabel("Distance")
                ylabel("Variogram-Covariance")
                title("Couple " + str(i) + "-" + str(j))
        #if test==0:
            #hold off

    ##%%%%%% Build the paramfit array of cells

    if onemodel==1: 
        if onevariable==1:
            param=[C[1],param0[2:np]]
        else:
            param=[C[1],param0[2:np]]
    else:
        for i in xrange(0,nm-1):
            np=len(param0[i])
            param[i]=[C[i],param0[i][2:np]]

    return param


def covario(c,Z,cl,method,options=[0]):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: not completed: plotting commands
    
    % covario                   - multivariate covariance function estimation (Jan 1,2001)
    %
    % Estimate the covariance functions and cross covariance
    % functions for a set of variables which are known at the
    % same set of coordinates. 
    %
    % SYNTAX :
    %
    % [d,C,o]=covario(c,Z,cl,method,options);
    %
    % INPUT : 
    %
    % c         n by d       matrix of coordinates for the locations where the
    %                        values are known. A line corresponds to the vector
    %                        of coordinates at a location, so the number of columns
    %                        is equal to the dimension of the space. There is no
    %                        restriction on the dimension of the space.
    % Z         n by nv      matrix of values for the variables. Each line is
    %                        associated with the corresponding line vector of
    %                        coordinates in the c matrix, and each column corresponds
    %                        to a different variable.
    % cl        nc+1 by 1    vector giving the limits of the distance classes that are
    %                        used for estimating the covariances and cross covariances.
    %                        The distance classes are open on the left and closed on
    %                        the right. The lower limit for the first class is >=0.
    % method    string       that contains the name of the method used for computing
    %                        the distances between pairs of locations. method='kron'
    %                        uses a Kronecker product, whereas method='loop' uses a loop
    %                        over the locations. Using the Kronecker product is faster
    %                        for a small number of locations but may suffer from memory
    %                        size limitations depending on the memory available as its
    %                        requires the storage of a distance matrix. The loop method
    %                        may be used whatever the number of data locations and must be
    %                        used if an Out of Memory error message is generated. Both
    %                        methods yield exactly the same estimates.
    % options  1 by 1,3 or 4 vector of optional parameters that can be used if default
    %                        values are not satisfactory (otherwise this vector can simply
    %                        be omitted from the input list of variables), where :
    %                        options(1) displays the estimated covariances if the value is
    %                        set to one (default value is 0),
    %                        options(2) and options(3) are the minimum and maximum values
    %                        for the angles to be considered, using the same conventions as
    %                        for the pairsplot.m function. Angles can only be specified for
    %                        planar coordinates, i.e. when the number of columns in c is
    %                        equal to two,
    %                        options(4) is equal to 0 if the mean is null and equal to 1 if
    %                        the mean is constant but non null (default value is 1).
    %
    % OUTPUT :
    %
    % d        nc by 1       vector giving the sorted values of the mean distance separating
    %                        the pairs of points that belong to the same distance class. 
    % C        nv by nv      symmetric array of cells that contains the covariance function
    %                        and cross covariance function estimates for the distance classes
    %                        specified in d. Diagonal cells contain the n by 1 vector of
    %                        covariance estimates, whereas off-diagonal cells contain the nc
    %                        by 1 vector of cross covariance estimates. If Z is a column
    %                        vector (only one variable), then C is simply a column vector having
    %                        same size as d.
    % o        nc by 1       vector giving the number of pairs of points that belong to the 
    %                        corresponding distance classes.
    %
    % NOTE :
    %
    % The d, C and o output variables can be used without modification as input
    % for the coregfit.m function. When a distance class do not contain any pairs
    % of points, the function output a warning message. The d and C elements for
    % the corresponding distance class are thus coded as NaN\'s, whereas the
    % corresponding o element is equal to 0.
    '''

    ##%%%%%% Initialize the parameters

    if not method.isalpha():
      return 'method should be a char string'
    
    cl = sort(cl)
    if cl[0]<0:
      return 'Minimum class distance must be >=0'
    
    n = c.shape[0]
    nc = len(cl)
    minim = cl[0]
    maxim = cl[-1]
    try:
        nv = Z.shape[1]
    except IndexError:
        Z = Z[:,NewAxis]
        nv=Z.shape[1]

    C = zeros((nv,nv,nc-1),Float32)

    noptions = len(options)

    if noptions>=3:
        a = options[1]*2*pi/360
        b = options[2]*2*pi/360
        if c.shape[1]!=2:
            return 'Angle limits are specified only for planar coordinates'
        if (a==b)|(amin([a,b])< -pi/2)|(amax([a,b])> pi/2):
            return 'Angle limits must be different and between or equal to -90 and 90'

    ##%%%%%% Substract the means from the data if needed

    if noptions==4:
        options4=options[3]
    else :
        options4=1

    if options4==1:
        for i in xrange(0,nv):
            Z[:,i]=Z[:,i]-mean(Z[:,i])

    if method =='kron':   ##%%%%% Uses a Kronecker product for computing distances

        ##%%% Compute the distances

        unit=ones((n,1))
        dc=BMEmatlab.kron(unit,c)-BMEmatlab.kron(c,unit)
        if dc.shape[1]==1:
            dist=abs(dc)
        else:
            dist=sqrt(transpose(sum(transpose(dc**2))))

        ## %%% Compute the angles

        if noptions>=3:
            finddc1null=find(dc[:,0]==0)
            finddc1notnull=find(dc[:,0]!=0)
            ang=zeros((dc.shape[0],1),Float32)
            put(ang,finddc1null,(pi/2)*sign(take(dc[1],finddc1null)))
            put(ang,finddc1notnull,atan(take(dc[1],finddc1notnull)/take(dc[0],finddc1notnull)))

        ## %%% Select couples for appropriate distances and angles

        cond=(dist>max([0,minim])) & (dist<=maxim)
        if noptions>=3:
            conda=(ang>a)
            condb=(ang<=b)
            if a<b:
                cond=cond & (conda & condb)
            else:
                cond=cond & (conda | condb)
        dist=take(dist,cond)
        m=len(dist)
        if m==0:
            return 'No couples of values within the specified classes'

        ##%%% Loop over the number of variables and compute (cross)covariogram

        isclass=[[] for i in xrange(nc)]
        d=zeros((nc,1),Float32)*NaN
        o=zeros((nc,1),Float32)
        for k in xrange(nc):         
            isclass[k]=find((dist>cl[k])&(dist<=cl[k+1]))
            o[k]=len(isclass[k])/2
            if o[k]!=0:
                d[k]=sum(dist(isclass[k]))/(2*o[k])

        for i in xrange(0,nv):
            for j in xrange(i,nv):
                zi=Z[:,i]
                zj=Z[:,j]
                product=BMEmatlab.kron(unit,zi)*BMEmatlab.kron(zj,unit)
                product=take(product,cond)
                c=zeros((nc,1))*NaN
                for k in xrange(0,nc):
                    if o[k]!=0:
                        c[k]=sum(product(isclass[k]))/(2*o[k])
                C[i,j]=c
                if i!=j:
                    C[j,i]=c
 
    else:                      ##%%%%% Uses a loop over the data for computing distances
        d=zeros((nc,1))
        o=zeros((nc,1))
        for i in xrange(0,nv):
            for j in xrange(0,nv):
                C[i,j]=zeros((nc,1))

        for i in xrange(0,n):
            for j in xrange(0,n):
                dist=sqrt(sum((c[i,:]-c[j,:])**2))
                cond=(dist>max([0,minim]))&(dist<=maxim)
                if noptions==3:
                    dc=c[i,0:1]-c[j,0:1]
                    if dc[0]==0:
                        ang=(pi/2)*sign(dc[1])
                    else:
                        ang=atan(dc[1]/dc[0])
                    conda=(ang>a)
                    condb=(ang<=b)
                    if a<b:
                        cond=cond & (conda & condb)
                    else:
                        cond=cond & (conda | condb)
                if cond==1:
                    index=sum(dist>cl)
                    if (index>=1) & (index<=nc):
                        put(d,index,take(d,index)+dist)
                        put(o,index,take(o,index)+1)
                        for k in xrange(0,nv):
                            for l in xrange(k,nv):
                                put(C[k,l],index,take(C[k,l],index)+Z[i,k]*Z[j,l])

        for i in xrange(0,nc):
            if o[i]==0:
                d[i]=NaN
                for j in xrange(0,nv):
                    for k in xrange(j,nv):
                        C[j,k][i]=NaN
                        C[k,j][i]=NaN
            else:
                d[i]=d[i]/o[i]
                for j in xrange(0,nv):
                    for k in xrange(j,nv):
                        C[j,k][i]=C[j,k][i]/o[i]
                        C[k,j][i]=C[j,k][i]
                o[i]=o[i]/2

    ##%%%%%% display the computed covariogram if options[0]=1

    if options[0]==1:
        test=(hold==True)
        for i in xrange(0,nv):
            for j in xrange(i,nv):
                minCij=min(C[i,j])
                maxCij=max(C[i,j])
                subplot(nv,nv,(i-1)*nv+j)
                plot(d,C[i,j],'.')
                hold(True)
                text(d,C[i,j],str(o),fontsize='medium')
                set(gca(),fontsize=6)
                axis([0,max(d),min([0,-1.1*sign(minCij)*minCij]),max([0,1.1*sign(maxCij)*maxCij])])
                plot([0,max(d)],[0,0],':')
                xlabel('Distance',fontsize=8)
                ylabel('Covariogram',fontsize=8)
                title('Couple '+str(i)+'-'+str(j),fontsize=8)
        if test==False:
            hold(False)

    ##%%%%%% C is a vector if there is only one variable

    if nv==1:
         C=C[0]

    ###%%%%%% Check if there are no NaN

    if len(find(isnan(d))):
         print 'Warning : some distance classes do not contain pairs of points'

    return d,C,o


def crosscovario(c1,c2,z1,z2,cl,method,options=[0]):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: not completed: plotting commands
    
    % crosscovario              - single cross covariance function estimation (Jan 1,2001)
    %
    % Single covariance function or cross covariance function
    % estimation. To the opposite of the covario.m function,
    % crosscovario.m is is able to deal with the case where the
    % values of the two variables are given for two sets of
    % partially or totally different locations. 
    %
    % SYNTAX :
    %
    % [d,c,o]=crosscovario(c1,c2,z1,z2,cl,method,options);
    %
    % INPUT : 
    %
    % c1        n1 by d      matrix of coordinates for the locations where the
    %                        values of the first variable are known. A line
    %                        corresponds to the vector of coordinates at a location,
    %                        so the number of columns is equal to the dimension of
    %                        the space. There is no restriction on the dimension of
    %                        the space.
    % c2        n2 by d      matrix of coordinates for the locations where the values
    %                        of the second variable are known, using the same conventions
    %                        as for c1.
    % z1        n1 by 1      column vector of values for the first variable at the
    %                        coordinates specified in c1.
    % z2        n2 by 1      column vector of values for the second variable at the
    %                        coordinates specified in c2.
    % cl        nc+1 by 1    vector giving the limits of the distance classes that are
    %                        used for estimating the cross covariance. The distance classes
    %                        are open on the left and closed on the right. The lower limit
    %                        for the first class is >=0.
    % method    string       that contains the name of the method used for computing
    %                        the distances between pairs of locations. method='kron'
    %                        uses a Kronecker product, whereas method='loop' uses a loop
    %                        over the locations. Using the Kronecker product is faster
    %                        for a small number of locations but may suffer from memory
    %                        size limitations depending on the memory available as its
    %                        requires the storage of a distance matrix. The loop method
    %                        may be used whatever the number of data locations and must be
    %                        used if an Out of Memory error message is generated. Both
    %                        methods yield exactly the same estimates.
    % options  1 by 1,3 or 4 vector of optional parameters that can be used if default
    %                        values are not satisfactory (otherwise this vector can simply
    %                        be omitted from the input list of variables), where :
    %                        options(1) displays the estimated cross covariance if the value
    %                        is set to one (default value is 0),
    %                        options(2) and options(3) are the minimum and maximum values
    %                        for the angles to be considered, using the same conventions as
    %                        for the pairsplot.m function. Angles can only be specified for
    %                        planar coordinates, i.e. when the number of columns in c1 and
    %                        c2 is equal to two,
    %                        options(4) is equal to 0 if the mean is null and equal to 1 if
    %                        the mean is constant but non null (default value is 1).
    %
    % OUTPUT :
    %
    % d        nc by 1       vector giving the sorted values of the mean distance separating
    %                        the pairs of points that belong to the same distance class. 
    % c        nc by 1       vector of estimated covariance or cross covariance function
    %                        values (if z1 and z2 are identical, the function computes the
    %                        covariance function).
    % o        nc by 1       vector giving the number of pairs of points that belong to the 
    %                        corresponding distance classes.
    %
    % NOTE :
    %
    % Note that using crosscovario(c1,c2,z1,z2,cl,method) yields the
    % same result as using crosscovario(c2,c1,z2,z1,cl,method) both for
    % method='kron'and method='loop'.
    '''

    ##%%%%%% Initialize the parameters
    
    if not method.isalpha():
        return 'method should be a char string'
    
    cl=sort(cl)
    if cl[0]<0:
        return 'Minimum class distance must be >=0'
    
    n1=c1.shape[0]
    n2=c2.shape[0]
    nc=len(cl)-1
    minim=cl[0]
    maxim=cl[-1]
    
    noptions=len(options)
    
    if noptions>=3:
        a=options[1]*2*pi/360
        b=options[2]*2*pi/360
        if (c1.shape[1]!=2 & c2.shape[1]!=2):
            return 'Angle limits are specified only for planar coordinates'
        if (a==b) | (min([a,b])<-pi/2) | (max([a,b])>pi/2):
            return 'Angle limits must be different and between or equal to -90 and 90'
    
    ##%%%%%% Substract the means
    
    if noptions==4:
        options4=options[4]
    else:
        options4=1
    
    if options4==1:
        z1=z1-mean(z1)
        z2=z2-mean(z2)
    
    if method.find('kron')!=-1:   ##%%%%% Uses a Kronecker product for computing distances
        
        ##%%% Compute the distances
        
        unit1=ones((n1,1))
        unit2=ones((n2,1))
        dc=BMEmatlab.kron(unit1,c2)-BMEmatlab.kron(c1,unit2)
        if dc.shape[1]==1:
            dist=abs(dc)
        else:
            dist=sqrt(transpose(sum(transpose(dc**2))))
        
        ##%%% Compute the angles
        
        if noptions>=3:
            finddc1null=find(dc[:,0]==0)
            finddc1notnull=find(dc[:,0]!=0)
            ang=zeros((dc.shape[0],1),Float32)
            ang[finddc1null]=(pi/2)*sign(take(dc[:,1],finddc1null))
            ang[finddc1notnull]=atan(take(dc[:,1],finddc1notnull)/take(dc[:,0],finddc1notnull))
        
        ##%%% Select couples for appropriate distances and angles
        
        cond=(dist>max([0,minim])) & (dist<=maxim)
        if noptions>=3:
            conda=(ang>a)
            condb=(ang<=b)
            if a<b:
                cond=cond & (conda & condb)
            else:
                cond=cond & (conda | condb)
        dist=take(dist,cond)
        m=len(dist)
        if m==0:
            return 'No couples of values within the specified classes'
        
        ##%%% Compute the cross-covariogram
        
        isclass=[[] for i in xrange(0,nc)]
        d=zeros((nc,1),Float32)*NaN
        o=zeros((nc,1),Float32)
        for k in xrange(0,nc):
            isclass[k]=find((dist>cl[k]) & (dist<=cl[k+1]))
            o[k]=len(isclass[k])/2
            if o[k]!=0:
                d[k]=sum(dist(isclass[k]))/(2*o[k])
        
        product=BMEmatlab.kron(unit1,z2)*BMEmatlab.kron(z1,unit2)
        product=take(product,cond)
        c=zeros((nc,1),Float32)*NaN
        for k in xrange(0,nc):
            if o[k]!=0:
                c[k]=sum(product(isclass[k]))/(2*o[k])

    else:                      ##%%%%% Uses a loop over the data for computing distances
        d=zeros((nc,1),Float32)
        o=zeros((nc,1),Float32)
        c=zeros((nc,1),Float32)

        for i in xrange(0,n1):
            for j in xrange(0,n2):
                dist=sqrt(sum((c1[i,:]-c2[j,:])**2))
                cond=(dist>max([0,minim])) & (dist<=maxim)
                if noptions==3:
                    dc=c1[i,0:1]-c2[j,0:1]
                    if dc[0]==0:
                        ang=(pi/2)*sign(dc[1])
                    else:
                        ang=atan(dc[1]/dc[0])
                    conda=(ang>a)
                    condb=(ang<=b)
                    if a<b:
                        cond=cond & (conda & condb)
                    else:
                        cond=cond & (conda | condb)
                if cond==1:
                    index=sum(dist>cl)
                    if (index>=1) & (index<=nc):
                        put(d,index,take(d,index)+dist)
                        put(o,index,take(o,index)+1)
                        put(c,index,take(c,index)+z1[i]*z2[j])
        
        for i in xrange(0,nc):
            if o[i]==0:
                d[i]=NaN
                c[i]=NaN
            else:
                d[i]=d[i]/o[i]
                c[i]=c[i]/o[i]
                o[i]=o[i]/2

    ##%%%%%% display the computed cross-covariogram if options[0]=1

    if options[0]==1:
        test=(hold==True)
        minc=min(take(c,not isnan(c)))
        maxc=max(take(c,not isnan(c)))
        maxd=max(take(d,not isnan(d)))
        plot(d,c,'.')
        hold(True)
        text(d,c,str(o),fontsize=8)
        set(gca(),fontsize=6)
        axis([0,maxd,min([0,-1.1*sign(minc)*minc]),max([0,1.1*sign(maxc)*maxc])])
        plot([0,maxd],[0,0],':')
        xlabel('Distance',fontsize=8)
        ylabel('Cross-covariogram',fontsize=8)
        if test==0:
            hold(False)
    
    ##%%%%%% Check if there are no NaN

    if len(find(isnan(d))):
        print 'Warning : some distance classes do not contain pairs of points'
            
    return d,c,o


def crosscovarioST(c1,c2,z1,z2,cls,clt,options=[0]):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: not completed: plotting commands
    
    % crosscovarioST            - space/time cross covariance estimation (Jan 1,2001)
    %
    % Single space/time covariance or cross covariance estimation.
    % The function can be used when the values of the two variables
    % are given for two sets of partially or totally different
    % space/time locations. 
    %
    % SYNTAX :
    %
    % [ds,dt,c,o]=crosscovarioST(c1,c2,z1,z2,cls,clt,options);
    %
    % INPUT : 
    %
    % c1       n1 by d+1     matrix of space/time coordinates for the locations
    %                        where the values of the first variable are known. A
    %                        line corresponds to the vector of space/time coordinates
    %                        at a location, so the number of columns is equal to the
    %                        dimension of the space plus one, where the last column
    %                        refers to the temporal coordinate. There is no restriction
    %                        on the dimension of the space.
    % c2       n2 by d+1     matrix of space/time coordinates for the locations where
    %                        the values of the second variable are known, using the same
    %                        conventions as for c1.
    % z1       n1 by 1       column vector of values for the first variable at the
    %                        coordinates specified in c1.
    % z2       n2 by 1       column vector of values for the second variable at the
    %                        coordinates specified in c2.
    % cls      ncs+1 by 1    vector giving the limits of the spatial distance classes that
    %                        are used for estimating the covariance or cross covariance. The
    %                        distance classes are open on the left and closed on the right.
    %                        The lower limit for the first class is >=0.
    % clt      nct+1 by 1    vector giving the limits of the temporal "distance" classes
    %                        that are used for estimating the covariance or cross covariance.
    %                        As for cls, the classes are open on the left and closed on the
    %                        right. The lower limit for the first class is >=0.
    % options  1 by 1,3 or 4 vector of optional parameters that can be used if default
    %                        values are not satisfactory (otherwise this vector can simply
    %                        be omitted from the input list of variables), where :
    %                        options(1) displays the estimated cross covariance if the value
    %                        is set to one (default value is 0),
    %                        options(2) and options(3) are the minimum and maximum values
    %                        for the angles to be considered, using the same conventions as
    %                        for the pairsplot.m function. Angles can only be specified for
    %                        spatial planar coordinates, i.e. when the number of columns in
    %                        c1 and c2 is equal to three,
    %                        options(4) is equal to 0 if the mean is null and equal to 1 if
    %                        the mean is constant but non null (default value is 1).
    %
    % OUTPUT :
    %
    % ds       ncs by nct    matrix giving the sorted values of the mean spatial distance
    %                        separating the pairs of points that belong to the same space/time
    %                        distance class. Each line of ds correspond to the same spatial
    %                        class, where each column of ds correspond to the same temporal class.
    % dt       ncs by nct    matrix giving the sorted values of the mean temporal distance
    %                        separating the pairs of points that belong to the same space/time
    %                        distance class, with same conventions as for ds.
    % c        ncs by nct    matrix of estimated space/time covariance or cross covariance
    %                        values (if z1 and z2 are identical, the function computes the
    %                        space/time covariance), with same dimensions as ds and dt.
    % o        ncs by nct    matrix giving the number of pairs of points that belong to the
    %                        corresponding space/time distance classes, with same dimensions
    %                        as ds and dt.
    %
    % NOTE :
    %
    % Due to the nature of the problem itself, it is generally expected that
    % a large number of locations are involved in the computation of the
    % space/time covariance or cross covariance. For that reason and to the
    % opposite of the crosscovario.m function, only the loop method has been
    % implemented in the function.
    '''

    ## %%%%%% Initialize the parameters

    ## cls=sort(cls);
    ## clt=sort(clt);
    ## if (cls(1)<0)|(clt(1)<0),
    ##   error('Minimum space/time class distances must be >=0');
    ## end;
    ## n1=size(c1,1);
    ## n2=size(c2,1);
    ## ncs=length(cls)-1;
    ## nct=length(clt)-1;
    ## minims=cls(1);
    ## maxims=cls(ncs+1);
    ## minimt=clt(1);
    ## maximt=clt(nct+1);
    ## dim=size(c1,2)-1;

    ## if nargin==6,
    ##   options(1)=0;
    ##   noptions=1;
    ## else
    ##   noptions=length(options);
    ## end;

    ## if noptions==3,
    ##   a=options(2)*2*pi/360;
    ##   b=options(3)*2*pi/360;
    ##   if dim~=2,
    ##     error('Angle limits are specified only for planar coordinates');
    ##   end;
    ##   if (a==b)|(min([a,b])<-pi/2)|(max([a,b])>pi/2),
    ##     error('Angle limits must be different and between or equal to -90 and 90');
    ##   end;
    ## end;

    ## %%%%%% Substract the means

    ## if noptions==4,
    ##   options4=options(4);
    ## else 
    ##   options4=1;
    ## end;

    ## if options4==1,
    ##   z1=z1-mean(z1);
    ##   z2=z2-mean(z2);
    ## end;

    ## %%%%% Uses a loop over the data for computing distances

    ## ds=zeros(ncs,nct);
    ## dt=zeros(ncs,nct);
    ## o=zeros(ncs,nct);
    ## c=zeros(ncs,nct);

    ## for i=1:n1,
    ##   for j=1:n2,
    ##     dists=sqrt(sum((c1(i,1:dim)-c2(j,1:dim)).^2));
    ##     distt=abs(c1(i,dim+1)-c2(j,dim+1));
    ##     conds=(dists>max([0 minims]))&(dists<=maxims);
    ##     condt=(distt>max([0 minims]))&(distt<=maxims);
    ##     cond=(conds & condt);
    ##     if noptions==3,
    ##       dc=c1(i,1:2)-c2(j,1:2);
    ##       if dc(1)==0,
    ##         ang=(pi/2)*sign(dc(2));
    ##       else
    ##         ang=atan(dc(2)/dc(1));
    ##       end;
    ##       conda=(ang>a);
    ##       condb=(ang<=b);
    ##       if a<b,
    ##         cond=cond & (conda & condb);
    ##       else
    ##         cond=cond & (conda | condb);
    ##       end;
    ##     end;
    ##     if cond==1,
    ##       indexs=sum(dists>cls);
    ##       indext=sum(distt>clt);
    ##       if ((indexs>=1) & (indexs<=ncs) & (indext>=1) & (indext<=nct)),
    ##         ds(indexs,indext)=ds(indexs,indext)+dists;
    ##         dt(indexs,indext)=dt(indexs,indext)+distt;
    ##         o(indexs,indext)=o(indexs,indext)+1;
    ##         c(indexs,indext)=c(indexs,indext)+z1(i)*z2(j);
    ##       end;
    ##     end;
    ##   end;
    ## end;

    ## for i=1:ncs,
    ##   for j=1:nct,
    ##     if o(i,j)==0,
    ##       ds(i,j)=NaN;
    ##       dt(i,j)=NaN;
    ##       c(i,j)=NaN;
    ##     else
    ##       ds(i,j)=ds(i,j)/o(i,j);
    ##       dt(i,j)=dt(i,j)/o(i,j);
    ##       c(i,j)=c(i,j)/(2*o(i,j));
    ##       o(i,j)=o(i,j)/2;
    ##     end;
    ##   end;
    ## end;

    ## %%%%%% display the computed cross-covariance if options(1)=1

    ## if options(1)==1,
    ##   test=(ishold==1);
    ##   minc=min(c(:));
    ##   maxc=max(c(:));
    ##   posx=(cls(1:ncs)+cls(2:ncs+1))/2;
    ##   posy=(clt(1:nct)+clt(2:nct+1))/2;
    ##   maxx=max(ds(:));
    ##   maxy=max(dt(:));

    ##   subplot(1,2,1);
    ##     obj1=contour(posx,posy,c');
    ##     clabel(obj1,'FontSize',8);
    ##     set(gca,'FontSize',8);
    ##     axis([0 maxx 0 maxy]);
    ##     xlabel('Space distance','FontSize',8);
    ##     ylabel('Time distance','FontSize',8);
    ##     axis('square');
    ##   subplot(1,2,2);
    ##     mesh(posx,posy,c');
    ##     set(gca,'FontSize',8);
    ##     axis([0 maxx 0 maxy min([0;-1.1*sign(minc)*minc]) max([0;1.1*sign(maxc)*maxc])]);
    ##     xlabel('Space distance','FontSize',8);
    ##     ylabel('Time distance','FontSize',8);
    ##     zlabel('Cross-covariance','FontSize',8);
    ##     axis('square');

    ##   if test==0,
    ##     hold off;
    ##   end;
    ## end;

    ## %%%%%% Check if there are no NaN

    ## if length(find(isnan([ds;dt])))~=0,
    ##   disp('Warning : some space/time classes do not contain pairs of points');
    ## end;

    

    ##%%%%%% Initialize the parameters

    cls=sort(cls)
    clt=sort(clt)
    if (cls[0]<0) | (clt[0]<0):
        return 'Minimum space/time class distances must be >=0'
    n1=c1.shape[0]
    n2=c2.shape[0]
    ncs=len(cls)-1
    nct=len(clt)-1
    minims=cls[0]
    maxims=cls[-1]
    minimt=clt[0]
    maximt=clt[-1]
    dim=c1.shape[1]-1

    noptions=len(options)

    if noptions==3:
        a=options[1]*2*pi/360
        b=options[2]*2*pi/360
        if dim!=2:
            return 'Angle limits are specified only for planar coordinates'
        if (a==b) | (min([a,b])<-pi/2) | (max([a,b])>pi/2):
            return 'Angle limits must be different and between or equal to -90 and 90'

    ##%%%%%% Substract the means

    if noptions==4:
        options4=options[3]
    else:
        options4=1

    if options4==1:
        z1=z1-mean(z1)
        z2=z2-mean(z2)

    ##%%%%% Uses a loop over the data for computing distances

    ds=zeros((ncs,nct),Float32)
    dt=zeros((ncs,nct),Float32)
    o=zeros((ncs,nct),Float32)
    c=zeros((ncs,nct),Float32)

    for i in xrange(0,n1):
        for j in xrange(0,n2):
            dists=sqrt(sum((c1[i,:dim]-c2[j,:dim])**2))
            distt=abs(c1[i,dim+1]-c2[j,dim+1])
            conds=(dists>max([0,minims])) & (dists<=maxims)
            condt=(distt>max([0,minims])) & (distt<=maxims)
            cond=(conds & condt)
            if noptions==3:
                dc=c1[i,:1]-c2[j,:1]
                if dc[0]==0:
                    ang=(pi/2)*sign(dc[1])
                else:
                    ang=atan(dc[1]/dc[0])
                conda=(ang>a)
                condb=(ang<=b)
                if a<b:
                    cond=cond & (conda & condb)
                else:
                    cond=cond & (conda | condb)
            if cond==1:
                indexs=sum(dists>cls)
                indext=sum(distt>clt)
                if ((indexs>=1) & (indexs<=ncs) & (indext>=1) & (indext<=nct)):
                    put(ds,[indexs,indext],take(ds,[indexs,indext])+dists)
                    put(dt,[indexs,indext],take(dt,[indexs,indext])+distt)
                    put(o,[indexs,indext],take(o,[indexs,indext])+1)
                    put(c,[indexs,indext],take(c,[indexs,indext])+z1[i]*z2[j])

    for i in xrange(0,ncs):
        for j in xrange(0,nct):
            if o[i,j]==0:
                ds[i,j]=NaN
                dt[i,j]=NaN
                c[i,j]=NaN
            else:
                ds[i,j]=ds[i,j]/o[i,j]
                dt[i,j]=d[i,j]/o[i,j]
                c[i,j]=c[i,j]/(2*o[i,j])
                o[i,j]=o[i,j]/2

    ##%%%%%% display the computed cross-covariance if options(1)=1

    if options[0]==1:
        test=(hold==True)
        minc=min(ravel(c))
        maxc=max(ravel(c))
        posx=(cls[:ncs]+cls[1:ncs+1])/2
        posy=(clt[:nct]+clt[1:nct+1])/2
        maxx=max(ravel(ds))
        maxy=max(ravel(dt))

        subplot(1,2,1)
        obj1=contour(transpose(c),posx,posy)    ## REQUIRES THE CONTOUR MODULE
        clabel(obj1,fontsize=8)
        set(gca(),fontsize=8)
        axis([0,maxx,0,maxy])
        xlabel('Space distance',fontsize=8)
        ylabel('Time distance',fontsize=8)
        axis('square')
        subplot(1,2,2)
        mesh(posx,posy,transpose(c))           ## MESH : FIND EQUIVALENT FUNCTION (in PyGist ?)
        set(gca(),fontsize=8)
        axis([0,maxx,0,maxy,min([0,-1.1*sign(minc)*minc]),max([0,1.1*sign(maxc)*maxc])])
        xlabel('Space distance',fontsize=8)
        ylabel('Time distance',fontsize=8)
        zlabel('Cross-covariance',fontsize=8)
        axis('square')

        if test==0:
            hold(False)

    ##%%%%%% Check if there are no NaN

    if len(find(isnan(concatenate(ds,dt))))!=0:
        print 'Warning : some space/time classes do not contain pairs of points'

    return ds,dt,c,o

def crossvario(c1,c2,z1,z2,cl,method,*options):
    """
    
    Translated by Dimitri D'Or - December 2004

    Status: not completed: plotting commands

    % crossvario                - single cross variogram estimation (Jan 1,2001)
    %
    % Single variogram or cross variogram estimation. To the opposite of
    % the vario.m function, crossvario.m is is able to deal with the case
    % where the values of the two variables are given for two sets of
    % partially different locations. When several variograms or cross
    % variograms have to be computed for variables known at the same set
    % of locations, using vario.m is computationally more efficient than
    % a repeated use of crossvario.m.
    %
    % SYNTAX :
    %
    % [d,v,o]=crossvario(c1,c2,z1,z2,cl,method,options);
    %
    % INPUT : 
    %
    % c1        n1 by d      matrix of coordinates for the locations where the
    %                        values of the first variable are known. A line
    %                        corresponds to the vector of coordinates at a location,
    %                        so the number of columns is equal to the dimension of
    %                        the space. There is no restriction on the dimension of
    %                        the space.
    % c2        n2 by d      matrix of coordinates for the locations where the values
    %                        of the second variable are known, using the same conventions
    %                        as for c1.
    % z1        n1 by 1      column vector of values for the first variable at the
    %                        coordinates specified in c1.
    % z2        n2 by 1      column vector of values for the second variable at the
    %                        coordinates specified in c2.
    % cl        nc+1 by 1    vector giving the limits of the distance classes that are
    %                        used for estimating the cross variogram. The distance classes
    %                        are open on the left and closed on the right. The lower limit
    %                        for the first class is >=0.
    % method    string       that contains the name of the method used for computing
    %                        the distances between pairs of locations. method='kron'
    %                        uses a Kronecker product, whereas method='loop' uses a loop
    %                        over the locations. Using the Kronecker product is faster
    %                        for a small number of locations but may suffer from memory
    %                        size limitations depending on the memory available as its
    %                        requires the storage of a distance matrix. The loop method
    %                        may be used whatever the number of data locations and must be
    %                        used if an Out of Memory error message is generated. Both
    %                        methods yield exactly the same estimates.
    % options   1 by 1 or 3  vector of optional parameters that can be used if default
    %                        values are not satisfactory (otherwise this vector can simply
    %                        be omitted from the input list of variables), where :
    %                        options(1) displays the estimated cross variogram if the value
    %                        is set to one (default value is 0),
    %                        options(2) and options(3) are the minimum and maximum values
    %                        for the angles to be considered, using the same conventions as
    %                        for the pairsplot.m function. Angles can only be specified for
    %                        planar coordinates, i.e. when the number of columns in c1 and
    %                        c2 is equal to two.
    %
    % OUTPUT :
    %
    % d         nc by 1      vector giving the sorted values of the mean distance separating
    %                        the pairs of points that belong to the same distance class. 
    % v         nc by 1      vector of estimated variogram or cross variogram values (if z1
    %                        and z2 are identical, the function computes the variogram).
    % o         nc by 1      vector giving the number of pairs of points that belong to the 
    %                        corresponding distance classes.
    %
    % NOTE :
    %
    % Note that using crossvario(c1,c2,z1,z2,cl,method) yields the same result
    % as using crossvario(c2,c1,z2,z1,cl,method) both for method='kron' and
    % method='loop'.
    """

    ## %%%%%% Select the subset of identical coordinates

    if not isinstance(method,str):
        return 'method should be a char string'

    index=findpairs(c1,c2)
    if len(index)==0:
        return 'There are no identical coordinates. Cannot compute cross-variogram'

    c=take(c1,index[:,0])
    z1=take(z1,index[:,0])
    z2=take(z2,index[:,1])

    ## %%%%%% Initialize the parameters

    cl=sort(cl)
    if cl(1)<0:
        return 'Minimum class distance must be >=0'
    
    n=c.shape[0]
    nc=len(cl)-1
    minim=cl[0]
    maxim=cl[nc]

    if not options:
        options[0]=0
        noptions=1
    else:
        noptions=len(options)


    if noptions==3:
        a=options[1]*2*pi/360
        b=options[2]*2*pi/360
        if (c1.shape[1]!=2 and c2.shape[1]!=2):
            return 'Angle limits are specified only for planar coordinates'
        if (a==b) or (min([a,b])<-pi/2) or (max([a,b])>pi/2):
            return 'Angle limits must be different and between or equal to -90 and 90'

    if method=='kron':   ##%%%%% Uses a Kronecker product for computing distances

        ## %%% Compute the distances

        unit=ones((n,1))
        dc=BMEmatlab.kron(unit,c)-BMEmatlab.kron(c,unit)
        if dc.shape[1]==1:
            dist=abs(dc)
        else:
            dist=sqrt(transpose(sum(transpose(dc**2))))

        ## %%% Compute the angles

        if noptions==3:
            finddc1null=find(dc[:,0]==0)
            finddc1notnull=find(dc[:,0]!=0)
            ang=zeros((dc.shape[0],1))
            put(ang,finddc1null,(pi/2)*sign(take(dc[1],finddc1null)))
            put(ang,finddc1notnull,atan(take(dc[1],finddc1notnull)/take(dc[0],finddc1notnull)))

        ##%%% Select couples for appropriate distances and angles

        cond=(dist>max([0,minim])) & (dist<=maxim)
        if noptions==3:
            conda=(ang>a)
            condb=(ang<=b)
            if a<b:
                cond=cond & (conda & condb)
            else:
                cond=cond & (conda | condb)
        dist=take(dist,cond)
        m=len(dist)
        if m==0:
            return 'No couples of values within the specified classes'

        ## %%% Compute the cross-variogram

        isclass=[[] for i in xrange(0,nc)]
        d=zeros((nc,1),Float32)*NaN
        o=zeros((nc,1),Float32)
        for k in xrange(nc):
            isclass[k]=find((dist>cl[k]) & (dist<=cl[k+1]))
            o[k]=len(isclass[k])/2
            if o[k]!=0:
                d[k]=sum(dist(isclass[k]))/(2*o[k])

        dz1=BMEmatlab.kron(unit,z1)-BMEmatlab.kron(z1,unit)
        dz2=BMEmatlab.kron(unit,z2)-BMEmatlab.kron(z2,unit)
        product=dz1*dz2
        product=take(product,cond)
        v=zeros((nc,1),Float32)*NaN
        for k in xrange(nc):
            if o[k]!=0:
                v[k]=sum(product(isclass[k]))/(4*o[k])

    else:                      ## %%%%% Uses a loop over the data for computing distances

        d=zeros((nc,1),Float32)
        o=zeros((nc,1),Float32)
        v=zeros((nc,1),Float32)

        for i in xrange(0,n):
            for j in xrange(0,n):
                dist=sqrt(sum((c[i,:]-c[j,:])**2))
                cond=(dist>max([0,minim]))&(dist<=maxim)
                if noptions==3:
                    dc=c[i,0:1]-c[j,0:1]
                    if dc[0]==0:
                        ang=(pi/2)*sign(dc[1])
                    else:
                        ang=atan(dc[1]/dc[0])
                    conda=(ang>a)
                    condb=(ang<=b)
                    if a<b:
                        cond=cond & (conda & condb)
                    else:
                        cond=cond & (conda | condb)
                if cond==1:
                    index=sum(dist>cl)
                    if (index>=1) & (index<=nc):
                        put(d,index,take(d,index)+dist)
                        put(o,index,take(o,index)+1)
                        for k in xrange(0,nv):
                            for l in xrange(k,nv):
                                put(v,index,take(v,index)+(z1(i)-z1(j))*(z2(i)-z2(j)))

        for i in xrange(nc):
            if o[i]==0:
                d[i]=NaN
                v[i]=NaN
            else:
                d[i]=d[i]/o[i]
                v[i]=v[i]/(2*o[i])

    ## %%%%%% display the computed cross-variogram if options(1)=1

    if options[0]==1:
        test=(ishold==True)
        minv=min(take(v,not isnan(v)))
        maxv=max(take(v,not isnan(v)))
        maxd=max(take(d,not isnan(d)))
        plot(d,v,'.')
        hold(True)
        text(d,v,str(o),fontfize='medium')
        set(gca,fontsize=6)
        axis([0, maxd, min([0,-1.1*sign(minv)*minv]), max([0,1.1*sign(maxv)*maxv])])
        plot([0, maxd],[0, 0],':')
        xlabel('Distance',fontsize=8)
        ylabel('Cross-variogram',fontsize=8)
        if test==False:
            hold(False)

    ## %%%%%% Check if there are no NaN


    if len(find(isnan(d)))!=0:
        print 'Warning : some distance classes do not contain pairs of points'
        
    return d,v,o


def crossvarioST(c1,c2,z1,z2,cls,clt,*options):
    """
    
    Translated by Dimitri D'Or - December 2004

    Status: not completed: plotting commands

    % crossvarioST              - space/time cross variogram estimation (Jan 1,2001)
    %
    % Single space/time variogram or cross variogram estimation.
    % The function can be used when the values of the two variables
    % are given for two sets of partially different space/time
    % locations. 
    %
    % SYNTAX :
    %
    % [ds,dt,v,o]=crossvarioST(c1,c2,z1,z2,cls,clt,options);
    %
    % INPUT : 
    %
    % c1        n1 by d+1    matrix of space/time coordinates for the locations
    %                        where the values of the first variable are known. A
    %                        line corresponds to the vector of space/time coordinates
    %                        at a location, so the number of columns is equal to the
    %                        dimension of the space plus one, where the last column
    %                        refers to the temporal coordinate. There is no restriction
    %                        on the dimension of the space.
    % c2        n2 by d+1    matrix of space/time coordinates for the locations where
    %                        the values of the second variable are known, using the same
    %                        conventions as for c1.
    % z1        n1 by 1      column vector of values for the first variable at the
    %                        coordinates specified in c1.
    % z2        n2 by 1      column vector of values for the second variable at the
    %                        coordinates specified in c2.
    % cls       ncs+1 by 1   vector giving the limits of the spatial distance classes that
    %                        are used for estimating the variogram or cross variogram. The
    %                        distance classes are open on the left and closed on the right.
    %                        The lower limit for the first class is >=0.
    % clt       nct+1 by 1   vector giving the limits of the temporal "distance" classes
    %                        that are used for estimating the variogram or cross variogram.
    %                        As for cls, the classes are open on the left and closed on the
    %                        right. The lower limit for the first class is >=0.
    % options   1 by 1 or 3  vector of optional parameters that can be used if default
    %                        values are not satisfactory (otherwise this vector can simply
    %                        be omitted from the input list of variables), where :
    %                        options(1) displays the estimated cross variogram if the value
    %                        is set to one (default value is 0),
    %                        options(2) and options(3) are the minimum and maximum values
    %                        for the angles to be considered, using the same conventions as
    %                        for the pairsplot.m function. Angles can only be specified for
    %                        planar coordinates, i.e. when the number of columns in c1 and
    %                        c2 is equal to three.
    %
    % OUTPUT :
    %
    % ds        ncs by nct   matrix giving the sorted values of the mean spatial distance
    %                        separating the pairs of points that belong to the same space/time
    %                        distance class. Each line of ds correspond to the same spatial
    %                        class, where each column of ds correspond to the same temporal class.
    % dt        ncs by nct   matrix giving the sorted values of the mean temporal distance
    %                        separating the pairs of points that belong to the same space/time
    %                        distance class, with same conventions as for ds.
    % v         ncs by nct   matrix of estimated space/time variogram or cross variogram
    %                        values (if z1 and z2 are identical, the function computes the
    %                        space/time variogram), with same dimensions as ds and dt.
    % o         ncs by nct   matrix giving the number of pairs of points that belong to the
    %                        corresponding space/time distance classes, with same dimensions
    %                        as ds and dt.
    %
    % NOTE :
    %
    % Due to the nature of the problem itself, it is generally expected that
    % a large number of locations are involved in the computation of the
    % space/time variogram or cross variogram. For that reason and to the
    % opposite of the crossvario.m function, only the loop method has been
    % implemented in the function.
    """

    ## %%%%%% Select the subset of identical coordinates

    ## index=findpairs(c1,c2);
    ## if isempty(index),
    ##   error('There are no identical coordinates. Cannot compute cross-variogram');
    ## end;
    ## c=c1(index(:,1),:);
    ## z1=z1(index(:,1));
    ## z2=z2(index(:,2));

    ## %%%%%% Initialize the parameters

    ## cls=sort(cls);
    ## clt=sort(clt);
    ## if (cls(1)<0)|(clt(1)<0),
    ##   error('Minimum space/time class distances must be >=0');
    ## end;
    ## n=size(c,1);
    ## ncs=length(cls)-1;
    ## nct=length(clt)-1;
    ## minims=cls(1);
    ## maxims=cls(ncs+1);
    ## minimt=clt(1);
    ## maximt=clt(nct+1);
    ## dim=size(c,2)-1;

    ## if nargin==6,
    ##   options(1)=0;
    ##   noptions=1;
    ## else
    ##   noptions=length(options);
    ## end;

    ## if noptions==3,
    ##   a=options(2)*2*pi/360;
    ##   b=options(3)*2*pi/360;
    ##   if dim~=2,
    ##     error('Angle limits are specified only for planar coordinates');
    ##   end;
    ##   if (a==b)|(min([a,b])<-pi/2)|(max([a,b])>pi/2),
    ##     error('Angle limits must be different and between or equal to -90 and 90');
    ##   end;
    ## end;

    ## %%%%% Uses a loop over the data for computing distances

    ## ds=zeros(ncs,nct);
    ## dt=zeros(ncs,nct);
    ## o=zeros(ncs,nct);
    ## v=zeros(ncs,nct);

    ## for i=1:n,
    ##   for j=i+1:n,
    ##     dists=sqrt(sum((c(i,1:dim)-c(j,1:dim)).^2));
    ##     distt=abs(c(i,dim+1)-c(j,dim+1));
    ##     conds=(dists>max([0 minims]))&(dists<=maxims);
    ##     condt=(distt>max([0 minims]))&(distt<=maxims);
    ##     cond=(conds & condt);
    ##     if noptions==3,
    ##       dc=c(i,1:2)-c(j,1:2);
    ##       if dc(1)==0,
    ##         ang=(pi/2)*sign(dc(2));
    ##       else
    ##         ang=atan(dc(2)/dc(1));
    ##       end;
    ##       conda=(ang>a);
    ##       condb=(ang<=b);
    ##       if a<b,
    ##         cond=cond & (conda & condb);
    ##       else
    ##         cond=cond & (conda | condb);
    ##       end;
    ##     end;
    ##     if cond==1,
    ##       indexs=sum(dists>cls);
    ##       indext=sum(distt>clt);
    ##       if ((indexs>=1) & (indexs<=ncs) & (indext>=1) & (indext<=nct)),
    ##         ds(indexs,indext)=ds(indexs,indext)+dists;
    ##         dt(indexs,indext)=dt(indexs,indext)+distt;
    ##         o(indexs,indext)=o(indexs,indext)+1;
    ##         v(indexs,indext)=v(indexs,indext)+(z1(i)-z1(j))*(z2(i)-z2(j));
    ##       end;
    ##     end;
    ##   end;
    ## end;

    ## for i=1:ncs,
    ##   for j=1:nct,
    ##     if o(i,j)==0,
    ##       ds(i,j)=NaN;
    ##       dt(i,j)=NaN;
    ##       v(i,j)=NaN;
    ##     else
    ##       ds(i,j)=ds(i,j)/o(i,j);
    ##       dt(i,j)=dt(i,j)/o(i,j);
    ##       v(i,j)=v(i,j)/(2*o(i,j));
    ##     end;
    ##   end;
    ## end;

    ## %%%%%% display the computed cross-variogram if options(1)=1

    ## if options(1)==1,
    ##   test=(ishold==1);
    ##   minv=min(v(:));
    ##   maxv=max(v(:));
    ##   posx=(cls(1:ncs)+cls(2:ncs+1))/2;
    ##   posy=(clt(1:nct)+clt(2:nct+1))/2;
    ##   maxx=max(ds(:));
    ##   maxy=max(dt(:));

    ##   subplot(1,2,1);
    ##     obj1=contour(posx,posy,v');
    ##     clabel(obj1,'FontSize',8);
    ##     set(gca,'FontSize',8);
    ##     axis([0 maxx 0 maxy]);
    ##     xlabel('Space distance','FontSize',8);
    ##     ylabel('Time distance','FontSize',8);
    ##     axis('square');
    ##   subplot(1,2,2);
    ##     mesh(posx,posy,v');
    ##     set(gca,'FontSize',8);
    ##     axis([0 maxx 0 maxy min([0;-1.1*sign(minv)*minv]) max([0;1.1*sign(maxv)*maxv])]);
    ##     xlabel('Space distance','FontSize',8);
    ##     ylabel('Time distance','FontSize',8);
    ##     zlabel('Cross-variogram','FontSize',8);
    ##     axis('square');

    ##   if test==0,
    ##     hold off;
    ##   end;
    ## end;

    ## %%%%%% Check if there are no NaN

    ## if length(find(isnan([ds;dt])))~=0,
    ##   disp('Warning : some space/time classes do not contain pairs of points');
    ## end;


    ## %%%%%% Select the subset of identical coordinates

    index=findpairs(c1,c2)
    if len(index)==0:
        return 'There are no identical coordinates. Cannot compute cross-variogram'

    c=take(c1,index[:,0])
    z1=take(z1,index[:,0])
    z2=take(z2,index[:,1])

    ## %%%%%% Initialize the parameters

    cls=sort(cls)
    clt=sort(clt)
    if (cls[0]<0) or (clt[0]<0):
        return 'Minimum space/time class distances must be >=0'
    n=c.shape[0]
    ncs=len(cls)-1
    nct=len(clt)-1
    minims=cls[0]
    maxims=cls[ncs+1]
    minimt=clt[0]
    maximt=clt[nct+1]
    dim=c.shape[1]-1

    if not options:
        options[0]=0
        noptions=1
    else:
        noptions=len(options)

    if noptions==3:
        a=options[1]*2*pi/360
        b=options[2]*2*pi/360
        if dim!=2:
            return 'Angle limits are specified only for planar coordinates'
        if (a==b) or (min([a,b])<-pi/2) or (max([a,b])>pi/2):
            return 'Angle limits must be different and between or equal to -90 and 90'

    ## %%%%% Uses a loop over the data for computing distances

    ds=zeros((ncs,nct),Float32)
    dt=zeros((ncs,nct),Float32)
    o=zeros((ncs,nct),Float32)
    v=zeros((ncs,nct),Float32)

    for i in xrange(n):
        for j in xrange(i,n):
            dists=sqrt(sum((c[i,:dim]-c[j,:dim])**2))
            distt=abs(c[i,dim+1]-c[j,dim+1])
            conds=(dists>max([0, minims])) and (dists<=maxims)
            condt=(distt>max([0, minims])) and (distt<=maxims)
            cond=(conds and condt)
            if noptions==3:
                dc=c[i,:2]-c[j,:2]
                if dc[0]==0:
                    ang=(pi/2)*sign(dc[1])
                else:
                    ang=atan(dc[1]/dc[0])
                conda=(ang>a)
                condb=(ang<=b)
                if a<b:
                    cond=cond & (conda & condb)
                else:
                    cond=cond & (conda | condb)
            if cond==1:
                indexs=sum(dists>cls)
                indext=sum(distt>clt)
                if ((indexs>=1) & (indexs<=ncs) & (indext>=1) & (indext<=nct)):
                    ds[indexs,indext]=ds[indexs,indext]+dists
                    dt[indexs,indext]=dt[indexs,indext]+distt
                    o[indexs,indext]=o[indexs,indext]+1
                    v[indexs,indext]=v[indexs,indext]+(z1[i]-z1[j])*(z2[i]-z2[j])

    for i in xrange(ncs):
        for j in xrange(nct):
            if o[i,j]==0:
                ds[i,j]=NaN
                dt[i,j]=NaN
                v[i,j]=NaN
            else:
                ds[i,j]=ds[i,j]/o[i,j]
                dt[i,j]=dt[i,j]/o[i,j]
                v[i,j]=v[i,j]/(2*o[i,j])

    ## %%%%%% display the computed cross-variogram if options(1)=1

    if options[0]==1:
        test=(ishold==True)
        minv=min(ravel(v))
        maxv=max(ravel(v))
        posx=(cls[:ncs]+cls[1:ncs+1])/2
        posy=(clt[:nct]+clt[1:nct+1])/2
        maxx=max(ravel(ds))
        maxy=max(ravel(dt))

        subplot(1,2,1)
        obj1=contour(posx,posy,transpose(v))
        clabel(obj1,fontsize=8)
        set(gca,fontsize=8)
        axis([0, maxx, 0, maxy])
        xlabel('Space distance',fontsize=8)
        ylabel('Time distance',fontsize=8)
        axis('square')

        subplot(1,2,2)
        mesh(posx,posy,transpose(v))
        set(gca,fontsize=8)
        axis([0, maxx, 0, maxy, min([0,-1.1*sign(minv)*minv]), max([0,1.1*sign(maxv)*maxv])])
        xlabel('Space distance',fontsize=8)
        ylabel('Time distance',fontsize=8)
        zlabel('Cross-variogram',fontsize=8)
        axis('square')

        if test==False:
            hold(False)

    ## %%%%%% Check if there are no NaN

    if len(find(isnan(concatenate((ds,dt)))))!=0:
        print 'Warning : some space/time classes do not contain pairs of points'

    return ds,dt,v,o


def decluster(c,model,param):
    """
    
    Translated by Dimitri D'Or - December 2004

    Status: debugged, not tested

    % decluster                 - spatial declustering weights (Jan 1,2001)
    %
    % Compute the spatial weights associated with values. These weights
    % reflect the information content of the values ; values which
    % have close coordinates will tend to exhibit similar values and 
    % will have smaller weights than isolated values. The weight
    % associated with an isolated value is thus higher than the weight
    % for a value located inside a cluster. Weights in the decluster.m
    % function are based on the kriging weights for estimating a constant
    % mean over the area. Negative kriging weights are taken in absolute
    % values and are normalized in order to sum to one. 
    %
    % SYNTAX :
    %
    % [w]=decluster(c,model,param);
    %
    % INPUT :
    %
    % c      n by d   matrix of coordinates for the variable. A line
    %                 corresponds to the vector of coordinates at a
    %                 location, so that the number of columns corresponds
    %                 to the dimension of the space. There is no restriction
    %                 on the dimension of the space.
    % model  string   that contains the name of the variogram or covariance
    %                 model which is used for the computation.
    % param  1 by k   vector of values for the parameters of model, according
    %                 to the conventions for the corresponding variogram or 
    %                 covariance model.
    %
    % OUTPUT :
    %
    % w      n by 1   vector of weights associated with the c coordinates.
    %                 Weights are positive and sum to one.
    %
    % NOTE :
    % 
    % For a detailed discussion about the coding of the model and param
    % variables for nested variogram or covariance models and for space/time
    % cases, the reader is referred to the detailed description given for
    % the kriging.m fuction. Only one variable can be used in this function.
    """

    n=c.shape[0]
    unit=ones((n,1))
    K=coord2K(c,c,model,param)
    K=concatenate((concatenate((K,unit),1),concatenate((transpose(unit),zeros((1,1))),1)))
    k=concatenate((zeros((n,1)),ones((1,1))))
    w=matrixmultiply(inverse(K),k)
    w=abs(w[:n])
    w=w/sum(w)

    return w

def designmatrix(c,order):
    """
    
    Translated by Dimitri D'Or - December 2004

    Status: debugged, not tested

    % designmatrix              - design matrix in a linear regression model (Jan 1,2001)
    %
    % Build the design matrix associated with a polynomial
    % mean of a given order in a linear regression model
    % of the form z=X*b+e.
    %
    % SYNTAX :
    %
    % [X,index]=designmatrix(c,order);
    %
    % INPUT :
    %
    % c       n by d       matrix of coordinates for the locations. A line
    %                      corresponds to the vector of coordinates at a
    %                      location, so the number of columns in c corresponds
    %                      to the dimension of the space. There is no restriction
    %                      on the dimension of the space.
    % order   scalar       order of the polynomial mean along the spatial axes
    %                      specified in c, where order>=0. When order=NaN, an empty
    %                      X matrix is returned.
    %
    % OUTPUT :
    %
    % X       n by k       design matrix, where each column corresponds to one
    %                      of the polynomial term, sorted in the first place 
    %                      with respect to the degree of the polynomial term,
    %                      and sorted in the second place with respect to the 
    %                      axis number. 
    %
    % index   1 or 2 by k  matrix associated with the columns of X. The first line
    %                      specifies the degree of the estimated polynomial term for
    %                      the corresponding column of X, and the second line specifies
    %                      the axis number to which this polynomial term belongs. The
    %                      axis are numbered according to the columns of c. E.g., the axis
    %                      2 corresponds to the second column of c. Note that the value 0
    %                      in the second line of index is associated with the polynomial
    %                      term of degree equal to 0 (i.e., the constant term) that is 
    %                      defined jointly for all the axes. In the singular case where c
    %                      is a column vector (i.e., the dimension of the space is equal
    %                      to 1), there is only one line for the index variable.
    %
    % NOTE :
    %
    % 1- It is also possible to process several variables at the same time
    % (multivariate case). It is needed to specify additionally tags in the
    % c matrix. These tags are provided as a vector of values that refers to
    % the variable, the values ranging from 1 to nv, where nv is the number
    % of variables. E.g., if there are 3 variables, the input index column vector
    % must be defined, and the elements in index are equal to 1, 2 or 3. The
    % c and index variables are grouped using the MATLAB cell array notation,
    % so that c={c, index}, is now the correct input variable. Using the same
    % logic, order is now a column vector specifying the order of the polynomial
    % mean for each variable. For the output variable index, there is an additional
    % first column that refers to the variable number associated with the
    % corresponding column of X.
    %
    % 2- For space/time data, the convention is that the last column of the c
    % matrix of coordinates corresponds to the time axis. Is is then possible to
    % specify a different order for the polynomial along the spatial axes and the
    % temporal axis. For the univariate case, order is a 1 by 2 vector, where
    % order(1) is the order of the spatial polynomial and order(2) is the order of
    % the temporal polynomial. For the multivariate case where nv different variables
    % are considered, order is a nv by 2 matrix, where the first and second columns
    % of order contain the order of the spatial and the temporal polynomial for
    % each of the nv variables, respectively. If in that case order is entered as
    % a 1 by 2 matrix, the same spatial order corresponding to order(1) will be used
    % for all the variables.
    """


    noindex= not isinstance(c,list)
    if noindex:
        [n,nd]=c.shape
    else:
        [n,nd]=c[0].shape
        nv=max(c[1])
        
    X=array([[] for i in xrange(0,n)])
    index=array([[] for i in xrange(0,n)])

    if noindex:
        if order.shape[1]==1:
            order=[order, order]
        if not (isnan(order[0]) & isnan(order[1])):
            X=ones((n,1))
            index=zeros((2,1))
        if not isnan(order[0]):
            for j in xrange(order[0]):
                for k in xrange(nd-1):
                    X=concatenate((X,c[:,k]**j),1)
                    index=concatenate((index,array([j,k])[:,NewAxis]),1)
        if not isnan(order[1]):
            for j in xrange(order[1]):
                X=concatenate((X,c[:,nd]**j),1)
                index=concatenate((index,array([j,nd])[:,NewAxis]),1)
        if (not len(index)==0) & (nd==1):
            index=index[0,:]
    else:
        if order.shape[0]==1:
            order=BMEmatlab.kron(order,ones((nv,1)))
        if order.shape[1]==1:
            order=[order, order]
        for l in xrange(nv):
            findvar=find(c[1]==l)
            if not len(findvar)==0:
                if not (isnan(order[l,0]) & isnan(order[l,1])):
                    Xvar=zeros((n,1))
                    put(Xvar,findvar,1)
                    X=concatenate((X,Xvar),1)
                    index=concatenate((index,array([l,0,0])[:,NewAxis]),1)
                if not isnan(order[l,0]):
                    for j in xrange(order[l,0]):
                        for k in xrange(nd-1):
                            Xvar=zeros((n,1))
                            put(Xvar,findvar,take(c[0][:,k],findvar)**j)
                            X=concatenate((X,Xvar),1)
                            index=concatenate((index,array([l,j,k])[:,NewAxis]),1)
                if not isnan(order[l,1]):
                    for j in xrange(order[l,1]):
                        Xvar=zeros((n,1))
                        put(Xvar,findvar,take(c[0][:,nd],findvar)**j)
                        X=concatenate((X,Xvar),1)
                        index=concatenate((index,array([l,j,nd])[:,NewAxis]),1)
        if (not len(index)==0) & (nd==1):
            index=index[0:2,:]

    return X,index


def entropy(z,pdf):
    """
    
    Translated by Dimitri D'Or - December 2004

    Status: debugged, not tested

    % entropy                   - Shannon's entropy of a pdf  (Jan 1,2001)
    %
    % Compute the Shannon's entropy value of a distribution from
    % a discrete definition of the probability distribution function.
    %
    % SYNTAX :
    %
    % [e]=entropy(z,pdf);
    %
    % INPUT :
    %
    % z      n by 1   vector of values.
    % pdf    n by 1   vector of values for the probability distribution
    %                 function at the z values.
    %
    % OUTPUT :
    %
    % e      scalar   Shannon's entropy of the distribution.
    """

    index=argsort(z)
    z=sort(z)
    pdf=take(pdf,index)

    cond=(pdf!=0)
    z=take(z,cond)
    pdf=take(pdf,cond)

    n=len(z)
    zmin=min(z)
    zmax=max(z)

    e=trapezint(z,-log(pdf)*pdf,zmin,zmax)

    return e

def exponentialcdf(z,param):
    """
    
    Translated by Dimitri D'Or - December 2004

    Status: debugged, not tested

    % exponentialcdf            - exponential cumulative distribution function (Jan 1,2001)
    %
    % Compute the values of the cumulative distribution
    % function for an exponential distribution with
    % specified parameter.
    %
    % SYNTAX :
    %
    % [p]=exponentialcdf(z,param);
    %
    % INPUT :
    %
    % z        n by k   matrix of values for which the probability
    %                   distribution function must be computed.
    % param    scalar   parameter of the exponential distribution, 
    %                   where param>0.
    %
    % OUTPUT :
    %
    % p        n by k   matrix of values for the cumulative probabilities
    %                   computed at the corresponding z values.
    """

    if param<=0:
        return 'parameter for the exponential distribution must be positive'

    p=1-exp(-param*z)
    put(p,z<0,0)

    return p

def exponentialinv(p,param):
    """

    Translated by Dimitri D'Or - December 2004

    Status: debugged, not tested

    % exponentialinv            - inverse exponential cumulative distribution function (Jan 1,2001)
    %
    % Compute the quantiles for an exponential distribution
    % with specified parameter.
    %
    % SYNTAX :
    %
    % [z]=exponentialinv(p,param);
    %
    % INPUT :
    %
    % p        n by k   matrix of values for the cumulative probabilities.
    % param    scalar   parameter of the exponential distribution, 
    %                   where param>0.
    %
    % OUTPUT :
    %
    % z        n by k   matrix of quantile values computed at the
    %                   corresponding p values.
    """

    if param<=0:
        return 'parameter for the exponential distribution must be positive'

    if max(ravel(p))>1 | min(ravel(p))<0:
        return 'cumulative probability values must be inside [0,1]'

    z=-log(1-p)/param

    return z


def exponentialpdf(z,param):
    """

    Translated by Dimitri D'Or - December 2004

    Status: debugged, not tested

    % exponentialpdf            - exponential probability distribution function (Jan 1,2001)
    %
    % Compute the values of the probability distribution
    % function for an exponential distribution with
    % specified parameter.
    %
    % SYNTAX :
    %
    % [pdf]=exponentialpdf(z,param);
    %
    % INPUT :
    %
    % z        n by k   matrix of values for which the probability
    %                   distribution function must be computed.
    % param    scalar   parameter of the exponential distribution, 
    %                   where param>0.
    %
    % OUTPUT :
    %
    % pdf      n by k   matrix of values for the probability distribution
    %                   function computed at the corresponding z values.
    """

    if param<0:
        return 'the parameter for the exponential distribution cannot be negative'

    pdf=param*exp(-param*z)
    put(pdf,z<0,0)

    return pdf


def exponentialstat(param):
    """

    Translated by Dimitri D'Or - December 2004

    Status: debugged, not tested

    % exponentialstat           - Mean, variance and median of the exponential distribution (Jan 1,2001)
    %
    % Return the mean,variance and median for an exponential
    % distribution with specified parameter.
    %
    % SYNTAX :
    %
    % [m,v,q50]=exponentialstat(param);
    %
    % INPUT :
    %
    % param    scalar   parameter of the exponential distribution, 
    %                   where param>0.
    %
    % OUTPUT :
    %
    % m        scalar   mean of the distribution.
    % v        scalar   variance of the distribution.
    % q50      scalar   mediane of the distribution.
    """
    
    m=1/param
    v=1/(param**2)
    q50=-log(0.5)/param

    return m,v,q50

## def fminsmodelfit(param,d,v,o,model,nc,nm,np):
##     """

##     Translated by Dimitri D'Or - December 2004

##     Status: debugged, not tested

##     % fminsmodelfit             - fminsearch subroutine for modelfit.m (Jan 1,2001)
##     %
##     % Objective function used by the modelfit.m function for
##     % minimizing a weighted sum of squares (see modelfit.m).
##     %
##     % SYNTAX :
##     %
##     % [f]=fminsmodelfit(param,d,v,o,model,nc,nm,np);
##     """

##     vtheor=zeros((nc,1))
##     index=0
##     for i in xrange(nm):
##         parami=param[index:index+np[i]]
##         vtheor=vtheor+model[i](d,parami)
##         index=index+np[i]

##     f=dot(o,(((v-vtheor)/vtheor)**2))

##     return f

def fminsmodelfit(param,d,v,o,model,nc,nm,np):
    '''
    
    Translated by Didrik Pinte - July 2004

    Status: debugged, not tested
   
    
    % fminsmodelfit             - fmins subroutine for modelfit.m (Jan 1,2001)
    %
    % Objective function used by the modelfit.m function for
    % minimizing a weighted sum of squares (see modelfit.m).
    %
    % SYNTAX :
    %
    % [f]=fminsmodelfit(param,d,v,o,model,nc,nm,np);
    '''

    vtheor=zeros(nc,Float32)
    index=0
    for i in xrange(0,nm):
        parami=param[index:index+np[i]].tolist()
        out = ravel(model[i](d,parami))
        vtheor = vtheor + out
        index=index+np[i]

    temp = (((v-vtheor)/vtheor)**2)
    f= dot(o,temp)
    return f

def fmintriangularinv(q,param,p):
    """

    Translated by Dimitri D'Or - December 2004

    Status: debugged, not tested

    % fmintriangularinv         - fmin subroutine for triangularinv.m (Jan 1,2001)
    %
    % Objective function used by the triangularinv.m function
    % for identifying the value of the quantile in a triangular
    % distribution and for a given probability value
    % (see triangularinv.m).
    %
    % SYNTAX :
    %
    % [f]=fmintriangularinv(q,param,p);
    """
    
    f=(p-triangularcdf(q,param))**2

    return f


def gauss2other(z,yfile,cdfyfile,method='linear'):
    """

    Translated by Dimitri D'Or - December 2004

    Status: debugged, not tested

    % gauss2other               - transform from a Gaussian pdf to an arbitrary pdf (Jan 1,2001)
    %
    % Do the monotonic transformation which is reciprocal to the
    % transformation made by other2gauss.m. The function tranforms
    % a vector of zero mean unit variance Gaussian distributed
    % values into a vector of arbitrary distributed values.
    %
    % SYNTAX :
    %
    % [y]=gauss2other(z,yfile,cdfyfile,method);
    %
    % INPUT :
    %
    % z          n by 1   vector of zero mean unit variance Gaussian
    %                     distributed values.
    % yfile      k by 1   vector of values used to define the cumulative
    %                     distribution function for the y values.
    % cdfyfile   k by 1   vector of the cumulative distribution function
    %                     values at the yfile values.
    % method     string   which is optional and specifies the interpolation
    %                     method to be used. See interp1.m for available
    %                     methods. Default value is 'linear'.
    %
    % OUTPUT :
    %
    % y          n by 1   vector of non Gaussian transformed values associated
    %                     with the z values.
    %
    % NOTE :
    %
    % 1- For the z values associated with values for the cumulative
    % distribution function that are outside the definition
    % given by (yfile,cdfyfile), the corresponding output y values
    % are coded as NaN's.
    %
    % 2- It is also possible to process several variables at the same time
    % (multivariate case). It is then needed to specify additional tags
    % for the z values. These tags are provided as a vector of values that
    % refer to the variable, the values ranging from 1 to nv, where nv is
    % the number of variables. E.g., if there are 3 variables, an indexz
    % column vector must be defined, having same number of elements than z
    % (the function will also accepts a single integer value if all the
    % values belong to the same variable). The z and indexz vectors are
    % grouped using the MATLAB cell array notation, so that z={z,indexz}
    % is now the correct input variable. The yfile and cdfyfile vectors are
    % now cell arrays too, where each cell corresponds to the appropriate
    % vector of values for the corresponding variable.
    """

    if not isinstance(z,list):
        cdfz=gausscdf(z,[0, 1])
        finterp1d=interp1d(cdfyfile,yfile,method)
        y=finterp1d(cdfz)
    else:
        index=z[1]
        if len(index)==1:
            cdfz=gausscdf(z[0],[0, 1])
            finterp1d=interp1d(cdfyfile[index],yfile[index],method)
            y=finterp1d(cdfz)
        else:
            n=len(z[0])
            y=zeros((n,1))*NaN
            nv=max(index)
            for i in xrange(nv):
                indexi=find(index==i)
                cdfz=gausscdf(take(z[0],indexi),[0, 1])
                finterp1d=interp1(cdfyfile[i],yfile[i],method)
                put(y,indexi,finterp1d(cdfz))

    return y


def gaussbicdf(z1,z2,param,*options):
    """

    Translated by Dimitri D'Or - December 2004

    Status: not completed: use of mvnAG1

    % gaussbicdf                - bivariate Gaussian cumulative distribution function (Jan 1,2001)
    %
    % Compute the values of a bivariate Gaussian cumulative
    % distribution function for a set of couples of values.
    % The gaussbicdf.m function uses a FORTRAN77 subroutine
    % for evaluating the corresponding double integrals. 
    %
    % SYNTAX :
    %
    % [cdf]=gaussbicdf(z1,z2,param,options);
    %
    % INPUT :
    %
    % z1       n by k   matrix of values for the first variable.
    % z2       n by k   matrix of values for the second variable, where each
    %                   element in z1 is associated with the corresponding
    %                   element in z2. Either z1 or z2 can be scalar. In such
    %                   a case, the corresponding variable is kept constant and
    %                   the output variable has the same size than the other matrix.
    % param    1 by 5   vector of parameters for the bivariate distribution, where :
    %                   param(1) and param(2) are the mean of the first and second
    %                   variables,
    %                   param(3) and param(4) are the variances of the first and second
    %                   variables,
    %                   param(5) is the correlation coefficient.
    % options  1 by 2   vector of optional parameters that can be used if default
    %                   values are not satisfactory (otherwise this vector can simply
    %                   be omitted from the input list of variables), where :
    %                   options(1) is the maximum number of evaluations for the integral
    %                              (default value is 50000),
    %                   options(2) is the relative error on the estimation of the integral
    %                              (default value is 1e-4).
    %
    % OUTPUT :
    %
    % cdf      n by k   matrix of values for the cumulative distribution function computed
    %                   at the couples of values specified in z1 and z2.
    """

    if not options:
        options[0]=50000
        options[1]=1e-4

    sz1=z1.shape
    sz2=z2.shape
    if sum(sz1)==2:
        z1=ones(sz2)*z1
    if sum(sz2)==2:
        z2=ones(sz1)*z2
    [n,k]=z1.shape

    m1=param[0]
    m2=param[1]
    v1=param[2]
    v2=param[3]
    rho=param[4]

    z1=(z1-m1)/sqrt(v1)
    z2=(z2-m2)/sqrt(v2)
    C=array([[1, rho],[rho, 1]])

    if ((rho!=1) & (rho!=-1) & (rho!=0)):
        cdf=zeros((n,k))
        for i in xrange(n):
            for j in xrange(k):
                p=mvnAG1([-15,-15],[z1[i,j],z2[i,j]],C,options[0],0,options[1])
                cdf[i,j]=p
    else:
        if rho==0:         ## %%% independance case
            cdf=gausscdf(z1,[0, 1])*gausscdf(z2,[0, 1])
        if rho==1:         ## %%% degenerate case : z2=z1
            z=min(z1,z2)
            cdf=gausscdf(z,[0, 1])
        if rho==-1:        ## %%% degenerate case : z2=-z1
            cdf=zeros(z1.shape)
            cond=(z2>-z1)
            zmax=max(take(z2,cond),take(z1,cond))
            zmin=min(take(z2,cond),take(z1,cond))
            put(cdf,cond,gausscdf(zmax,[0, 1])-gausscdf(zmin,[0, 1]))

    return cdf


def gaussbipdf(z1,z2,param):
    """

    Translated by Dimitri D'Or - December 2004

    Status: debugged, not tested

    % gaussbipdf                - bivariate Gaussian probability distribution function (Jan 1,2001)
    %
    % Compute the values of the bivariate Gaussian probability
    % distribution function for a set of couples of values. 
    %
    % SYNTAX :
    %
    % [pdf]=gaussbipdf(z1,z2,param);
    %
    % INPUT :
    %
    % z1       n by k   matrix of values for the first variable.
    % z2       n by k   matrix of values for the second variable, where each
    %                   element in z1 is associated with the corresponding
    %                   element in z2. Either z1 or z2 can be scalar. In such
    %                   a case, the corresponding variable is kept constant and
    %                   the output variable has the same size than the other matrix.
    % param    1 by 5   vector of parameters for the bivariate distribution, where :
    %                   param(1) and param(2) are the mean of the first and second
    %                   variables,
    %                   param(3) and param(4) are the variances of the first and second
    %                   variables,
    %                   param(5) is the correlation coefficient.
    %
    % OUTPUT :
    %
    % pdf      n by k   matrix of values for the probability distribution function
    %                   computed at the couples of values specified in z1 and z2.
    """

    m1=param[0]
    m2=param[1]
    v1=param[2]
    v2=param[3]
    rho=param[4]

    z1=(z1-m1)/sqrt(v1)
    z2=(z2-m2)/sqrt(v2)

    if abs(rho)!=1:
        K=1/(2*pi*sqrt(v1*v2)*sqrt(1-rho**2))
        h=(1/(1-rho**2))*(z1c**2-2*rho*z1c*z2c+z2c**2)
        pdf=K*exp(-h/2)
    else:
        return 'correlation coefficient cannot be equal to +/- 1'

    return pdf


def gausscdf(z,param):
    """

    Translated by Dimitri D'Or - December 2004

    Status: debugged, not tested especially for erf

    % gausscdf                  - Gaussian cumulative distribution function (Jan 1,2001)
    %
    % Compute the values of the cumulative distribution
    % function for a Gaussian distribution with
    % specified mean and variance parameters.
    %
    % SYNTAX :
    %
    % [p]=gausscdf(z,param);
    %
    % INPUT :
    %
    % z        n by k   matrix of values for which the probability
    %                   distribution function must be computed.
    % param    1 by 2   parameters of the Gaussian distribution, where :
    %                   param(1) is the mean of the distribution,
    %                   param(2) is the variance of the distribution.
    %
    % OUTPUT :
    %
    % p        n by k   matrix of values for the cumulative probabilities
    %                   computed at the corresponding z values.
    """

    m=param[0]
    v=param[1]

    if v<0:
        return 'a variance cannot be negative'

    if v!=0:
        z=(z-m)/sqrt(v)
        p=((special.erf(z/sqrt(2))+1)/2)
    else:
        p=zeros(z.shape)
        index=find(z>=m)
        put(p,index,1)

    return p


def gaussinv(p,param):
    """

    Translated by Dimitri D'Or - December 2004

    Status: debugged, not tested especially for erfinv

    % gaussinv                  - inverse Gaussian cumulative distribution function (Jan 1,2001)
    %
    % Compute the quantiles for a Gaussian distribution with
    % specified mean and variance parameters.
    %
    % SYNTAX :
    %
    % [z]=gaussinv(p,param);
    %
    % INPUT :
    %
    % p        n by k   matrix of values for the cumulative probabilities.
    % param    1 by 2   parameters of the Gaussian distribution, where :
    %                   param(1) is the mean of the distribution,
    %                   param(2) is the variance of the distribution.
    %
    % OUTPUT :
    %
    % z        n by k   matrix of quantile values computed at the
    %                   corresponding p values.
    %
    % NOTE :
    %
    % For values of p equal to 0 and 1, instead of the +/- Inf values,
    % the function returns the values param(1) +/- 8.3*sqrt(param(2)).
    """

    m=param[0]
    v=param[1]

    if v<0:
        return 'a variance cannot be negative'

    if max(ravel(p))>1 | min(ravel(p))<0:
        return 'cumulative probability values must be inside [0,1]'

    z=sqrt(2)*erfinv(2*p-1)
    findisinf=find(isinf(z))
    signisinf=sign(take(z,findisinf))
    put(z,findisinf,signisinf*8.3)
    z=z*sqrt(v)+m

    return z


def gaussmultcdf(z,m,C,*options):
    """

    Translated by Dimitri D'Or - December 2004

    Status: not completed: use of mvnAG1

    % gaussmultcdf              - multivariate Gaussian cumulative distribution function (Jan 1,2001)
    %
    % Compute the value of the multivariate Gaussian cumulative
    % distribution functions for a vector of values for these
    % variables. The gaussmultcdf.m function uses a FORTRAN77
    % subroutine for evaluating the corresponding multiple integrals.
    %
    % SYNTAX :
    %
    % [cdf]=gaussmultcdf(z,m,C,options);
    %
    % INPUT :
    %
    % z        n by 1   vector of values for the Gaussian distributed variables.
    % m        n by 1   vector of means for the variables.
    % C        n by n   symmetric covariance matrix for the variables.
    % options  1 by 2   vector of optional parameters that can be used if default
    %                   values are not satisfactory (otherwise this vector can simply
    %                   be omitted from the input list of variables), where :
    %                   options(1) is the maximum number of evaluations for the integral
    %                              (default value is 50000),
    %                   options(2) is the relative error on the estimation of the integral
    %                              (default value is 1e-4).
    %
    % OUTPUT :
    %
    % cdf      scalar   value of the cumulative distribution function computed
    %                   for the z vector of values.
    """

    if not options:
        options[0]=50000
        options[1]=1e-4

    D=diag(diag(C))
    z=(z-m)/sqrt(diag(D))
    C=(D**(-0.5))*C*(D**(-0.5))
    n=len(z)
    unit=ones((n,1))
    cdf=mvnAG1.sadmvn(-15*unit,z,C,options[0],0,options[1])

    return cdf


def gaussmultpdf(z,m,C):
    """

    Translated by Dimitri D'Or - December 2004

    Status: debugged, not tested

    % gaussmultpdf              - multivariate Gaussian probability distribution function (Jan 1,2001)
    %
    % Compute the value of the multivariate Gaussian probability
    % distribution function for a vector of values for these
    % variables. 
    %
    % SYNTAX :
    %
    % [pdf]=gaussmultpdf(z,m,C);
    %
    % INPUT :
    %
    % z      n by 1   vector of values for the Gaussian distributed variables.
    % m      n by 1   vector of means for the variables.
    % C      n by n   symmetric covariance matrix for the variables.
    %
    % OUTPUT :
    %
    % pdf    scalar   value of the probability distribution function computed
    %                 for the z vector of values.
    """

    n=z.shape[0]
    A=((2*pi)**(n/2))*sqrt(linalg.det(C))
    pdf=exp(-0.5*dot(dot((z-m),inv(C)),(z-m)))/A

    return pdf


def gausspdf(z,param):    
    '''
    
    Translated by Didrik Pinte - July 2004

    Status: debugged, not tested
    
    Verified and tested
    % gausspdf                  - Gaussian probability distribution function (Jan 1,2001)
    %
    % Compute the values of the probability distribution
    % function for a Gaussian distribution with
    % specified mean and variance parameters.
    %
    % SYNTAX :
    %
    % [pdf]=gausspdf(z,param);
    %
    % INPUT :
    %
    % z        n by k   matrix of values for which the probability
    %                   distribution function must be computed.
    % param    1 by 2   parameters of the Gaussian distribution, where :
    %                   param(1) is the mean of the distribution,
    %                   param(2) is the variance of the distribution.
    %
    % OUTPUT :
    %
    % pdf      n by k   matrix of values for the probability distribution
    %                   function computed at the corresponding z values.
    '''
    
    m=param[0]
    v=param[1]

    if v < 0 :
        error('a variance cannot be negative')

    if v != 0:
        A=1/sqrt(2*pi*v)
        pdf = A * exp( -0.5 * ((z-m)/sqrt(v))**2)
    else:
        pdf = zeros((z.shape), Float32)
        index=find(z==m)
        pdf[index]=Inf
    
    return pdf


def gaussplot(z,*w):
    """

    Translated by Dimitri D'Or - December 2004

    Status: debugged, not tested  ## How to handle figure plotting when no outputs are set to the function?

    % gaussplot                 - Gaussian probability plot (Jan 1,2001)
    %
    % Plot values versus the corresponding quantiles of a zero
    % mean unit variance Gaussian distribution. Quantiles values
    % are obtained from an estimate of the cumulative distribution
    % function. For Gaussian distributed values, the graph should
    % appear close to a straigth line. 
    %
    % SYNTAX :
    %
    % [zq,qgauss]=gaussplot(z,w); 
    %
    % INPUT :
    %
    % z        n by 1   vector of values.
    % w        n by 1   optional vector of weights associated with the
    %                   values in z. These weights must be positive and
    %                   sum to one. If w is not specified, the weights are
    %                   all taken as equal. Weights may come, e.g., from a
    %                   declustering technique (see, e.g., decluster.m).
    %
    % OUTPUT :
    %
    % zq       m by 1   vector of sorted z values where the duplicate z values
    %                   have been removed (m<=n).
    % qgauss   m by 1   vector of the Gaussian quantile values associated with
    %                   the values in zq.
    %
    % NOTE :
    %
    % When output variables are specified, the graphic is not displayed and
    % gaussplot.m simply returns the values for these variables instead. 
    """

    if not w:
        (zq,cdfzq)=cdfest(z)
    else:
        (zq,cdfzq)=cdfest(z,w)

    qgauss=gaussinv(cdfzq,[0,1])

    ##if nargout==0:        ## How to handle figure plotting when no outputs are set to the function?
    plot(zq,qgauss,'o')
    xlabel('Values')
    ylabel('Gaussian quantiles')
    ##end;

    return zq,qgauss


def gaussstat(param):
    """

    Translated by Dimitri D'Or - December 2004

    Status: debugged, not tested

    % gaussstat                 - Mean, variance and median of the Gaussian distribution (Jan 1,2001)
    %
    % Return the mean,variance and median for a Gaussian
    % distribution with specified mean and variance parameters.
    %
    % SYNTAX :
    %
    % [m,v,q50]=gaussstat(param);
    %
    % INPUT :
    %
    % param    1 by 2   parameters of the Gaussian distribution, where :
    %                   param(1) is the mean of the distribution,
    %                   param(2) is the variance of the distribution.
    %
    % OUTPUT :
    %
    % m        scalar   mean of the distribution.
    % v        scalar   variance of the distribution.
    % q50      scalar   mediane of the distribution.
    """

    m=param[0]
    v=param[1]
    q50=param[0]

    return m,v,q50


def histline(z,nbins,bounds=None,w=None):
    """

    Translated by Dimitri D'Or - December 2004

    Status: debugged, not tested

    % histline                  - density-scaled histogram plotted using a line (Jan 1, 2001)
    %
    % Plot with line a density scaled histogram using regularly spaced bins
    % between specified minimum and maximum values. The density
    % scaled histogram is an histogram that has been rescaled so
    % that the histogram has the meaning of a probability distribution
    % function, i.e. the sum of the displayed rectangular surfaces is
    % equal to 1. 
    %
    % SYNTAX :
    %
    % [hdle,zline,fline,n,x]=histline(z,nbins,bounds,w);
    %
    % INPUT :
    %
    % z        n by 1      column vector of values.
    % nbins    scalar      number of bins to be used for drawing the histogram.
    % bounds   1 by 2      optional vector of values for the lower and upper limits of
    %                      the first and last bins, respectively.
    %                      default is bounds=[min(z) max(z)];
    % w        n by 1      optional vector of weights associated with the values
    %                      in z (see decluster.m).
    %
    % OUTPUT :
    %
    % hdle     scalar      handle for the histogram line, obtained as follow: hdle=plot(zline,fline);
    % zline    1 by k      vector of z values for the histogram line
    % fline    1 by k      vector of f values for the histogram line
    % n        nbins by 1  vector for the values of the density scaled histogram.
    % x        nbins by 1  vector of the bins centers for which the values in n
    %                      have been computed.
    %
    % NOTE :
    %
    % 1- Not-a-Number (NaN) values for z are automatically stripped out when
    % drawing the histogram.
    %
    % 2- Values in z that are outside the interval specified by the bounds
    % vector are automatically counted in the first or last bin if they
    % are below the lower limit or above the upper limit, respectively.
    """

    z=ravel(z)
    if bounds==None:
        bounds=[min(z), max(z)]
    if w==None:
        w=ones((z.shape))/len(z)
    zmin=bounds[0]
    zmax=bounds[0]

    (frequency,zcenter)=histscaled(z,nbins,bounds,w)
    nbins=len(frequency)
    diffzcenter=diff(zcenter)
    zbegin[0]=zcenter[0]-diffzcenter[0]/2
    zbegin[1:nbins]=zcenter[1:nbins]-diffzcenter/2
    zend=zcenter[:nbins-1]+diffzcenter/2
    zend[nbins]=zcenter[nbins]+diffzcenter[nbins-1]/2
    indexBegin=xrange(0,2*nbins,2) + 1
    indexEnd=xrange(0,2*nbins,2) + 1
    zline[0]=zbegin[0]
    Cbin[0]=0
    fline[0]=0
    put(zline,indexBegin,zbegin)
    put(zline,indexEnd,zend)
    put(fline,indexBegin,frequency)
    put(fline,indexEnd,frequency)
    put(zline,2*nbins+1,zend[nbins])
    Cbin[2*nbins+1]=0
    fline[2*nbins+1]=0
    hdle=plot(zline,fline)

    n=frequency
    x=zcenter

    return hdle,zline,fline,n,x



def histscaled(z,nbins,bounds=None,w=None):
    '''
    
    Translated by Didrik Pinte - July 2004

    Status: debugged, not tested
    
    % histscaled                - density-scaled histogram  plotted using bars (Jan 1,2001)
    %
    % Plot with bars a density scaled histogram using regularly spaced bins
    % between specified minimum and maximum values. The density
    % scaled histogram is an histogram that has been rescaled so
    % that the histogram has the meaning of a probability distribution
    % function, i.e. the sum of the displayed rectangular surfaces is
    % equal to 1. 
    %
    % SYNTAX :
    %
    % [n,x]=histscaled(z,nbins,bounds,w);
    %
    % INPUT :
    %
    % z        n by 1      column vector of values.
    % nbins    scalar      number of bins to be used for drawing the histogram.
    % bounds   1 by 2      optional vector of values for the lower and upper limits of
    %                      the first and last bins, respectively.
    %                      default is bounds=[min(z) max(z)];
    % w        n by 1      optional vector of weights associated with the values
    %                      in z (see decluster.m).
    %
    % OUTPUT :
    %
    % n        nbins by 1  vector for the values of the density scaled histogram.
    % x        nbins by 1  vector of the bins centers for which the values in n
    %                      have been computed.
    %
    % NOTE :
    %
    % 1- Not-a-Number (NaN) values for z are automatically stripped out when
    % drawing the histogram.
    %
    % 2- Values in z that are outside the interval specified by the bounds
    % vector are automatically counted in the first or last bin if they
    % are below the lower limit or above the upper limit, respectively.
    %
    % 3- When output variables are specified, the graphic is not displayed
    % and histscaled.m simply returns the values for these variables instead.
    % The optional ouput variables have the same meaning than for the hist.m
    % function. '''

    if bounds==None:
        bounds=[amin(z),amax(z)]
    
    if w == None:
        w=ones(z.shape,float)/len(z)

    w = take(w,find(~ isnan(w)))
    z = take(z,find(~ isnan(z)))
    
    histClass = numpy.arange(bounds[0],bounds[1]+(bounds[1]-bounds[0])/float(nbins),(bounds[1]-bounds[0])/float(nbins))
    x= (histClass[1:nbins+1]+histClass[0:nbins])/2;
    n = zeros(x.shape,float);
    for i in xrange(0,nbins):
        index= find((z>histClass[i]) & (z<=histClass[i+1]))
        n[i] = sum(take(w,index))

    n[0] = n[0] + sum(compress(z<=histClass[0],w))
    n[nbins-1]=n[nbins-1] + sum(compress(z>histClass[nbins-1],w))
    n=n/(x[1]-x[0])

    #if True:
    pyplot.bar(x,n)
    
    return (n,x)

def kerneldensity(z,zfile,v,w=None):
    '''
    
    Translated by Didrik Pinte - July 2004

    Status: approved
    
    % kerneldensity             - pdf estimation using a Gaussian kernel (Jan 1,2001)
    %
    % Implementation of the traditional kernel method for estimating
    % a nonparameteric univariate probability distribution function
    % from a set of values. The Gaussian kernel is used in this
    % function ; it is characterized by a smoothing parameter that
    % corresponds to the variance of the Gaussian distribution.
    %
    % SYNTAX :
    %
    % [pdfzfile]=kerneldensity(z,zfile,v,w);
    %
    %
    % INPUT :
    %
    % z         n by 1   vector of values.
    % zfile     nk by 1  vector of values for which the density must be estimated.
    % v         scalar   variance of the Gaussian kernel.
    % w         n by 1   optional column vector of weights for the z values.
    %                    These weights must be positive and must sum to one
    %                    (see, e.g., decluster.m above). If w is not specified,
    %                    all the weights are taken as equal.
    %
    % OUTPUT :
    %
    % pdfzfile  nk by 1  vector of estimated values for the kernel smoothed probability
    %                    distribution function computed at the zfile values. 
    '''
    
    n=len(z)
    m=len(zfile)
    pdfzfile = zeros((m),typecode=Float32,savespace=1)

    if w==None:
        w=ones((n),Float32)/n

    for i in xrange(0,n):
        pdfzfile += w[i] * gausspdf(zfile,[z[i],v])

    return pdfzfile


def kernelregression(ck,c,z,v,order,options=0):
    """

    Translated by Dimitri D'Or - December 2004

    Status: debugged, not tested

    % kernelregression          - prediction using a Gaussian kernel regression method (Jan 1,2001)
    %
    % Implementation of the regression.m function in a moving
    % neighbourhood context. Regression is conducted locally
    % at a set of coordinates using a least squares estimation
    % procedure for a linear regression model, where the
    % deterministic part of the linear model is a polynomial of
    % arbitrary order. Instead of using ordinary least squares,
    % the function uses a diagonal covariance matrix, where the
    % variances are inversely proportional to the weights provided
    % by a Gaussian kernel, so that weights are monotonically 
    % decreasing with the distance from the estimation location. 
    %
    % SYNTAX :
    %
    % [zk]=kernelregression(ck,c,z,v,order,options); 
    %
    % INPUT :
    %
    % ck        nk by d   matrix of coordinates for the estimation locations.
    %                     A line corresponds to the vector of coordinates at
    %                     an estimation location, so the number of columns
    %                     corresponds to the dimension of the space. There is
    %                     no restriction on the dimension of the space.
    % c         n by d    matrix of coordinates for the locations of the values,
    %                     with the same conventions as for ck.
    % z         n by 1    vector of values at the c coordinates.
    % v         scalar    variance of the Gaussian kernel along the spatial axes.
    % order     scalar    order of the polynomial mean along the spatial axes
    %                     specified in c, where order>=0.
    % options   scalar    optional parameter that can be used if default value
    %                     is not satisfactory (otherwise this parameter can simply
    %                     be omitted from the input list of variables). options is
    %                     equal to 1 or 0, depending if the user wants or does not
    %                     want to display the order number of the location which is
    %                     currently processed by the function (default value is 0).
    %
    % OUTPUT :
    %
    % zk        nk by 1   vector of estimated mean values at the estimation locations.
    %                     A value coded as NaN means that no estimation has been
    %                     performed at that location due to the lack of available data.
    %
    % NOTE :
    %
    % 1- It is worth noting that when order=0, the kernelregression.m
    % function is doing a simple moving average estimation where the
    % weights are proportionally inverse to the values of the Gaussian kernel.
    %
    % 2- To the opposite of the regression.m function, it is only possible to
    % process a single variable using this function, but space/time data are
    % allowed. In that case, the convention is that the last column of the c
    % and ck matrices of coordinates corresponds to the time axis. The order
    % variable is then a 1 by 2 vector, where order(1) is the order of the
    % spatial polynomial and order(2) is the order of the temporal polynomial.
    % The v variable is then a 1 by 2 vector, where v(1) and v(2) are the
    % variances of the Gaussian kernel along the spatial axes and the temporal
    % axis, respectively.
    """

    nk=ck.shape[0]

    zk=zeros((nk,1))*NaN
    nmax=Inf
    isST=(len(v)==2)
    if isST==1:
        dmax[0]=4*sqrt(v[0])
        dmax[1]=4*sqrt(v[1])
        dmax[2]=1             ## % this value is arbitrary as all locations are kept in the neighbourhood
    else:
        dmax=4*sqrt(v)

    for i in xrange(nk):
        (csub,zsub,dsub,nsub)=neighbours(ck[i,:],c,z,nmax,dmax)
        if nsub>0:
            csub=csub-ck[i]
            if isST==1:
                w=gausspdf(dsub[:,0],[0, v[0]])*gausspdf(dsub[:,1],[0, v[1]])
            else:
                w=gausspdf(dsub,[0, v])
            w=1./w
            K=diag(w)
            (best,Vbest,dum1,dum2)=regression(csub,zsub,order,K)
            zk[i]=best[0]
        if options !=0:
            print i,'/',nk

    return zk

def kurtosis(Z):
    """

    Translated by Dimitri D'Or - December 2004

    Status: debugged, not tested

    % kurtosis                  - experimental kurtosis coefficient (Jan 1,2001)
    %
    % Compute the experimental kurtosis coefficients for
    % a set of variables. The kurtosis coefficient is a
    % measure of the ""flatteness"" of a distribution,
    % compared to a Gaussian distribution. The theoretical
    % kurtosis coefficient is defined in a way that it is
    % equal to zero for Gaussian distributed variables, and
    % lower or greater than zero for distributions which are
    % flatter or sharper than the Gaussian distribution,
    % respectively. 
    %
    % SYNTAX :
    %
    % [k]=kurtosis(Z); 
    %
    % INPUT :
    %
    % Z    n by nv    matrix of values for the different variables, where
    %                 each column corresponds to the values of one variable.
    %
    % OUTPUT :
    %
    % k    1 by nv    vector of estimated kurtosis coefficients.
    """

    [m,n]=Z.shape
    u1=mean(Z)
    u2=(std(Z)**2)
    dev=(Z-BMEmatlab.kron(ones((m,1)),u1))
    u4=sum(dev**4)/(m-1)
    k=(u4.__float__()/(u2**2))-3

    return k


def loggausscdf(z,param):
    """

    Translated by Dimitri D'Or - December 2004

    Status: debugged, not tested

    % loggausscdf               - log-Gaussian cumulative distribution function (Jan 1,2001)
    %
    % Compute the values of the cumulative distribution
    % function for a log-Gaussian distribution with
    % specified mean and variance parameters for the
    % associated Gaussian distribution.
    %
    % SYNTAX :
    %
    % [p]=loggausscdf(z,param);
    %
    % INPUT :
    %
    % z        n by k   matrix of values for which the probability
    %                   distribution function must be computed.
    % param    1 by 2   parameters of the Gaussian distribution, where :
    %                   param(1) is the mean of the distribution,
    %                   param(2) is the variance of the distribution.
    %
    % OUTPUT :
    %
    % p        n by k   matrix of values for the cumulative probabilities
    %                   computed at the corresponding z values.
    """
    m=param[0]
    v=param[1]

    if v<0:
        return 'a variance cannot be negative'

    isnegative=find(z<=0)
    z=where(z<=0,z,Nan)
    z=log(z)

    if v!=0:
        z=(z-m)/sqrt(v)
        p=((erf(z/sqrt(2))+1)/2)
        put(p,isnegative,0)
    else:
        p=zeros((z.shape))
        where(z>=m,z,1)

    return p


def loggaussinv(p,param):
    """

    Translated by Dimitri D'Or - December 2004

    Status: debugged, not tested

    % loggaussinv               - inverse log-Gaussian cumulative distribution function (Jan 1,2001)
    %
    % Compute the quantiles for a log-Gaussian distribution with
    % specified mean and variance parameters for the associated
    % Gaussian distribution.
    %
    % SYNTAX :
    %
    % [z]=loggaussinv(p,param);
    %
    % INPUT :
    %
    % p        n by k   matrix of values for the cumulative probabilities.
    % param    1 by 2   parameters of the Gaussian distribution, where :
    %                   param(1) is the mean of the distribution,
    %                   param(2) is the variance of the distribution.
    %
    % OUTPUT :
    %
    % z        n by k   matrix of quantile values computed at the
    %                   corresponding p values.
    %
    % NOTE :
    %
    % For values of p equal to 1, the function returns the values
    % exp(param(1)+8.3*sqrt(param(2))) instead of the +Inf value
    % (see gaussinv.m).
    """

    m=param[0]
    v=param[1]

    if v<0:
        return 'a variance cannot be negative'

    if max(ravel(p))>1 or min(ravel(p))<0:
        return 'cumulative probability values must be inside [0,1]'

    z=sqrt(2)*erfinv(2*p-1)
    findisinf=find(isinf(z))
    signisinf=sign(take(z,findisinf))
    put(z,findisinf,signisinf*8.3)
    z=z*sqrt(v)+m
    z=exp(z)
    if len(signisinf)!=0:
        index=find(signisinf==-1)
        put(z,take(findisinf,index),0)

    return z


def modelfit(d,v,o,model,param0,options=None):
    '''
    
    Translated by Didrik Pinte - July 2004

    Status: debugged, not tested
    
    % modelfit                  - single variogram/covariance least squares fitting (Jan 1,2001)
    %
    % Iterated non linear weighted least squares estimation
    % procedure for selecting the parameters of a single (cross)
    % variogram or covariance function model specified by the
    % user. Weights are taken as proportional to the number of
    % pairs of points that were used in the estimation for each
    % distance class.
    %
    % SYNTAX :
    %
    % [param]=modelfit(d,v,o,model,param0,options);
    %
    % INPUT :
    %
    % d          nc by 1      vector giving the sorted values of the mean distance separating
    %                         the pairs of points that belong to the same distance class. 
    % v          nc by 1      vector of estimated (cross) variogram or covariance values.
    % o          nc by 1      vector giving the number of pairs of points that belong to the 
    %                         corresponding distance classes.
    % model      string       that contains the name of the variogram or covariance model for
    %                         which parameters are sought (see the MODELS directory for possible
    %                         names of variograms and covariance functions). 
    % param0     1 by k       initial guess values for the parameters of model, according to the
    %                         convention for the corresponding variogram or covariance model.
    % options    1 by 1 or 4  vector of optional parameters that can be used if default values
    %                         are not satisfactory (otherwise this vector can simply be omitted
    %                         from the input list of variables), where :
    %                         options(1) displays the estimated and fitted model if the value is
    %                         set to one (default value is 0),
    %                         options(2) is the termination tolerance for the estimation of the
    %                         parameters (default value is 1e-4),  xtol
    %                         options(3) is the termination tolerance for the quadratic objective
    %                         function (default value is 1e-4), ftol
    %                         options(4) is the maximum number of iterations (default value is
    %                         400 * # of models). maxiter
    %
    % OUPUT :
    %
    % param      1 by k       estimated values of the parameters for model, according to a
    %                         weighted least square criterion.
    %
    % NOTE :
    %
    % For a detailed discussion about the coding of the model and param0
    % variables for nested models, the reader is referred to the detailed
    % description given for the kriging.m function.
    '''

    ##%%%%%% Initialize the parameters

    nc=len(d);
    if not isinstance(model,list):
        nm=1;
        np=len(param0)
        model=[model]
    else:
        nm=len(model)
        np=zeros(nm)
        for i in xrange(0,nm):
            if isinstance(param0[i],list):
                np[i]=len(param0[i])
            else:
                np[i]=1

    if options==None:
        options1=0
        xtol=1e-4
        ftol=1e-4
        maxiter=400*nm
    else :
        options1=options[0]
        xtol=options[1]
        ftol=options[2]
        maxiter=options[3]

     ##%%%%%% Reset param0 from cell array to vector 

    paramtemp=[]
    if nm > 1:
        for i in xrange(0,nm):
            if isinstance(param0[i],list):
                paramtemp+=param0[i]
            else : paramtemp.append(param0[i])
        param0=paramtemp

    ##%%%%%% Parameter estimation by least squares

    args = (d,v,o,model,nc,nm,np)
    (param, fopt, iter, funcalls, warnflag) = optimize.fmin(fminsmodelfit,param0,args,xtol,ftol,maxiter,full_output=1);
    #others = {fopt, iter, funcalls, warnflag, allvecs}

    ##%%%%%% Reset param from vector to list

    if nm>1 :
        paramtemp=[]
        index=0
        for i in xrange(0,nm):
            tparam = param[index:index+np[i]]
            if isinstance(tparam,list):
                paramtemp += tparam
            else:
                paramtemp.append(param[index:index+np[i]])
            index=index+np[i]
        param=paramtemp

    ##%%%%%% Display the fitted variogram/covariance 

    if options1==1 :
        print 'Computation stopped after ' + str(iter) + '/' + str(options[3]) + ' iterations'

        vtheor=0
        if nm==1:
            vtheor=vtheor+ model[1](d,param)
        else:
            for i in xrange(0,nm):
                vtheor=vtheor+ model[i](d,param[i])

        kwargs={'fontsize':6}
        plot(d,v,'*')
        plot(d,vtheor,)
        xlabel('Distance',{'fontsize':8})
        ylabel('Variogram/Covariance',{'fontsize':8})
        title('Estimated and fitted variograms/covariances',{'fontsize':8})
        minv = amin(ravel(array([v,vtheor])))
        maxv = amax(ravel(array([v,vtheor])))
        axis([0, amax(d), amin([0,-1.1*sign(minv)*minv]),amax([0,1.1*sign(maxv)*maxv])])
        plot(array([0,amax(d)]),array([0,0]),':')

    return param


def pdf2cdf(z,pdf):
    '''
    
    Translated by Didrik Pinte - July 2004

    Status: debugged, not tested
    
    % pdf2cdf                   - compute the cdf from the pdf (Jan 1,2001)
    %
    % Computes the values of the cumulative distribution function
    % based on a discrete definition of the corresponding probability
    % distribution function. The integration of the probability
    % distribution function is realized using a trapezoidal integration
    % formula. 
    %
    % SYNTAX :
    %
    % [cdf]=pdf2cdf(z,pdf);
    %
    % INPUT :
    %
    % z      n by 1   vector of values.
    % pdf    n by 1   vector of the probability distribution function values
    %                 at the z values.
    %
    % OUPUT :
    %          
    % cdf    n by 1   vector of the cumulative distribution function values
    %                 at the z values.
    %
    % NOTE :
    %
    % As the routine uses a trapezoidal integration scheme, it is recommended
    % to have a finely discretized definition of the probability distribution
    % function, with values that are equal to or close to zero for the lowest
    % z values, so that the neglected probability below the lowest value in z
    % is null or close to zero. The value of the cumulative distribution
    % function at the lowest z value is set equal to half the value of the
    % cumulative distribution function at the following z value.
    '''
    
    n=len(z)

    index1 = argsort(z)
    z = sort(z)                       ##%%% Keep track of the rank for the sorted values
    pdf=take(pdf,index1)              ##%%% sort the pdf values accordingly

    dz=diff(z)                        ##%%% Compute the z increments
    pdf=(pdf[0:n-1]+pdf[1:n])/2.0

    cdf = cumsum(dz*pdf)              ##%%% Integrates using trapezoidal formula

    temp = cdf.tolist()
    temp.insert(0,cdf[0]/2.0)
    cdf = array((temp))               ##%%% Add first value taken as half of the second one

    ##%%% Check if the maximum cdf value>1
    
    maxcdf = amax(cdf)                ##%%% maxcdf=MA.maximum(cdf)
    
    if maxcdf>1:                      ##%%% and reset the cdf values if so
        cdf=cdf/maxcdf

    index2 = argsort(index1)
    index1=sort(index1)               ##%%% Resort the values in the original order
    cdf=take(cdf,index2)

    return cdf


def vario(c,Z,cl,method,options=[1]):
    '''  
    
    Translated by Didrik Pinte - July 2004

    Status: debugged, tested (the kron method is not working correctly)    

    %
    % vario                     - multivariate variogram estimation (Jan 1,2001)
    %
    % Estimate the variograms and cross variograms for a set of
    % variables which are known at the same set of coordinates. 
    %
    % SYNTAX :
    %
    % [d,V,o]=vario(c,Z,cl,method,options);
    %
    % INPUT : 
    %
    % c         n by d       matrix of coordinates for the locations where the
    %                        values are known. A line corresponds to the vector
    %                        of coordinates at a location, so the number of columns
    %                        is equal to the dimension of the space. There is no
    %                        restriction on the dimension of the space.
    % Z         n by nv      matrix of values for the variables. Each line is
    %                        associated with the corresponding line vector of
    %                        coordinates in the c matrix, and each column corresponds
    %                        to a different variable.
    % cl        nc+1 by 1    vector giving the limits of the distance classes that are
    %                        used for estimating the variograms and cross variograms.
    %                        The distance classes are open on the left and closed on
    %                        the right. The lower limit for the first class is >=0.
    % method    string       that contains the name of the method used for computing
    %                        the distances between pairs of locations. method='kron'
    %                        uses a Kronecker product, whereas method='loop' uses a loop
    %                        over the locations. Using the Kronecker product is faster
    %                        for a small number of locations but may suffer from memory
    %                        size limitations depending on the memory available, as it
    %                        requires the storage of a distance matrix. The loop method
    %                        may be used whatever the number of data locations is and
    %                        must be used if an Out of Memory error message is generated.
    %                        Both  methods yield exactly the same estimates.
    % options   1 by 1 or 3  vector of optional parameters that can be used if default
    %                        values are not satisfactory (otherwise this vector can simply
    %                        be omitted from the input list of variables), where :
    %                        options(1) displays the estimated variograms if the value is
    %                        set to one (default value is 0),
    %                        options(2) and options(3) are the minimum and maximum values
    %                        for the angles to be considered, using the same conventions as
    %                        for the pairsplot.m function. Angles can only be specified for
    %                        planar coordinates, i.e. when the number of columns in c is
    %                        equal to two.
    %
    % OUTPUT :
    %
    % d         nc by 1      vector giving the sorted values of the mean distance separating
    %                        the pairs of points that belong to the same distance class. 
    % V         nv by nv     symmetric array of cells that contains the variograms and cross
    %                        variograms estimates for the distance classes specified in d.
    %                        Diagonal cells contain the nc by 1 vector of variogram estimates,
    %                        whereas off-diagonal cells contain the nc by 1 vector of cross
    %                        variogram estimates. If Z is a column vector (only one variable),
    %                        then V is simply a column vector having same size as d.
    % o         nc by 1      vector giving the number of pairs of points that belong to the 
    %                        corresponding distance classes.
    %
    % NOTE :
    %
    % 1- The d, V and o output variables can be used without modification as input
    % for the coregfit.m function.
    %
    % 2- When a distance class do not contain any pairs of points, the function
    % output a warning message. The d and V elements for the corresponding distance
    % class are thus coded as NaN s, whereas the corresponding o element is equal to 0.
    '''
    
    ##%%%%%% Initialize the parameters    
    
    if not method.isalpha():
      return 'method should be a char string'

    cl = sort(cl)
    if cl[0]<0:
      return 'Minimum class distance must be >=0'
    
    n = c.shape[0]
    nc = len(cl)
    minim = cl[0]
    maxim = cl[-1]
    try:
        nv = Z.shape[1]
    except IndexError:
        Z = Z[:,NewAxis]
        nv=Z.shape[1]

    V = zeros((nv,nv,nc-1),Float32)

    noptions = len(options)

    if noptions==3:
        a = options[1]*2*pi/360
        b = options[2]*2*pi/360
        if c.shape[1]!=2:
            print 'Angle limits are specified only for planar coordinates'
            return [None,None,None]
        if (a==b)|(amin([a,b])< -pi/2)|(amax([a,b])> pi/2):
            print 'Angle limits must be different and between or equal to -90 and 90'
            return [None,None,None]

    if method.find('kron')!=-1:
        print 'Using the kron method'
        ##%%%%% Uses a Kronecker product for computing distances
        ##%%% Compute the distances
        unit=ones([n,1],Float32)
        dc=BMEmatlab.kron(unit,c)-BMEmatlab.kron(c,unit);
        if dc.shape[1]==1:
            dist=abs(dc)
        else:
            dist=sqrt(sum(dc**2,1))


        ##print "Compute the angles"

        if noptions>=3:
            finddc1null=where(dc[:,1]==0,1,0)
      
            ang=zeros((dc.shape[0]),Float32);

            putmask(ang, finddc1null, (pi/2)*sign(dc[:,1]))
            putmask(ang, ~finddc1null, arctan(dc[:,1]/dc[:,0]))
       
        ##print "Select couples for appropriate distances and angles"

        cond = where((dist > amax([0,minim]))& (dist <= maxim),1,0)
        
        if noptions>=3:
            conda=(ang>a)
            condb=(ang<=b)
            if a < b:
                cond = cond & (conda & condb)
            else:
                cond = cond & (conda | condb)

        dist = compress(cond,dist)
        m = len(dist)
        if m==0:
            return 'No couples of values within the specified classes'

        ##print "Loop over the number of variables and compute (cross)variogram"

        isclass = [];
        d = zeros((nc-1),Float32)
        d[:] = NAN
        o = zeros((nc-1),Float32)
        for k in range(0,nc-1):
            isclass.append(find((dist>cl[k])&(dist<=cl[k+1])))
            o[k] = len(isclass[k])/2
            if not o[k]==0:
                d[k] = sum(take(dist,isclass[k]))/(2*o[k])

        for i in xrange(0,nv):
            for j in xrange(i,nv):
                zi = Z[:,i]
                zj = Z[:,j]
                dzi = BMEmatlab.kron(unit,zi[:,NewAxis]) - BMEmatlab.kron(zi[:,NewAxis],unit)
                dzj = BMEmatlab.kron(unit,zj[:,NewAxis])- BMEmatlab.kron(zj[:,NewAxis],unit)
                product = transpose(multiply(transpose(dzi),transpose(dzj)))
                product = compress(cond, ravel(product))
                v = zeros((nc-1),Float32)
                for k in xrange(0,nc-1):
                    if not o[k]==0:
                        comProd = compress(isclass[k],product)                                         
                        v[k]= sum(comProd)/(4*o[k])                     
                V[i,j]= v
                if i!=j:
                    V[j,i]=v
    
    else :
    ## %%%%% Uses a loop over the data for computing distances
        print "Using the loop method"
        d=zeros((nc-1),Float32)
        o=zeros((nc-1),Float32)
        for i in xrange(0,nv):
            for j in xrange(0,nv):
                v = zeros((nc-1),Float32)
                V[i,j] = v
        for i in xrange(0,n):
            for j in range(i+1,n):
                dist=sqrt(sum((c[i,:]-c[j,:])**2))
                cond = where(dist>amax([0,minim])&(dist<=maxim),1,0)
                if noptions>=3:
                    dc=c[i,:]-c[j,:]
                    
                    if dc[0] == 0 :
                        ang = (pi/2)*sign(dc[1])
                    else: 
                        ang = arctan(dc[1]/dc[0])
               
                    conda=(ang>a)
                    condb=(ang<=b)
                    if a<b:
                        cond=cond & (conda & condb)
                    else:
                        cond=cond & (conda | condb)
                        
                if cond==1:
                    index=sum(dist>cl)
                    if (index>=0) & (index<nc):
                        d[index-1]=d[index-1]+dist
                        o[index-1]=o[index-1]+1
                        for k in range(0,nv):
                            for l in range(k,nv):
                                V[k,l][index-1]=V[k,l][index-1]+(Z[i,k]-Z[j,k])*(Z[i,l]-Z[j,l])
        for i in xrange(0,nc-1):
          if o[i]==0:
              d[i]=NAN
              for j in xrange(0,nv):
                  for k in xrange(j,nv):
                      V[j,k][i]=NAN
                      V[k,j][i]=NAN

          else:
              d[i]=d[i]/o[i]
              for j in xrange(0,nv):
                  for k in xrange(j,nv):
                      V[j,k][i]=V[j,k][i]/(2*o[i])
                      V[k,j][i]=V[j,k][i]                  

    ##%%%%%% display the computed variograms if options(1)=1
    ## To be checked and corrected.
    if options[0]== 1:
        graph=0
        for i in xrange(0,nv):
            for j in xrange(i,nv):
                minVij=amin(V[i,j])
                maxVij=amax(V[i,j])
                subplot(nv,nv,(i-1)*nv+j)
                graph+=1
                digit = str(nv)+str(nv)+str(graph)
                #print digit
                subplot(digit)
                #print str(V[i,j].shape)
                plot(d,V[i,j],'o')
                #set(gca,'FontSize',6)
                axes = [0,amax(d),amin([0,-1.1*sign(minVij)*minVij]),amax([0,1.1*sign(maxVij)*maxVij])]
                #print str(mean(V[i,j]))
                #print(axes)
                axis(axes)
                plot([0,amax(d)[0]],[0,0],':')
                xlabel("Distance")
                ylabel("Variogram")
                if nv >1 :
                    title(("Couple "+ str(i+1)+"-"+str(j+1)))

    ##%%%%%% V is a vector if there is only one variable

    if nv==1:
        V=V[0,0]

    ##%%%%%% Check if there are no NaN

    if len(find(isnan(d))) > 0 :
      print('Warning : some distance classes do not contain pairs of points')

    return (d,V,o)
