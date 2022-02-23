from matplotlib.matlab import *
from scipy import *

def exponentialC(D,param):
    '''
    
    Translated by Dimitri D''Or - November 2004

    Status: debugged, not tested
    
    % exponentialC              - exponential covariance model (Jan 1,2001)
    %
    % Compute the exponential covariance matrix from a distance matrix.
    % 
    % SYNTAX :
    %
    % [C]=exponentialC(D,param);
    %
    % INPUT :
    %
    % D       n by m   distance matrix, having real values >=0.
    % param   1 by 2   vector of parameters such that param[0] is the sill
    %                  and param[1] is the range of the model (i.e. the
    %                  distance to reach 5% of the sill value).
    %
    % OUTPUT :
    %
    % C       n by m   covariance matrix with same size as D.
    '''

    C=param[0]*exp(-3*D/param[1])    
    return C


def exponentialV(D,param):

    '''
    
    Translated by Didrik Pinte - July 2004

    Status: debugged, not tested
    
    % exponentialV              - exponential variogram model (Jan 1,2001)
    %
    % Compute the exponential variogram matrix from a distance matrix.
    %
    % SYNTAX :
    %
    % [V]=exponentialV(D,param);
    %
    % INPUT :
    %
    % D       n by m   distance matrix, having real values >=0.
    % param   1 by 2   vector of parameters such that param[0] is the sill
    %                  and param[1] is the range of the model (i.e. the
    %                  distance to reach 95% of the sill value).
    %
    % OUTPUT :
    %
    % V       n by m   variogram matrix with same size as D.
    '''
    
    V = param[0]*(1-exp(-3*D/param[1]))
    return V


def gaussianC(D,param):
    """
    
    Translated by Dimitri D''Or - November 2004

    Status: debugged, not tested
    
    % gaussianC                 - Gaussian covariance model (Jan 1,2001)
    %
    % Compute the Gaussian covariance matrix from a distance matrix.
    % 
    % SYNTAX :
    %
    % [C]=gaussianC(D,param);
    %
    % INPUT :
    %
    % D       n by m   distance matrix, having real values >=0.
    % param   1 by 2   vector of parameters such that param[0] is the sill
    %                  and param[1] is the range of the model (i.e. the
    %                  distance to reach 5% of the sill value).
    %
    % OUTPUT :
    %
    % C       n by m   covariance matrix with same size as D.
    """

    C=param[0]*exp(-(sqrt(3)*D/param[1])**2)
    return C

def gaussianCST(Ds,Dt,param):
    """
    
    Translated by Dimitri D'Or - November 2004

    Status: debugged, not tested
    
    % gaussianCST               - Gaussian non separable S/T covariance model (Jan 1, 2001)
    %
    % Gaussian non separable space-time covariance model with a
    % space/time metric such that the space/time distance=
    % spatial distance+k*temporal distance, where k>0 refers to
    % the space/time metric. This function is provided as an example
    % of non separable space/time covariance model.
    %
    % SYNTAX :
    %
    % [C]=gaussianCST(Ds,Dt,param);
    %
    % INPUT :
    %
    % Ds     n by m   distance matrix in space, having real values >=0.
    % Dt     n by m   distance matrix in time with same dimensions as Ds
    % param  1 by 3   vector of parameters such that :
    %                 param[0] is the sill of the space-time model,
    %                 param[1] is the space/time range, i.e. the space/time
    %                          distance to reach 5% of the sill,
    %                 param[2] refers to the space/time metric, such that
    %                          space/time distance=Ds+param[2]*Dt.
    %
    % OUTPUT :
    %
    % C      n by m   space/time covariance matrix with same size as Ds and Dt.
    """

    C=gaussianC(Ds+param[2]*Dt,[param[0],param[1]])
    return C


def gaussianV(D,param):
    """
    
    Translated by Dimitri D'Or - November 2004

    Status: debugged, not tested
    
    % gaussianV                 - Gaussian variogram model (Jan 1,2001)
    %
    % Compute the Gaussian variogram matrix from a distance matrix.
    % 
    % SYNTAX :
    %
    % [V]=gaussianV(D,param);
    % 
    % INPUT :
    %
    % D       n by m   distance matrix, having real values >=0.
    % param   1 by 2   vector of parameters such that param[0] is the sill
    %                  and param[1] is the range of the model (i.e. the
    %                  distance to reach 95% of the sill value).
    %
    % OUTPUT :
    %
    % V       n by m   variogram matrix with size as D.
    """

    V=param[0]*(1-exp(-(sqrt(3)*D/param[1])**2))
    return V


def holecosC(D,param):
    """
  
    Translated by Dimitri D'Or - November 2004

    Status: debugged, not tested

    % holecosC                  - cosinusoidal hole effect covariance model (Jan 1,2001)
    %
    % Compute the cosinusoidal hole effect covariance matrix
    % from a distance matrix.
    %
    % SYNTAX :
    %
    % [C]=holecosC(D,param);
    %
    % INPUT :
    %
    % D       n by m   distance matrix, having real values >=0.
    % param   1 by 2   vector of parameters such that param[0] is the
    %                  half amplitude and param[1] is half the periodicity
    %                  of the hole effect.
    %
    % OUTPUT :
    %
    % C       n by m   covariance matrix with same size as D.
    %
    % NOTE :
    %
    % This model is only valid in a 1-dimensional space.
    """

    C=param[0]*cos(pi*D/param[1])
    return C


def holecosV(D,param):
    """
  
    Translated by Dimitri D'Or - November 2004

    Status: debugged, not tested

    % holecosV                  - cosinusoidal hole effect variogram model (Jan 1,2001)
    %
    % Compute the cosinusoidal hole effect variogram matrix
    % from a distance matrix.
    %
    % SYNTAX :
    %
    % [V]=holecosV(D,param);
    %
    % INPUT :
    %
    % D       n by m   distance matrix, having real values >=0.
    % param   1 by 2   vector of parameters such that param[0] is the
    %                  half amplitude and param[1] is half the periodicity
    %                  of the hole effect.
    %
    % OUTPUT :
    %
    % V       n by m   variogram matrix with same size as D.
    %
    % NOTE :
    %
    % This model is only valid in a 1-dimensional space.
    """

    V=param[0]*(1-cos(pi*D/param[1]))
    return V


def holesinC(D,param):
    """
  
    Translated by Dimitri D'Or - November 2004

    Status: debugged, not tested

    % holesinC                  - sinusoidal hole effect covariance model (Jan 1,2001)
    %
    % Compute the sinusoidal hole effect covariance matrix
    % from a distance matrix.
    %
    % SYNTAX :
    %
    % [C]=holesinC(D,param);
    %
    % INPUT :
    %
    % D       n by m   distance matrix, having real values >=0.
    % param   1 by 2   vector of parameters such that param[0] is
    %                  the sill and param[1] is the periodicity of
    %                  the hole effect.
    %
    % OUTPUT :
    %
    % C       n by m   covariance matrix with same size as D.
    %
    % NOTE :
    %
    % This model is only valid in a d-dimensional space with d<=3.
    """

    isnull=find(D==0)
    put(D,isnull,NaN)
    C=param[0]*sin(pi*D/param[1])/((pi*D)/param[1])
    put(C,isnull,param[0])

    return C

def holesinV(D,param):
    """
  
    Translated by Dimitri D'Or - November 2004

    Status: debugged, not tested

    % holesinV                  - sinusoidal hole effect variogram model (Jan 1,2001)
    %
    % Compute the sinusoidal hole effect variogram matrix
    % from a distance matrix.
    %
    % SYNTAX :
    %
    % [V]=holesinV(D,param);
    %
    % INPUT :
    %
    % D       n by m   distance matrix, having real values >=0.
    % param   1 by 2   vector of parameters such that param[0] is
    %                  the sill and param[1] is the periodicity of
    %                  the hole effect
    %
    % OUTPUT :
    %
    % V       n by m   variogram matrix with same size as D.
    %
    % NOTE :
    %
    % This model is only valid in a d-dimensional space with d<=3.
    """

    isnull=find(D==0)
    put(D,isnull,NaN)
    V=param[0]*(1-sin(pi*D/param[1])/((pi*D)/param[1]))
    put(V,isnull,0)

    return V


def linearV(D,param):
    """
  
    Translated by Dimitri D'Or - November 2004

    Status: debugged, not tested

    % linearV                   - linear variogram model (Jan 1,2001)
    %
    % Compute the linear variogram matrix from a distance matrix.
    %
    % SYNTAX :
    %
    % [V]=linear(D,param);
    %
    % INPUT :
    %
    % D       n by m   distance matrix, having real values >=0.
    % param   scalar   positive slope of the linear model.
    %
    % OUTPUT :
    %
    % V       n by m   variogram matrix with same size as D.
    """

    V=param*D
    return V


def model2indic(D,model,param,p):
    """
  
    Translated by Dimitri D'Or - November 2004

    Status: debugged, not tested

    % model2indic               - indicator covariance model from Gaussian distribution (Jan 1,2001)
    %
    % Compute the values of an indicator variogram or covariance function at
    % specified distances. The indicator variogram or covariance is computed
    % under the assumption of a Gaussian distribution for one or two given
    % value(s) of the cumulative distribution function and for a given
    % theoretical variogram or covariance function model with specified parameters.
    % 
    % SYNTAX :
    %
    % [V]=model2indic(D,model,param,p);
    %
    % INPUT :
    %
    % D       n by m       matrix of distances for which the indicator
    %                      variogram or covariance function model is evaluated.
    % model   string       that contains the name of the variogram or covariance
    %                      model. Only variogram models that have a covariance
    %                      counterpart can be used here. 
    % param   1 by k       vector of parameters for the specified variogram or
    %                      covariance model, where param[0] is the sill of the
    %                      model and is also the variance of the Gaussian distribution.
    % p       1 by 1 or 2  cumulative distribution function value(s) for which the
    %                      indicator variogram or covariance function must be computed.
    %                      If p is a scalar, the function computes the indicator variogram
    %                      or covariance function for the probability of being lower than
    %                      the Gaussian quantile associated with the p value. If p is a
    %                      1 by 2 vector, the function computes the indicator variogram or
    %                      covariance function for the probability of belonging to the
    %                      interval defined by the Gaussian quantiles associated with the
    %                      values in p, where 0<p[0]<p[1]<1.
    %
    % OUTPUT :
    %
    % V      indicator covariance/variogram matrix with same size as D.
    %
    % NOTE :
    %
    % If nested models must be specified, the same conventions as for modelplot.m are used.
    """

    ## %%% Initialize the parameters

    [m,n]=D.shape
    covparam[0]=1
    np=len(p)

    if not isinstance(model,list):
        nm=1
        model=[model]
        param=[param]
    else:
        nm=len(model)
        
    ## %%% Compute the covariance matrix for variogram
    ## %%% or covariance function model with variance
    ## %%% reset to 1

    C=zeros((m,n))
    sill=0
    for i in xrange(nm):
        sill=sill+param[i][0]
        C0=model[i](0,param[i])
        if C0!=0:
            C=C+model[i](D,param[i])
        else:
            C=C+param[i][0]-model[i](D,param[i])
    C=C/sill

    ## %%% Case when p is scalar

    if np==1:

        q=gaussinv(p,[0, 1])
        for i in xrange(m):
            for j in xrange(n):
                if C[i,j]==1:
                    C[i,j]=p-p**2
                else:
                    C[i,j]=gaussbicdf(q,q,[0, 0, 1, 1, C[i,j]])-p**2

    ## %%% Case when p is vector

    if np==2:

        if p[0]>=p[1]:
            return 'p[0] must be strictly lower than p[1]'

        pinf=p[0]
        psup=p[1]
        qinf=gaussinv(pinf,[0, 1])
        qsup=gaussinv(psup,[0, 1])
        p=psup-pinf

        for i in xrange(m):
            for j in xrange(n):
                if C[i,j]==1:
                    C[i,j]=p-p**2
                else:
                    C1=gaussbicdf(qsup,qsup,[0, 0, 1, 1, C[i,j]])
                    C2=gaussbicdf(qsup,qinf,[0, 0, 1, 1, C[i,j]])
                    C3=gaussbicdf(qinf,qsup,[0, 0, 1, 1, C[i,j]])
                    C4=gaussbicdf(qinf,qinf,[0, 0, 1, 1, C[i,j]])
                    C[i,j]=C1-C2-C3+C4-p**2

    ## %%% Set V for variogram or covariance model

    if C0==0:
        V=(p-p**2)-C
    else:
        V=C

    return V



def modelplot(d,model,param,Property=None,Value=None):
    """
    
    Translated by Didrik Pinte - July 2004

    Status: debugged, not tested
    
    % modelplot                 - plot variogram or covariance models (Jan 1,2001)
    %
    % Plot the values of single variogram and covariance models,
    % as well as for nested variogram and covariance models.
    % Nested models are defined as linear combinations of several
    % basic models.
    %
    % SYNTAX :
    %
    % [v]=modelplot(d,model,param,Property,Value)
    %
    % INPUT :
    %
    % d         n by 1   vector of sorted distance values for which the
    %                    variogram or covariance model must be computed
    %                    or displayed.
    % model     string   that contains the name of the variogram/covariance 
    %                    model.
    % param     1 by k   vector of parameters for model, according to the
    %                    conventions for the corresponding variogram or 
    %                    covariance model.
    % Property  1 by p   cell array where each cell cell is a string that contains
    %                    a legal name of a plot object property. This variable is
    %                    optional, as default values are used if Property is missing
    %                    from the input list of variables. Execute get(H), where H is
    %                    a plot handle, to see a list of plot object properties and
    %                    their current values. Execute set(H) to see a list of plot
    %                    object properties and legal property values. See also the help
    %                    for plot.m.
    % Value     1 by p   cell array where each cell is a legal value for the corresponding
    %                    plot object property as specified in Property.
    %
    % OUPUT :
    %
    % v         n by 1   optional vector of estimated variogram or covariance values at
    %                    the distances specified in d.
    %
    % NOTE :
    %
    % 1- For example, when Property=['Color','Linewidth'] and Value=['[0.5 0.5 1]',1],
    % the model will be displayed as a purple broken line with a width of 1 pixel.
    % By default, modelplot.m will use the default properties for plot.m.
    %
    % 2- When the output variable is specified, the graphic is not displayed and
    % modelplot.m simply returns the values for this variable instead. 
    %
    % 3- If nested models must be displayed, the model variable is a cell array
    % such that each cell is a string that contains the name of a model. E.g.,
    % for a nested variogram model including a nugget effect model and an 
    % exponential model, model=['nuggetV','exponentialV']. In that case, the param 
    % variable is a cell array too, where each cell contains the parameters of
    % the corresponding model. E.g., if the nugget effect is equal to 0.2 and the
    % exponential model has a sill equal to 0.8 and a range equal to 1,
    % param=[0.2,[0.8 1]].
    """
    ##%%%%%% Initialize the parameters

    if not isinstance(model,list):
        nm=1
        model=[model]
        param=[param]
    else:
        nm=len(model)


    if not Property==None:
        if not isinstance(Property,list):
            Property=[Property]
            Value=[Value]
            noptions=1
        else:
            noptions=len(Property)
    else:
        noptions=0


    ##%%%%%% Compute the variogram/covariance model values

    v=zeros(d.shape)
    for i in xrange(0,nm):
        out = model[i](d,param[i])
        v = v + out

    ##%%%%%% Display the variogram/covariance model values if required

   ##if nargout==0,
  ##test=(ishold==1)
    
    a=plot(d.flat,v.flat)
    for j in xrange(0,noptions):
        set(a,Property[j],Value[j])
    xlabel('Distance')
    ylabel('Variogram/Covariance')


def modelsyntax():
    """
    
    Translated by Dimitri D'Or - November 2004

    Status: debugged, not tested, not adapted to python grammar
    
    % modelsyntax               - Syntaxical help for using and creating models (Jan 01,2001)
    %
    % modelsyntax provides a help file to explain the syntax of the covariance
    % models in BMELIB works.
    %
    % SYNTAX :
    %
    % help modelsyntax;
    %
    % or, more simply, you may just type:
    %
    % modelsyntax;
    %
    %
    % CONTENT OF HELP FILE:
    %
    % 1- Basic Covariance models
    % --------------------------
    % 
    % A covariance model in BMELIB allows to compute the covariance matrix 
    % based on the distances between two sets of coordinates. The covariance
    % model is speficied by its name covmodel, and its parameters covparam.
    %
    % Example 1.1:
    %  Consider an exponential covariance for a 1-Dimentional domain with
    %  sill cc=1 and range aa=1. A plot of this covariance is created as follow
    %   cc=1;
    %   aa=1;
    %   covmodel='exponentialC';
    %   covparam=[cc aa];
    %   c1=[0];
    %   c2=[0:0.1:2]';
    %   [K]=coord2K(c1,c2,covmodel,covparam);
    %   figure;
    %   plot(c2,K);
    %   xlabel('Spatial lag r');
    %   ylabel('Covariance C(r)');
    %
    % Note that in BMELIB, the covariance is calculated using the function
    %   [K]=coord2K(c1,c2,covmodel,covparam);
    % This function calculates the distance matrix between the set of 
    % coordinates c1 and c2 as follow
    %   [D]=coord2dist(c1,c2);
    % and then it computes the covariance matrix by using covmodel as follow
    %   [K]=eval([covname '(D,covparam)']);
    %
    % Example 1.2:
    %  Consider the same exponential covariance as above. The covariance
    %  matrix can be directly computed from the distance matrix as follow
    %   [D]=coord2dist(c1,c2);
    %   [K]=exponentialC(D,covparam);
    %
    % It is worthwhile emphasizing as seen in these examples that 
    % the covariance model supported by BMELIB are calculated based on 
    % the distance between two sets of points. It is conceivable that the
    % a user may change the covariance model by directly providing her/his
    % own coord2K function. Such a function should then directly calculate
    % the covariance values based on the set of coordinates c1 and c2.
    % However as long as working with distances is acceptable, one should
    % take advantage of the BMELIB built in covariance models. The 
    % covariance models included with BMELIB are the following:
    %
    % nuggetC                   - nugget effect covariance model
    % sphericalC                - spherical covariance model
    % exponentialC              - exponential covariance model
    % gaussianC                 - Gaussian covariance model
    % holecosC                  - cosinusoidal hole effect covariance model
    % holesinC                  - sinusoidal hole effect covariance model
    %
    % For each of the above covariance model the syntax to calculate 
    % a covariance matrix is covname(D,covparam), where covname is the 
    % name of the covariance function, D is the matrix of distances, and
    % covparam is a vector of parameter. This first element of the 
    % parameter vector, i.e. covparam[0], is the sill of the covariance
    % model. In the case of the spherical, exponential and gaussian
    % covariance models, covparam [1] is the range of the covariance, i.e.
    % the distance to reach 95% of the sill value.
    %
    %
    % 2 - space/time covariance models
    % --------------------------------
    %
    % A coordinate set c1 is expressed as a n by nd matrix, where 
    % n is the number of point and nd is the dimension of the spatial
    % domain. For example in the 2D domain, the origin is given by
    %   c1=[0 0];
    % This concept is extended to the space time domain by specifying
    % an additional column which holds the time. For e.g. the origin
    % of a 2D spatial domain at time t=1 is given by
    %   c1=[0 0 1];
    % Whether a point is spatial only or space/time is determined 
    % by the covariance model. In BMELIB the space/time covariance
    % models comes in two flavors: separable, and non-separable.
    % A space/time separable covariance model is coded as 
    % 'covmodelS/covmodelT', where covmodelS refers to the spatial 
    % covariance function and covmodelT refers to the temporal 
    % covariance function. A non-separable covariance model is coded 
    % as 'covmodelST'
    %
    % Example 2.1
    %  Consider an space/time separable covariance model with a sill
    %  of cc=1, a spatial gaussian component with spatial range as=1
    %  in a 2D spatial domain, and a temporal exponential domain with 
    %  a temporal range of at=2. A space/time color plot of this 
    %  covariance is created as follow
    %   cc=1;
    %   as=1;
    %   at=2;
    %   covmodel='gaussianC/exponentialC';
    %   covparam=[cc as at];
    %   c1=[0 0 0];
    %   sg=(0:0.1:1.5)';
    %   nsg=length(sg);
    %   tg=(0:0.1:3)';
    %   ntg=length(tg);
    %   c2=[kron(ones(ntg,1),sg) zeros(nsg*ntg,1) kron(tg,ones(nsg,1))];
    %   [K]=coord2K(c1,c2,covmodel,covparam);
    %   sgMat=reshape(c2(:,1),nsg,ntg);
    %   tgMat=reshape(c2(:,end),nsg,ntg);
    %   Kmat=reshape(K,nsg,ntg);
    %   figure;
    %   pcolor(sgMat,tgMat,Kmat);
    %   colorbar;
    %   shading interp;
    %   title('Covariance C(r,t)');
    %   xlabel('Spatial lag r');
    %   ylabel('Temporal lag t');
    %
    % In this example the '/' character in the name of the covariance
    % indicates that the model is space/time separable. As a result
    % the last row of coordinates is taken as time, while the other rows
    % are for the spatial position. Hence the two sets of space/time 
    % coordinates c1 and c2, it is possible to calculate a matrix of 
    % spatial distances Ds, as well as a matrix of temporal distances Dt.
    % In other words, the calculation of the covariance matrix as follow
    %   [K]=coord2K(c1,c2,'gaussianC/exponentialC',[cc aa at])
    % is equivalent to calculating the spatial and time distance matrices,
    % and then calling the covariance in space and time as follow
    %   [Ds]=coord2dist(c1(:,1:end-1),c2(:,1:end-1));
    %   [Dt]=coord2dist(c1(:,end),c2(:,end));
    %   [K]=cc*gaussianC(Ds,[1 aa]).*exponentialC(Dt,[1 at]);
    %
    % Example 2.2
    %  A space/time non-separable gaussian covariance model with a sill
    %  of cc=1, and a space/time range of ast=1, where the space/time
    %  distance is defined as follow
    %    space/time distance = spatial distance + rst * temporal distance,
    %  where rst=0.5 is a space/time metric, is given by the following 
    %  covariance model:
    %   cc=1;
    %   ast=1;
    %   rst=0.5;
    %   covmodel='gaussianCST';
    %   covparam=[cc ast rst];
    %
    % In this example the 'ST' characters in the name of the covariance
    % indicates that the model is space/time separable. In this case 
    % a spacial distance matrix Ds and a temporal distance matrix Dt
    % are also calculated, but the covariance matrix is obtained as
    %   [K]=gaussianCST(Ds,Dt,covparam);
    %
    %
    % 3 - Nested covariance model
    % --------------------------------
    %
    % When more than one covariance structure has to
    % be specified, covmodel is a cell array where each
    % cell is a string that contains the name of a model.
    % Accordingly, covparam is a cell array where each cell
    % contains the parameters of the corresponding covariance
    % model.
    %
    % Example 3.1
    %  To add a nugget effect with variance 0.2 to the 
    %  space/time separable covariance model described in Example 2.1
    %  we just use
    %   cNugget=0.2;
    %   cc=1;
    %   as=1;
    %   at=2;
    %   covmodel={'nuggetC/nuggetC','gaussianC/exponentialC'};
    %   covparam={ [cNugget], [cc as at] };
    %
    %
    % 4 - Cross-Covariance model for several variables
    % ------------------------------------------------
    %
    % It is possible to specify an additional index for the coordinates
    % c1 and c2 taking integer values from 1 to nv. This index
    % specifies which variable is known at each coordinate.
    % In that case, c1 and c2 are cell arrays, where the first
    % cell is the matrix of coordinates and the second cell is
    % the vector of index values. Each covparam vector is then a
    % cell array where the first cell is a symmetric covariance
    % matrix for that model, and the second cell array is a
    % vector that contains the other parameters of the models.
    % The nhmax and nsmax variables are then vectors of nv
    % elements, where sum(nsmax)<=20.
    %
    % Example 4.1
    %  Consider two space/time variables X(s,t) and Y(s,t) in a
    %  spatial domain of dimension 2, i.e. s=[s1 s2].
    %  Consider the separable covariance model described in Example 3.1,
    %  which has a nugget component nested with a exponentialC/gaussianC
    %  separable component having a spatial range of as=1 and temporal 
    %  range of at=2. Let's assume that the sills (cross-variances) of 
    %  the nugget component for the variables X and Y is given by 
    %    cNugget=[0.4 0.1; ...
    %             0.1 0.1];
    %  and the sill (cross-variances) for the exponentialC/gaussianC
    %  component is given by
    %    cc=[1  0.5; ... 
    %        0.5 2];
    %  A space/time color plot of the space/time cross-covariance 
    %  between X and Y is created as follow
    %   cNugget=[0.4 0.1;0.1 0.1];
    %   cc=[1  0.5; 0.5 2];
    %   as=1;
    %   at=2;
    %   covmodel={'nuggetC/nuggetC','gaussianC/exponentialC'};
    %   covparam={ {cNugget,[]}, {cc, [as at]} };
    %   c1=[0 0 0];
    %   index1=[1];
    %   p1={c1,index1};
    %   sg=(0:0.1:1.5)';
    %   nsg=length(sg);
    %   tg=(0:0.1:3)';
    %   ntg=length(tg);
    %   c2=[kron(ones(ntg,1),sg) zeros(nsg*ntg,1) kron(tg,ones(nsg,1))];
    %   n2=length(c2);
    %   index2=2*ones(n2,1);
    %   p2={c2,index2};
    %   [K]=coord2K(p1,p2,covmodel,covparam);
    %   sgMat=reshape(c2(:,1),nsg,ntg);
    %   tgMat=reshape(c2(:,end),nsg,ntg);
    %   Kmat=reshape(K,nsg,ntg);
    %   figure;
    %   pcolor(sgMat,tgMat,Kmat);
    %   colorbar;
    %   shading interp;
    %   title('Cross covariance C_{X,Y}(r,t)');
    %   xlabel('Spatial lag r');
    %   ylabel('Temporal lag t');
    """

    help(modelsyntax)


def nparammodels():
    """
    
    Translated by Dimitri D'Or - November 2004

    Status: Approved
    
    % nparammodels              - name and # of parameters for covariance models (Jan 1, 2001)
    %
    % Return a cell array that contains the name of the available
    % covariance model and the number of parameters associated with
    % each one of these models.
    %
    % This function is used by coord2K in the space-time case when
    % the space-time covariance model is separable.
    %
    % SYNTAX :
    %
    % [nparam,models]=nparammodels;
    %
    % OUTPUT :
    %
    % models    k by 1  cell array where each cells is a string that contains
    %                   the name of an available covariance model.
    % nparam    k by 1  cell array that contains the number of parameters for
    %                   the corresponding covariance models.
    """

    models=[
      'exponentialC',
      'gaussianC',
      'holecosC',
      'holesinC',
      'nuggetC',
      'sphericalC',
    ]

    nparam=[
       2,
       2,
       2,
       2,
       1,
       2,
    ]

    return models,nparam


def nuggetC(D,param):
    """
    
    Translated by Dimitri D'Or - November 2004

    Status: debugged, not tested
    
    % nuggetC                   - nugget effect covariance model (Jan 1,2001)
    %
    % Compute the nugget effect covariance matrix from a distance matrix.
    % 
    % SYNTAX :
    %
    % [C]=nuggetC(D,param);
    %
    % INPUT :
    %
    % D       n by m   distance matrix, having real values >=0.
    % param   scalar   sill of the model.
    %
    % OUTPUT :
    %
    % C       n by m   covariance matrix with same size as D.
    """

    C=zeros(D.shape)
    if prod(D.shape)!=0:    #not isempty(D):
      C=param[0]*(D==0)

    return C


def nuggetCST(Ds,Dt,param):
    """
    
    Translated by Dimitri D'Or - November 2004

    Status: debugged, not tested
    
    % nuggetCST                 - Nugget non separable S/T covariance model (Jan 1, 2001)
    %
    % Nugget space-time covariance model (i.e. model with nugget
    % effect in space and time)
    %
    % SYNTAX :
    %
    % [C]=nuggetCST(Ds,Dt,param);
    %
    % INPUT :
    %
    % Ds     n by m   distance matrix in space, having real values >=0.
    % Dt     n by m   distance matrix in time with same dimensions as Ds
    % param  1 by 1   vector of parameters such that :
    %                 param[0] is the sill of the space-time model,
    %
    % OUTPUT :
    %
    % C      n by m   space/time covariance matrix with same size as Ds and Dt.
    """

    Cs=nuggetC(Ds,1)
    Ct=nuggetC(Dt,1)
    C=param[0]*Cs*Ct

    return C

def nuggetV(D,param):
    """
    
    Translated by Didrik Pinte - July 2004

    Status: debugged, not tested
    
    % nuggetV                   - nugget effect variogram model (Jan 1,2001)
    %
    % Compute the nugget effect variogram matrix from a distance matrix.
    %
    % SYNTAX :
    %
    % [V]=nuggetV(D,param)
    %
    % INPUT :
    %
    % D       n by m   distance matrix, having real values >=0.
    % param   scalar   sill of the model.
    %
    % OUTPUT :
    %
    % V       n by m   variogram matrix with same size as D.
    """
    if isinstance(param,(list,arraytype)):
        param = param[0]
    if not isinstance(param,(float,int,long)):
        print str(type(param))
        raise ValueError, "Element must be scalar. Parameter is : " + str(param)
    D = asarray(D)
    V = where(D!=0,param,0)
    return V


def power(D,param):
    """
    
    Translated by Dimitri D'Or - November 2004

    Status: debugged, not tested
    
    % powerV                    - power variogram model (Jan 1,2001)
    %
    % Compute the power variogram matrix from a distance matrix.
    %
    % SYNTAX :
    %
    % [V]=powerV(D,param);
    %
    % INPUT :
    %
    % D       n by m   distance matrix, having real values >=0.
    % param   1 by 2   vector of parameters such that :
    %                  param[0] is a slope factor, with 0<param[0]<2.
    %                  param[1] is the exponent for the power model.
    %
    % OUTPUT :
    %
    % V       n by m   variogram matrix with same size as D.
    """

    V=param[0]*(D**param[1])

    return V

def sphericalC(D,param):
    """
    
    Translated by Dimitri D'Or - November 2004

    Status: debugged, not tested
    
    % sphericalC                - spherical covariance model (Jan 1,2001)
    %
    % Compute the spherical covariance matrix from a distance matrix.
    % 
    % SYNTAX :
    %
    % [C]=sphericalC(D,param);
    %
    % INPUT :
    %
    % D       n by m   distance matrix, having real values >=0.
    % param   1 by 2   vector of parameters such that :
    %                  param[0] is the sill of the model,
    %                  param[1] is the range of the model.
    %
    % OUTPUT :
    %
    % C       n by m   covariance matrix with same size as D.
    %
    % NOTE :
    %
    % This model is only valid in a d-dimensional space with d<=3.
    """

    C=param[0]*(1-((3/2)*(D/param[1])-(1/2)*(D/param[1])**3))
    index=find(D>param[1])
    put(C,index,0)

    return C


def sphericalV(D,param):
    """
    
    Translated by Dimitri D'Or - November 2004

    Status: debugged, not tested
    
    % sphericalV                - spherical variogram model (Jan 1,2001)
    %
    % Compute the spherical variogram matrix from a distance matrix.
    % 
    % SYNTAX :
    %
    % [V]=sphericalV(D,param);
    %
    % INPUT :
    %
    % D       n by m   distance matrix, having real values >=0.
    % param   1 by 2   vector of parameters such that :
    %                  param[0] is the sill of the model,
    %                  param[1] is the range of the model.
    %
    % OUTPUT :
    %
    % V       n by m   variogram matrix with same size as D.
    %
    % NOTE :
    %
    % This model is only valid in a d-dimensional space with d<=3.
    """

    V=param[0]*((3/2)*(D/param[1])-(1/2)*(D/param[1])**3)
    index=find(D>param[1])
    put(V,index,param[0])

    return V
