from scipy import *
import scipy
# from matplotlib.matlab import *

def automesh(*args):
    '''
    
    Translated by Dimitri D''Or - November 2004

    Status: completed, not tested for x,y,z case
        
    This function has the same effect as the MATLAB function AUTOMESH:
    
    % AUTOMESH True if the inputs should be automatically meshgridded.
    %    AUTOMESH(X,Y) returns true if X and Y are vectors of
    %    different orientations.
    %
    %    AUTOMESH(X,Y,Z) returns true if X,Y,Z are vectors of
    %    different orientations.
    %
    %    AUTOMESH(...) returns true if all the inputs are vectors of
    %    different orientations.

    %   Copyright 1984-2002 The MathWorks, Inc.
    %    $Revision: 1.7 $ $Date: 2004/11/24 16:07:26 $
    '''

    for i in xrange(0,len(args)):
        if len(args[i].shape)==1:
            args[i]=args[i][:,NewAxis]
        
    ns=[]
    isvec=[]
    nd=[]
    
    for i in xrange(0,len(args)):
        ns.append(map(not_equal,args[i].shape,ones(len(args[i].shape))))    #% Location of non-singleton dimensions
        isvec.append(sum(ns[i])<=1)                                        #% Is vector.
        nd.append(len(args[i].shape))                                      #% Number of dimensions.

    nd=asarray(nd)

    ##% True if inputs are 2-D, all vectors, and their non-singleton
    ##% dimensions aren't along the same dimension.
    if alltrue(nd==2) & alltrue(isvec) & (not isequal(ns)):
        return True
    else:
        return False


def clearall():
    '''
    
    Translated by Dimitri D''Or - November 2004

    Status: not working
        
    This function has the same effect as the MATLAB function clear(all):
    It clears all the variables created by the user in the namespace.
    '''

    [vars().__delitem__(_k) for _k in vars().keys() if not _k.startswith('_')]

    

def kron(a,b):
    '''
    Thanks to Nils Wagner Wed Sep 24 16:52:28 CDT 2003 (see the Scipy mailing list)
    
    Imported by Didrik Pinte - September 2004

    Status: approved
    '''
    
    if not a.iscontiguous():
        a = reshape(a, a.shape)
    if not b.iscontiguous():
        b = reshape(b, b.shape)
    o = outerproduct(a,b)
    o.shape = a.shape + b.shape
    return concatenate(concatenate(o, axis=1), axis=1)

def isequal(vals):
    ''' 
    Written under the name allTheSame by Kent Johnson (on tutor@python.org) - 19 November 2004

    Status: Approved

    Test if vals is an iterable whose elements are all equal.
    Returns True if all the elements in tuple (or list or array or string) vals
    are equal and False, else.

    (Idem as Matlab function ISEQUAL)
    '''

    import itertools

    i = iter(vals)  # Raises TypeError if vals is not a sequence

    try:
        first = i.next()
    except StopIteration:
        # vals is an empty sequence
        return True

    for item in i:
        if first != item:
            return False

    return True


def ismember(a,s,flag=None):
    ''' 
    Written by Dimitri D''Or- 19 November 2004

    Status: debugged, not fully tested

    Test for each element in a if it is a member of s.
    Returns a list of booleans with True if a[i] is a member of s and False, else.

    (Idem as Matlab function ISMEMBER)
    '''

    # check the sizes and types of the elements
    
    if type(a[0])!=type(s[0]):
        return 'Elements to compare must be of the same type'
    if (not isinstance(a[0],int)) and (not isinstance(a[0],float)):
        if len(a[0])!=len(s[0]):
            return 'Elements to compare must have the same length'
        if isinstance(a[0],arraytype) and a[0].shape!=s[0].shape:
            return 'arrays to compare must have the same shape'

    # check for each element of a if it is a member of s

    if flag==None and isinstance(a[0],arraytype):
        b=ismember(ravel(a),ravel(s))
        b=reshape(b,a.shape)
    else:
        b=ones(len(a),)*1000
        for i in xrange(len(a)):
            c=ones(len(s),)*1000
            for j in xrange(len(s)):
                c[j]= (isequal((a[i],s[j])))
                b[i]= (sum(c)>=1)

    return b


def issorted(a,flag=None):
    '''
    
    Translated by Dimitri D''Or - November 2004

    Status: approved
        
    This function has the same effect as the MATLAB function ISSORTED:
    
    % ISSORTED True for sorted vector.
    %   ISSORTED(A) when A is a vector returns 1 if the elements of A are
    %   in sorted order (in other words, if A and SORT(A) are identical)
    %   and 0 if not.
    %
    %   For character arrays, ASCII order is used.
    %
    %   ISSORTED(A,'rows') when A is a matrix returns 1 if the rows of A 
    %   are in sorted order (if A and SORTROWS(A) are identical) and 0 if not.
    %
    '''

    isrows = flag == 'rows'

    tf = True

    if flag == None:

        numelA = len(a)

        if prod(a.shape) != numelA:
            return 'A must be a vector or ''rows'' must be specified.'

        ##% Convert to double and columns
        a = ravel(a).astype('d')

        ##% FOR loop breaks if it encounters something out of order
        for i in xrange(0,numelA-1):
            if not (a[i] <= a[i+1]):
                tf = False
                break

        ##% Check for NaN''s; sorted if everything after breakpoint is NaN
        if (not tf) and a[i] > a[i+1]:
            return tf

        ##% Proceed with this only if NaN found in position (i+1).
        elif isnan(a[-1]):                    ##% Check final element for NaN
            if alltrue(isnan(a[i+1:-1])):  ##% Check all others for NaN
                tf = True
    else:
        ##%% 'rows' implementation
        if not isrows:
            return 'Unknown flag.'

        nrows = a.shape[0]

        for c in xrange(0,nrows-1):
            diffvec = a[c+1,:] - a[c,:];
            for i in xrange(0,len(diffvec)):
                if diffvec[i] > 0:
                    break
                elif diffvec[i] < 0:
                    tf = False
                    return tf
                elif isnan(diffvec[i]):
                    if isnan(a[c+1,i]) and not isnan(a[c,i]):
                        break
                    elif isnan(a[c,i]) and  not isnan(a[c+1,i]):
                        tf = False
                        return tf

    return tf


def meshgrid(x,y,z=None):
    """

    Status: approved
    
    Modified from the matplotlib.mlab version by Dimitri D'Or on 05/10/2004
    to take into account the 3D case
    
    For vectors x, y with lengths Nx=len(x) and Ny=len(y), return X, Y
    where X and Y are (Ny, Nx) shaped arrays with the elements of x
    and y repeated to fill the matrix

    EG,

      [X, Y] = meshgrid([1,2,3], [4,5,6,7])

       X =
         1   2   3
         1   2   3
         1   2   3
         1   2   3


       Y =
         4   4   4
         5   5   5
         6   6   6
         7   7   7
  """

    # 2D Case

    if z==None :
        x = array(x)
        y = array(y)
        numRows, numCols = len(y), len(x)  # yes, reversed
        x.shape = 1, numCols
        X = repeat(x, numRows)

        y.shape = numRows,1
        Y = repeat(y, numCols, 1)
        return X, Y

    # 3D case

    else:
        nx=prod(x.shape)
        ny=prod(y.shape)
        nz=prod(z.shape)
        xx=reshape(x,(nx,1,1))
        yy=reshape(y,(1,ny,1))
        zz=reshape(z,(1,1,nz))
        X=repeat(xx,nx,1)
        X=repeat(X,nz,2)
        Y=repeat(yy,ny,0)
        Y=repeat(Y,nz,2)
        Z=repeat(zz,ny,0)
        Z=repeat(Z,nx,1)
        return X, Y, Z

def testMeshGrid():
    print "2D Case"
    x=array([0,1,2])
    y=x
    print meshgrid(x,y)
    print "3D Case"
    x=array([0,1])
    y=x
    z=x
    print meshgrid(x,y,z)    
    

def  setdiff(a,b,flag=None):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: Approved, not tested for strings
    
    setdiff - set difference (Oct 2004)
    
    This function has the same effect as the MATLAB function SETDIFF:
    
    SETDIFF Set difference.
       SETDIFF(A,B) when A and B are vectors returns the values
       in A that are not in B.  The result will be sorted.
    
       SETDIFF(A,B,'rows') when A are B are matrices with the same
       number of columns returns the rows from A that are not in B.
    
       [C,I] = SETDIFF(...) also returns an index vector I such that
       C = A(I) (or C = A(I,:)).
    
    '''

    isrows = (flag=='rows')

    if len(a.shape)==1:      ## a is then a (x,) array (= vector)
        a = a[:,NewAxis]                 ## a is converted into a (x,1) array
    if len(b.shape)==1:      ## a is then a (x,) array (= vector)
        b = b[:,NewAxis]                 ## a is converted into a (x,1) array

    rowsA = a.shape[0]
    colsA = a.shape[1]
    rowsB = b.shape[0]
    colsB = b.shape[1]

    cond1 = (rowsA > 1 and colsB <= 1)
    cond2 = (rowsB > 1 and colsA <= 1)

    rowvec = not(cond1 or cond2 or isrows)

    if flag==None:
        numelA = len(a)
        numelB = len(b)
        
        if (prod(a.shape) != numelA) or (prod(b.shape) != numelB) :
            return 'A and B must be vectors or ''rows'' must be specified.'
        
        ##% Handle empty arrays.
        
        if (numelA == 0):
        ##% Predefine outputs to be of the correct type.
            c = array([])
            ia = array([])
        ##% Ambiguous if no way to determine whether to return a row or column.
            ambiguous = (rowsA==0 and colsA==0) and ((rowsB==0 and colsB==0) or numelB == 1)
            if ambiguous!=True:
                c = c[:,NewAxis]
                ia =  ia[:,NewAxis]
        elif (numelB == 0):
        ##% If B is empty, invoke BMEmatlab UNIQUE to remove duplicates from A.
            [c,ia,dumb] = unique(a)
            return c, ia
        
        ##% Handle scalar: one element.  Scalar A done only.  
        ##% Scalar B handled within ISMEMBER and general implementation.
        
        elif (numelA == 1):
            if ismember(a,b,flag)!=True:
                c = a
                ia = 1
            else:
                c = array([])
                ia = array([])
            return c,ia
        
        ##% General handling.
        
        else:
            
        ##% Convert to columns.
            a = ravel(a)
            a = a[:,NewAxis]
            b = ravel(b)
            b = b[:,NewAxis]
            
        ##% Convert to double arrays, which sort faster than other types.
            
            typec=a.typecode()
            isdouble=(typec=='d')
            
            if a.typecode()!='d':
                a=a.astype('d')
            
            if b.typecode()!='d':
                b=b.astype('d')
                
        ##% Call ISMEMBER to determine list of non-matching elements of A.
            tf = (ismember(a,b)!=True)
            c = take(a,find(ravel(tf)))
                     
        ##% Call BMEmatlab UNIQUE to remove duplicates from list of non-matches.
            [c,ndx,dumb] = unique(c)
            
        ##% Find indices by using TF and NDX.
            wher = find(ravel(tf))
            if isinstance(ndx,int):
                ndx=[ndx]
            ia = take(wher,ndx)
            
        ##% Re-convert to correct output data type using FEVAL.
            if isdouble!=True:
                c = c.astype(typec)
                
        ##% If row vector, return as row vector.
        if rowvec:
            c = ravel(c)
            ia = ravl(ia)
                         
    else:         ##% 'rows' case
        if isrows!=True:
            return 'Unknown flag.'

##       ##% Automatically pad strings with spaces
##       if logical_and(isstr(a),isstr(b)):
##         if colsA > colsB:
##           b = [b repmat(' ',rowsB,colsA-colsB)]
##         elif colsA < colsB :
##           a = [a repmat(' ',rowsA,colsB-colsA)]
##           colsA = colsB;
##         end
##       elseif colsA ~= colsB
##         error('A and B must have the same number of columns.');
##       end

        ##% Handle empty arrays
        if rowsA == 0:
            c = zeros((rowsA,colsA),'Float32')
            ia = array([])
        elif (colsA == 0) and (rowsA > 0):
            c = zeros((1,0),'Float32')
            ia = rowsA
        ##% General handling
        else:
          ##% Remove duplicates from A; get indices only if needed
            [a,ia,dumb] = unique(a,flag)

          ##% Create sorted list of unique A and B; want non-matching entries
            [c,ndx] = sortrows(concatenate((a,b)))
            [rowsC,colsC] = c.shape
            if logical_and(rowsC > 1,colsC != 0):
          ##% d indicates the location of non-matching entries
                d = sum((c[:rowsC-1,:] != c[1:,:]),1)
            else:
                d = zeros((rowsC-1,0),'Float32')
            d = (d>=1) 
            d = concatenate((d,array([1])))                      #% Final entry always included.

          ##% d = 1 now for any unmatched entry of A or of B.
            n = a.shape[0]
            d = logical_and(d,(ndx < n))           #% Now find only the ones in A.
            c = take(c,find(d==1))
            ia = take(ia,take(ndx,find(d==1)))

##     % Automatically deblank strings
##     if isstr(a)
##       c = deblank(c);
##     end

    return c,ia

def sortrows(x,col=None):
    '''
    
    Written by Dimitri D''Or - October 2004

    Status: approved
    
    sortrows - Multiple sort of the rows of an array using as keys the col vector (Oct 2004)
    
    This function has the same effect as the MATLAB function SORTROWS:

    Y = SORTROWS(X) sorts the rows of the matrix X in ascending order as a
    group.  X is a 2-D numeric or char matrix.  For a char matrix containing
    strings in each row, this is the familiar dictionary sort.  When X is
    complex, the elements are sorted by ABS(X). Complex matches are further
    sorted by ANGLE(X).  X can be any numeric or char class.  Y is the same
    size and class as X.
   
    SORTROWS(X,COL) sorts the matrix based on the columns specified in the
    vector COL.  For example, SORTROWS(X,[2 3]) sorts the rows of X by the
    second and third columns of X.

    [Y,I] = SORTROWS(X) also returns an index matrix I such that Y = X(I,:).

    SYNTAX: [Y,ndx]=sortrows(X,col)

    INPUT:
    
    X   nxm   array of values to sort
    col ncx1  array of integers for the column to be used as keys. Must contain
              integers between 0 and m.
    
    OUTPUT:
    
    Y   nxm   array of sorted values
    NDX nx1   vector.  Double. Contains row indices into X.

    '''
    if col==None :
        col=arrayrange(0,x.shape[1])

    if isinstance(col,int):
        col=array([col])

    if isinstance(col,list):
        col=asarray(col)

    nrows=x.shape[0]

    if len(col)==1:
        ndx=argsort(x[:,col[0]])
        y=take(x,ndx)    # Sorting against the first column
        return y, ndx
    else :
        ndx=arange(0,nrows)[:,NewAxis]
        ndx = ndx.astype(x.typecode())
        x=concatenate((x,ndx),1)
        [y,grouplist]=sortgroups(x,col[0])
        for j in xrange(0,len(grouplist)):
            y_sub=concatenate((y[min(grouplist[j]):max(grouplist[j])+1,0:col[0]],y[min(grouplist[j]):max(grouplist[j])+1,col[0]+1:]),1)
            colsub=where(col>col[0],col-1,col)
            (y_sub,dumb)=sortrows(y_sub,colsub[1:])
            mini=(min(grouplist[j])-min(grouplist[j]))
            maxi=(max(grouplist[j])-min(grouplist[j])+1)
            y[min(grouplist[j]):max(grouplist[j])+1,0:col[0]]=y_sub[mini:maxi,0:col[0]]
            y[min(grouplist[j]):max(grouplist[j])+1,col[0]+1:]=y_sub[mini:maxi,col[0]:]
        ndx=y[:,-1].astype(int)
        y=y[:,:-1]
        
    return y, ndx

def sortgroups(x,col):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: approved
    
    Sort an array against integer col and gives the groups of rows for which the values on fist column are equal
    Subfunction of SORTROWS
    
    '''

    nrows=x.shape[0]

    y=take(x,argsort(x[:,col]))    # Sorting against the column col

    # Searching for groups for which the value on column col is the same
    j=0
    grouplist=[]
    group=array([j])
    while j<nrows-1 :
        if y[j,col]==y[j+1,col]:
            group=concatenate((group,array([j+1])))
        else :
            grouplist.append(group)
            del(group)
            group=array([j+1])
        j=j+1
    grouplist.append(group)

    return y,grouplist

def unique(a,flag=None):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: approved
    
    unique - set unique (Oct 2004)
    
    This function has the same effect as the MATLAB function UNIQUE:

    UNIQUE Set unique.
       UNIQUE(A) for the array A returns the same values as in A but
       with no repetitions.  A will also be sorted.  
    
       UNIQUE(A,'rows') for the matrix A returns the unique rows of A.
    
       [B,I,J] = UNIQUE(...) also returns index vectors I and J such
       that B = A(I) and A = B(J) (or B = A(I,:) and A = B(J,:)).
       
    '''

    rows = a.shape[0]
    rowvec = 0
    if len(a.shape)==1:      ## a is then a (x,) array (= vector)
        a = a[:,NewAxis]
        rowvec = 1
    cols = a.shape[1]

    numelA = prod(a.shape)

    if flag==None:

      ##% Handle empty: no elements.

      if (numelA == 0):
        ##% Predefine b to be of the correct type.
        b = array([])
        if max(a.shape) > 0:
          b = reshape(b,(0,1))
          ndx = zeros((0,1))
          pos = zeros((0,1))
        else:
          ndx = array([])
          pos = array([])
        return b,ndx,pos

      elif (numelA == 1):
        ##% Scalar A: return the existing value of A.
        b = a
        ndx = 0
        pos = 0
        return b,ndx,pos

        ##% General handling.  
      else:

        ##% Convert to columns
        a = ravel(a)

        ##% Convert to double array for purposes of faster sorting.
        ##% Additionally, UNIQUE calls DIFF, which requires double input.

        whichclass = a.typecode()    
        isdouble = (whichclass=='d')

        if isdouble!=True:
          a = a.astype('d')

        ##% Sort if unsorted.  Only check this for long lists.

        checksortcut = 1000

        if (numelA <= checksortcut) or (issorted(a)!=True):
            b = sort(a)
            ndx = argsort(a)
        else:
          b = a
          ndx = reshape(arrayrange(0,numelA),(numelA,1))  ##% If presorted, indices are 1,2,3,...

        ##% d indicates the location of non-matching entries.

        db = diff(b)

        ##% Since DIFF returns NaN in both Inf and NaN cases, 
        ##% use slower method of detection if NaN\'s detected in DIFF(b).
        ##% After sort, Infs or NaNs will be at ends of list only.

        if (isnan(db[0]) or isnan(db[-1])):
          d = b[0:numelA-1] != b[1:numelA]
        else:
          d = (db != 0)

        d = concatenate((d,[1]))                      ##% Final element is always member of unique list.

        index=matplotlib.matlab.find(d==1)
        b = take(b,index)                             ##% Create unique list by indexing into sorted list.

        pos = cumsum(concatenate(([0],d)))            ##% Lists position, starting at 0. ## REMARK D D''OR Nov 2004: in original matlab code d is full(d)
        pos = pos[:numelA]                            ##% Remove extra element introduced by d.
        pos = take(pos,ndx)                           ##% Re-reference POS to indexing of SORT.

        ##% Create indices if needed.
        ndx = take(ndx,index)

        ##% Re-convert to correct output data type using FEVAL.
        if isdouble!=True:
          b = b.astype(whichclass)

      ##% If row vector, return as row vector.
      if rowvec==True:
        b = reshape(b,(len(b),))
        ndx = reshape(ndx,(len(ndx),))
        pos = reshape(pos,(len(pos),))

    else:    ##% 'rows' case
        if flag!='rows':
            return 'Unknown flag.'

        ##% Handle empty: no rows.

        if (rows == 0):
            ##% Predefine b to be of the correct type.
            b = array([])
            ndx = []
            pos = []
            b = reshape(b,(0,cols))
            if cols > 0:
                ndx = reshape(ndx,(0,1))
            return b,ndx,pos

            ##% Handle scalar: one row.

        elif (rows == 1):
            b = a
            ndx = 1
            pos = 1
            return b,ndx,pos

        ##% General handling.
        ##% Conversion to double not done: SORTROWS is slower for doubles
        ##% than other types.

        (b,ndx) = sortrows(a)

        ##% d indicates the location of non-matching entries.

        d = b[0:rows-1,:]!=b[1:rows,:]

        ##% d = 1 if differences between rows.  d = 0 if the rows are equal.

        d = sometrue(d,1)
        d = concatenate((d,[1]))                ##% Final element is always member of unique list.

        index=matplotlib.matlab.find(d==1)
        b = take(b,index)                      ##% Create unique list by indexing into sorted list.

        ##% Create position mapping vector using CUMSUM.

        pos = cumsum(concatenate(([0],d)))     ##% Lists position, starting at 1. ## REMARK D D''OR Nov 2004: in original matlab code d is full(d)
        pos = pos[:rows]                     ##% Remove extra element introduced by d.
        pos = take(pos,ndx)                    ##% Re-reference POS to indexing of SORT.

        ##% Create indices if needed.
        ndx = take(ndx,index)

    return b,ndx,pos



def xyzchk(*args):
    '''
    
    Translated by Dimitri D''Or - October 2004

    Status: Approved
    
    XYZCHK Check arguments to 3-D data routines.
        
    This function has the same effect as the MATLAB function XYZCHK:

      [MSG,X,Y,Z,C] = XYZCHK(Z), or
      [MSG,X,Y,Z,C] = XYZCHK(Z,C), or
      [MSG,X,Y,Z,C] = XYZCHK(X,Y,Z), or
      [MSG,X,Y,Z,C] = XYZCHK(X,Y,Z,C), or
      [MSG,X,Y,Z,XI,YI] = XYZCHK(X,Y,Z,XI,YI) checks the input aguments
      and returns either an error message in MSG or valid X,Y,Z (and
      XI,YI) data.
    
    Dimitri D\'Or 20 Oct 2004
    '''

    def isvector(x):
        #%ISVECTOR True if x has only one non-singleton dimension.
        tf = (len(x) == prod(x.shape))
        return tf

    if len(args) not in xrange(1,6):
        return 'Number of arguments must be between 1 and 6'

    msg = []
    out5 = []
    out6 = []
    
    if len(args)==1:  # % xyzchk(z)
        z = args[0]
        if isinstance(z,str):
            msg = 'Input arguments must be numeric.'
            return msg,[],[],[],[]
        [m,n] = z.shape
        [x,y] = meshgrid(arrayrange(0,n),arrayrange(0,m))
        out5 = z      # % Default color matrix
        return msg, x, y, z, out5

    elif len(args)==2: # % xyzchk(z,c)
        z = args[0]
        c = args[1]
        [m,n] = z.shape
        [x,y] = meshgrid(arrayrange(0,n),arrayrange(0,m))
        if z.shape != c.shape:
            msg = 'Matrix C must be the same size as Z.';
            return msg,[],[],[],[]
        out5 = c
        return msg, x, y, z, out5

    elif len(args)>=3:  # % xyzchk(x,y,z,...)
        x = args[0]
        y = args[1]
        z = args[2]
        if len(z.shape)>1:
            [m,n] = z.shape
        elif len(z.shape)==1:
            m=z.shape[0]
            n=1
        if not isvector(z): # % z is a matrix
            #% Convert x,y to row and column matrices if necessary.
            if isvector(x) & isvector(y):
                [x,y] = meshgrid(x,y)
                if x.shape[1]!=n and y.shape[0]!=m:
                    msg = 'The lengths of X and Y must match the size of Z.'
                    return msg,[],[],[],[]
                elif x.shape[1]!=n:
                    msg = 'The length of X must match the number of columns of Z.'
                    return msg,[],[],[],[]
                elif y.shape[0]!=m:
                    msg = 'The length of Y must match the number of rows of Z.'
                    return msg,[],[],[],[]
            elif isvector(x) or isvector(y):
                msg = 'X and Y must both be vectors or both be matrices.'
                return msg,[],[],[],[]
            else:
                if x.shape != z.shape or y.shape != z.shape:
                    msg = 'Matrices X and Y must be the same size as Z.'
                    return msg,[],[],[],[]
        else: #% z is a vector
            if (not isvector(x)) or (not isvector(y)):
                msg = 'When Z is a vector, X and Y must also be vectors.';
                return msg,[],[],[],[]
            elif (len(x)!=len(z) or len(y)!=len(z)) and ( not ((len(x)==n) and (len(y)==m))):
                msg = 'X and Y must be same length as Z or the lengths of X and Y must match the size of Z.'
                return msg,[],[],[],[]

    if len(args)==4: #% xyzchk(x,y,z,c)
        c = args[3]
        if z.shape != c.shape:
            msg = 'Matrix C must be the same size as Z.';
            return msg,[],[],[],[]
        out5 = c
        return msg, x, y, z, out5

    if len(args)==5: #% xyzchk(x,y,z,xi,yi)
        xi = args[3]
        yi = args[4]

        if automesh(xi,yi):
            [xi,yi] = meshgrid(xi,yi)
        elif xi.shape != yi.shape:
            msg = 'XI and YI must be the same size or vectors of different orientations.'
            return msg,[],[],[],[]
        out5 = xi
        out6 = yi

    return msg,x,y,z,out5,out6
