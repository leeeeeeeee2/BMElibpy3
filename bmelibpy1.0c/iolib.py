from scipy import *

def readGeoEAS(datafile, colnum=[]):
    '''
    Translated by Didrik Pinte - July 2004

    Status: approved

    % readGeoEAS                - Reads data files in GEO format (Jan 1, 2001)
    %
    % Reads specified columns of a data files in GEO format.
    %
    % SYNTAX :
    %
    % [val,valname,filetitle]=readGeoEAS(datafile,colnum);
    %
    % INPUT :
    %
    % datafile   string   data file in Geo EAS format. This format has a title in the
    %                     the first line, the number of column nc in the second line
    %                     followed by nc lines with the name of each column,
    %                     followed by the data lines, each data lines having nc value
    % colnum     vector   optional parameter, specifying the columns to read. 
    %                     max(colnum) cannot be greater than nc. Default value is 1:nc
    %
    % OUTPUT :
    %
    % val        n by nv  matrix of value read for the nv variables, where nv=sum(colnum>0)
    % valname    cells    array of nv cells with the name of the nv variables (for each colum
    %                     of val)
    % filetitle  string   char array with the title of the file
    %
    % EXAMPLE :
    % 
    % Following is an example of a valid file in GeoAES format:
    %
    %BeginningOfFile--------------------------------------------------
    %BME hard data
    %3
    %s1
    %s2
    %primary variable
    %10.1    11.4    1.4 
    %12.2    13.6    1.6 
    %16.7    19.1    1.1 
    %10.9    16.9    0.9 
    %EndOfFile--------------------------------------------------------
    '''

    ## Opens the datafile read-only
    f=open(datafile, 'r')
    try:       
        ## Reads the file title
        filetitle = f.readline().rstrip()
        ## Reads the number of column
        nc = int(f.readline())
        ## Checks if nc exist
        if not nc:
            print "Problem reading the number of variables on line 2 of " + filetitle
            print nc
            return
        
        # Checks if the parameters colnum is there. If not, it creates a vector [0..nc]
        if len(colnum)<1:
            colnum = range(int(nc));
            
            valname = []
            val=[]    
            
            ## Test if the colnum index are all greater than 0
            if sum(where(colnum<0,1,0))>1:
                print('All the column number specified in colnum must be greater than 0')
                f.close()
                return(val,valname,filtetitle)
            
           ## Test if length of colnum is equal to 0. If true, return empty parameters
            if len(colnum)==0:
                return (val, valname, filetitle)
            
            ## Test if nc if greater than max(colnum)
            if nc < max(colnum):
                print datafile + ' has a request for column '
                print max(colnum)
                print ', but it has only '
                print nc
                print ' columns'
                return (val,valname,filetitle)

            ## Parses the column names in the file        
            datafilevalname = [];
            for ic in xrange(nc):               ##% Read the name of the data file columns
                colName= f.readline().rstrip()
                datafilevalname.append(colName)

            ## Assigns the read names to the valname list
            for i in xrange(len(colnum)):      ##% Assign the name of the specified columns to valname
                ##print datafilevalname[colnum[i]]
                valname.append(datafilevalname[colnum[i]])

            ## Parsing the content
            linenumber=nc+2
            for line in f.readlines():
                ## Splits the data using the whitespaces (uses strip() first to remove the starting and leading whitespaces)                       
                tempval = line.strip().split()         
                linenumber +=1
                ## Test if line has enough values to read
                if len(tempval)!=nc:
                    print 'Line ' + str(linenumber) + ' of file '+ datafile + ' should have '+ str(nc) + ' values'
                    print tempval
                    return (val,valname,filetitle)
                ## Appending the output readen list
                val.append(createOutputList(tempval, colnum))         
                        
            ## Returns the parameters converting the val list to a Numeric.array
            print 'Data loaded'
            return (array(val,Float32),valname, filetitle)
    finally:
        f.close()
                        
    
def createOutputList(mylist, columns):
    '''
    Written by Didrik Pinte - July 2004
    Subfunction for function readGeoEAS
    '''
    output = [float(mylist[columns[i]])  for i in range(len(columns))]
    return output


def writeGeoEAS(val,valname,filetitle,datafile,colnum=None,outofrangeval=0):
    '''
    Translated by Dimitri D''or - November 2004

    Status: not completed

    % writeGeoEAS               - Writes data files in Geo EAS format (Jan 1, 2001)
    %
    % Writes data in specified columns of a Geo EAS data files.
    %
    % SYNTAX :
    %
    % writeGeoEAS(val,valname,filetitle,datafile,colnum,outofrangeval);
    %
    % INPUT :
    %
    % val           n by nv    matrix with the values to write.
    % valname       cell       1 by nv cell array with the name of the 
    %                          variables (one for each colum of val)
    % filetitle     string     title of the file
    % datafile      string     data file in Geo EAS format. This format has a title in the
    %                          the first line, the number of column nc in the second line
    %                          followed by nc lines with the name of each column,
    %                          followed by the data lines, each data lines having nc value
    % colnum        vector     optional paramater specifying the column assigned to each variable.
    %                          length(colnum) must be equal to the number of column of val.
    %                          max(colnum) is the number of column nc of the datafile
    %                          If columns of datafile that are not assigne to a variable
    %                          are filled with outofrangeval.
    %                          Default value is colnum=1:size(val,2)
    % outofrangeval scalar     option value to use when writing NaNs. (Default=0)
    %
    % NOTE :
    %
    % See also readGeoEAS
    '''
    
    ## if nargin<5, colnum=1:size(val,2); end
    ## if nargin<6, outofrangeval=0; end

    ## if sum(colnum<=0)>1  
    ##   error('All the column number specified in colnum must be greater than 0');
    ## end
    ## if length(colnum)~=size(val,2)
    ##   error('length(colnum) must be equal to the number of column of val');
    ## end;

    ## nc=max(colnum);
    ## nvar=size(val,2);

    ## fid=fopen(datafile,'w');       % Open the data file
    ## if fid==-1,
    ##    fclose('all');
    ##    error(sprintf('Problem opening file %s',datafile));
    ## end;

    ## fprintf(fid,'%s\r',filetitle);  % Write the title of the data file
    ## fprintf(fid,'%d\r',nc);         % Write the number of columns
    ## for ic=1:nc,                    % Write the name of the data file columns
    ##   ivar=find(colnum==ic);
    ##   if isempty(ivar)
    ##     fprintf(fid,'%s\r','column not used');
    ##   else   
    ##     fprintf(fid,'%s\r',valname{ivar});
    ##   end
    ## end;

    ## nval=size(val,1);               % Reassign val to valfile
    ## valfile=outofrangeval*ones(nval,nc);
    ## valfile(1:nval,colnum)=val;
    ## valfile(isnan(valfile))=outofrangeval;

    ## for i=1:nval
    ##   str=sprintf('%12.10g ',valfile(i,:)); 
    ##   fprintf(fid,'%s\r',str);      % Write each data line
    ## end

    ## st=fclose(fid);
    ## if st==-1,
    ##   warning(sprintf('could not close %s properly',datafile));
    ## end;

    if colnum==None :
        xrange(1,val.shape[1])
    if sum(colnum<=0)>1 : 
      return 'All the column number specified in colnum must be greater than 0'
    if len(colnum)!=val.shape[1]:
      return'length(colnum) must be equal to the number of column of val'
  
    nc=max(colnum)
    nvar=val.shape[1]

    fid=open(datafile,'w')          #% Open the data file
    ## if fid==-1,
    ##    fclose('all');
    ##    error(sprintf('Problem opening file %s',datafile));
    ## end;

    fid.write(filetitle,'\r')            #% Write the title of the data file
    fid.write(nc,'\r')                   #% Write the number of columns
    for ic in xrange(0,nc) :             #% Write the name of the data file columns
        ivar=find(colnum==ic)
        if len(ivar)==0:
            fid.write('column not used\r')
        else :
            fid.write(valname[ivar],'\r')

    nval=val.shape[0]                    #% Reassign val to valfile
    valfile=outofrangeval*ones((nval,nc))
    valfile[0:nval,colnum]=val
 
    valfile[find(isnan(valfile))]=outofrangeval

    for i in xrange(0,nval) :
        fid.write(valfile[i,:],'\r')     #% Write each data line    

    fclose()
