'''
% GENLIBtutorial            - tutorial for the genlib directory (Jan 1,2001)

Translated by Dimitri D''Or - October 2004

Status: not completed

'''

import matplotlib
matplotlib.use('GTK')
from matplotlib.matlab import *
from scipy import *
from scipy.stats import *
from scipy.interpolate import *

from iolib import *
from genlib import *
from graphlib import *
from modellib import *
from statlib import *


_id=0

def _fig(*id):
    global _id
    _id+= 1
    return _id


##%%% Load the data file called 'Falmagne.txt' and display 
##%%% its content (see the Falmagne.txt file for the meaning
##%%% of these variables).

(val,valname,filetitle) = readGeoEAS('../BMELIB2.0b/tutorlib/Falmagne.txt')

ch=val[:,0:2]
sand=val[:,2]
silt=val[:,3]
clay=val[:,4]
code=val[:,5]

##%%% Open a new graphic window and display as circles the locations
##%%% given in the file with correct axes proportions for the graph,
##%%% then freeze the graphic.

figure(_fig())
plot(ch[:,0],ch[:,1],'o')
axis('equal')
title('Data points and estimation grid')
xlabel('Easting (m)')
ylabel('Northing (m)')
hold(True)

##%%% Create a two dimensional grid with regular spacing along the
##%%% axes and superimpose this grid over the area. The grid spacing
##%%% is 1000 meters along both axes.

minc=array([178000,90000])
dc=array([1000,1000])
nc=array([17,19])
ck=creategrid(minc,dc,nc)

plot(ck[:,0],ck[:,1],'+r')
axis([178000,194000,90000,108000])

##%%% Perform an estimation of the sand content at the nodes of the
##%%% grid using two methods :
##%%%
##%%% i)  a nearest neighbour estimation ;
##%%%
##%%% ii) a kernel smoothing with a variance of the Gaussian kernel
##%%%     equal to 1e6.
##%%%
##%%% In both cases, the values taken into account for the estimations
##%%% are the 20 closest values, as long as their distance from the
##%%% estimation location does not exceed 50000 meters. Only the values
##%%% of the sand content over the area are taken into account for the
##%%% estimation.

nhmax=20
dmax=50000

zkID=invdist(ck,ch,sand,Inf,nhmax,dmax)
zkKS=kernelsmoothing(ck,ch,sand,1e6,nhmax,dmax)

##%%% Open a new graphic window and display the result of the estimation
##%%% for both methods using pcolor.m, after a conversion for the data
##%%% format has been performed with col2mat.m

figure(_fig())
[ck1,ck2,ZkID]=col2mat(ck,zkID)
[CK1,CK2]=meshgrid(ck1,ck2)
pcolor(CK1,CK2,ZkID,cmap=cm.jet)
colorbar(tickfmt='%1.1f')
axis('equal')
axis([178000,194000,90000,108000])
title('Nearest neighbour')
xlabel('Easting (m)')
ylabel('Northing (m)')

figure(_fig())
[ck1,ck2,ZkKS]=col2mat(ck,zkKS)
pcolor(CK1,CK2,ZkKS,cmap=cm.jet)
colorbar(tickfmt='%1.1f')
axis('equal')
axis([178000,194000,90000,108000])
title('Kernel smoothing')
xlabel('Easting (m)')
ylabel('Northing (m)')
show()

print('End of the GENLIB tutorial')





