'''
% GRAPHLIBtutorial          - tutorial for the graphlib directory (Jan 1,2001)

Translated by Dimitri D''Or - November 2004

Status: not completed: colorplot

'''

import matplotlib
matplotlib.use('GTK')
from matplotlib.matlab import *
from scipy import *
from scipy.stats import *

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


## %%% Clear memory content and set echo to on.

## clear;
## clc;
## echo on;

## %%% Load the data file called 'Falmagne.txt' and display 
## %%% its content (see the Falmagne.txt file for the meaning
## %%% of these variables). ch has two columns which are the  
## %%% the spatial coordinates (in meters). The next three
## %%% variables have the values for the sand, silt and clay contents (in %).
## %%% The last variable is a numeric code for the corresponding soil type.


(val,valname,filetitle) = readGeoEAS('../BMELIB2.0b/tutorlib/Falmagne.txt')

ch=val[:,0:2]
sand=val[:,2]
silt=val[:,3]
clay=val[:,4]
code=val[:,5]
dir()


## %%% Type any key for continuing...

## pause;
## clc;

## %%% Open a new graphic window and display the histograms
## %%% and the scatter plots for the sand; silt and clay
## %%% contents

figure(_fig())
histscatterplot(concatenate((reshape(sand,(len(sand),1)), reshape(silt,(len(silt),1)), reshape(clay,(len(clay),1))),1))

## %%% Type any key for continuing...

## pause;
## clc;

## %%% Open a new graphic window and display the values for
## %%% sand content using colorplot.m, with the 'hot' color
## %%% map and squares having a black edge color. Axes are
## %%% set to correct proportions

## figure(_fig())
## kwargs={'marker':'s','markersize':12,'markeredgecolor':[0,0,0]}
## zrange=[min(sand),max(sand)]
## colorplot(ch,sand,'hot',zrange,kwargs)
## axis('equal')

## %%% Type any key for continuing...

## pause;
## clc;

## %%% Open a new graphic window and display the values for
## %%% silt content using markerplot.m, with circles having
## %%% a blue edge color and a transparent face color. Axes
## %%% are set to correct proportions.

figure(_fig())
sizelim=[3, 20]
kwargs={'marker':'o','markeredgecolor':[0,0,1],'markerfacecolor':'none'}
zrange=[min(silt),max(silt)]
markerplot(ch,silt,sizelim,zrange,kwargs)
axis('equal')

## %%% Type any key for continuing...

## pause;
## clc;

## %%% Open a new graphic window and display :
## %%%
## %%% i)  in the left hand side the values for clay content using
## %%%     poleplot.m, with aquamarine circles having a black
## %%%     border and the default size ;
## %%%
## %%% ii) in the right hand side the value 0 for clay content below
## %%%     15 % and the value 1 for clay content above 15 % using a
## %%%     blue Times font with letter size equal to 9 points
## %%%
## %%% The graphic window can be manually resized for improving the
## %%% vizualization

figure(_fig())
subplot(1,2,1)
kwargs={'marker':'o','markeredgecolor':[0,0,0],'markerfacecolor':[127/255, 1, 212/255]}
poleplot(ch,clay,kwargs)
subplot(1,2,2)
kwargs={'color':[0.5, 0, 0.9],'fontname':'Times','fontsize':9}
valplot(ch,clay>15,kwargs)

show()

## %%% Type any key for continuing...

## pause;
## clc;

## %%% Save for further use the variables from the Falmagne.txt file
## %%% into a MATLAB binary data file called Falmagne.mat

## save Falmagne ch sand silt clay code

## %%% End of the GRAPHLIB tutorial

## echo off;

