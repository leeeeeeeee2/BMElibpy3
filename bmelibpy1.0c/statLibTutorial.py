'''
% STATLIBtutorial          - tutorial for the statlib directory (Jan 1,2001)

Translated by Didrik Pinte - November 2004

Status: not completed: coregfit

'''

import matplotlib
matplotlib.use('TkAgg')
# from matplotlib.matlab import *
from matplotlib.pylab import *
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

(val,valname,filetitle) = readGeoEAS('../BMELIB2.0b/tutorlib/Falmagne.txt')

ch=val[:,0:2]
sand=val[:,2]
silt=val[:,3]
clay=val[:,4]
code=val[:,5]

who()

##raw_input('Pause - Press on enter to continue')

##%%% Compute and plot the density scaled histograms for the
##%%% sand, silt and clay contents. 

nbins=20
bounds=(0,100)

figure(_fig())
subplot(221)
histscaled(sand,nbins,bounds)
title('Sand')
subplot(222)
histscaled(silt,nbins,bounds)
title('Silt')
subplot(223)
histscaled(clay,nbins,bounds)
title('Clay')

show()

##%%% Compute the basic statistics for the sand, silt and clay
##%%% contents, respectively.

print "mean" ,
print map(mean,[sand,silt,clay])
print "variance",
print  map(var, [sand,silt,clay])
print "skewness",
print map(skew, [sand,silt,clay])
print "kurtosis",
print map(kurtosis,[sand,silt,clay])

##%%% Compute and plot the kernel
##%%% density estimates evaluated from 0 to 100 % by step of 1%.
##%%% The variance of the Gaussian kernel is equal to 20.

sandfileK=arrayrange(-15,55.5,0.5,Float32)
siltfileK=arrayrange(25,100.5,0.5,Float32)
clayfileK=arrayrange(-10,56,1,Float32)
v = 20.0


pdfsandfileK=kerneldensity(sand,sandfileK,v);
pdfsiltfileK=kerneldensity(silt,siltfileK,v);
pdfclayfileK=kerneldensity(clay,clayfileK,v);


subplot(221)
plot(sandfileK,pdfsandfileK,'r')
subplot(222)
plot(siltfileK,pdfsiltfileK,'r')
subplot(223)
plot(clayfileK,pdfclayfileK,'r')


##%%% Estimate and plot the cumulative distribution function for
##%%% the sand, silt and clay contents. Superimpose on the graph
##%%% the cumulative distribution function obtained by integrating
##%%% the kernel density estimates

sandfileE,cdfsandfileE=cdfest(sand);
siltfileE,cdfsiltfileE=cdfest(silt);
clayfileE,cdfclayfileE=cdfest(clay);

cdfsandfileK=pdf2cdf(sandfileK,pdfsandfileK);
cdfsiltfileK=pdf2cdf(siltfileK,pdfsiltfileK);
cdfclayfileK=pdf2cdf(clayfileK,pdfclayfileK);

axesDef = [0, 100, 0, 1]

figure(_fig())
subplot(221)
plot(sandfileE,cdfsandfileE)
plot(sandfileK,cdfsandfileK,'r')
title('Sand')
axis(axesDef)
subplot(222)
plot(siltfileE,cdfsiltfileE)
plot(siltfileK,cdfsiltfileK,'r')
title('Silt')
axis(axesDef)
subplot(223)
plot(clayfileE,cdfclayfileE)
plot(clayfileK,cdfclayfileK,'r')
title('Clay')
axis(axesDef)


##%%% Compute the directional variograms for the sand content along
##%%% the N-S (stars) and the W-E (circles) directions, using an
##%%% angular tolerance of 90 degrees. Distances classes are 500
##%%% meters wide and are ranging from 0 to 10000 meters.


print "Computing the directional variograms"

cl= range(0,10500,500)

optionsNS=[0, 45, -45]
dsandNS,vsandNS,osandNS = vario(ch,sand,cl,'loop',optionsNS)

optionsWE=[1,-45,45]
dsandWE,vsandWE,osandWE=vario(ch,sand,cl,'loop',optionsWE)

##%% Compute the omnidirectional variograms for the sand content.
##%%% Distances classes are the same than above.

dsand,vsand,osand=vario(ch,sand[:,NewAxis],cl,'loop')

figure(_fig())

plot(dsandNS,vsandNS,'o')
plot(dsandWE,vsandWE,'s')
plot(dsand,vsand,'+b')
plot(dsand,vsand)
axis([0, 100, 0, cl[-1]])

 
##%%% Plot a nested variogram models, composed of a nugget effect
##%%% with a sill equal to 30, and an exponential model with a sill
##%%% equal to 60 and a range equal to 3000.

modelsand=[nuggetV,exponentialV]
paramsand0=[30.0,[60.0,3000.0]]
modelplot(dsand,modelsand,paramsand0,'Color',[0,0,1])

##%%% Type any key for continuing...

##%%% Fit by a weighted least squares method the estimated
##%%% variogram using a nested model that includes a nugget
##%%% effect and an exponential model, and display the results.
##%%% The initial values for the parameters are as above.
paramsand=modelfit(dsand,vsand,osand,modelsand,paramsand0);
modelplot(dsand,modelsand,paramsand,'Color',[1,0,0]);
print "Optimized parameters"
print paramsand

##%%% Compute the whole set of omnidirectional variograms and
##%%% cross variograms for the raw values of the sand, silt
##%%% and clay contents. Distance classes are the same than above.

options=[1];
figure(_fig())
[d,V,o]=vario(ch,concatenate((sand[:,NewAxis],silt[:,NewAxis],clay[:,NewAxis]),1),cl,'loop',options);


##%%% Fit by a weighted least squares method the whole set of
##%%% variograms and cross variograms, using a nested model that
##%%% includes a nugget effect and an exponential model with range
##%%% equal to 3000. Display the estimated sill coefficients for
##%%% the nugget effect and the exponential models

show()

model=[nuggetV,exponentialV]
param0=[[30.0],[30.0,3000.0]]
options=[1]
[param]=coregfit(d,V,o,model,param0,options)

#show()
