""" FIDUCEO FCDR generation 
    Author: Arta Dilo, NPL M&M
    Date created: 07-03-2016
    Last update: 28-10-2016 
    
FCDR ensemble generation from a flat truth image and simple error structure, 
and visualise images and graphs of image values. 
Uses classes from classFCDR and (visualisation) functions from visFCDR. """

from classFCDR import *
from visFCDR import *


""" Generate measured image from an assumed truth, and plot values. """
ctruth = 100 # set a constant value for truth from range [15, 113]
iSize = [100,100] # set image size: no of rows & columns
tImg = Timage(iSize[0], iSize[1], ctruth) # create truth image
print 'Truth image values: ', tImg.getValArr()

cntU = 0.5  # set counts uncertainty
tmpU = 0.08 # set temperature uncertainty 
mErr = Error(cntU, tmpU) # create error (structure)
print 'Error structure (pixel err, scanline err):', mErr.getPxErr(), mErr.getSlnErr()

mImg = Mimage(tImg, mErr) # generate measured image 
print 'Measured image values: ', mImg.getValArr()
plotMeasured(tImg, mImg, mErr) # visualise measured image and graphs of values


""" Generate ensemble from the measured image with a given ensemble size, 
pick a number of instances and pixels for displaying images and graphs. """
#eSize = int(input("Enter the ensemble size: ")) # give ensemble size
eSize = 100 # set ensemble size, i.e. number of instances to generate
if eSize < 10: # limit size to range [10, 100]
    eSize = 10
elif eSize > 100:
    eSize = 100
print '\nGenerating an ensemble with', eSize, 'instances...'

ensImgs = Ensemble(eSize, mImg, mErr) # generate ensemble

# Pick a random pixels and show the distribution of values
ix = random.choice(range(iSize[0])) # row index
iy = random.choice(range(iSize[1])) # column index
plotPxDist(ensImgs, ix, iy) # display histograms of values and errors

# Visualise images of ensemble instances and graphs of values
#Ino = range(eSize)
pick = 5 # set number of instances to pick from the ensemble
Ino = random.sample(range(eSize), pick) # get random ids for instance picking 
plotEnsemble(ensImgs, Ino)