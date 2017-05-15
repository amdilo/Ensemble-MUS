""" Ensemble generation for a simplified Made-Up sensor
    Project: H2020 FIDUCEO 
    Author: Arta Dilo, NPL M&M
    Date created: 07-03-2016
    Last update: 15-05-2017 
    
FCDR ensemble generation from a flat truth image and simple error structure, 
and visualise images and graphs of image values. 
Uses classes from classFCDR and visualisation functions from visFCDR. """

from classFCDR import *
from visFCDR import *


""" Generate measured image from an assumed truth, and plot values. """
cTruth = 80 # set a constant value for truth from range [15, 113]
Cict = 100 # set a constant value for true ICT count
Tict = 200 # set a constant value for true ICT temperature

imgNoR = 100 # number of rows/scanlines in the image
imgNoC = 100 # number of columns in the image
iSize = [imgNoR,imgNoC] # set image size: no of rows & columns

# perturb ICT counts and temperature
rnd = np.random.random(imgNoR)
CictT = np.ones(imgNoR) * Cict + rnd
TictT = np.ones(imgNoR) * Tict + 2 * rnd

# generate flat True image, i.e. constant Earth count value
tImg = Timage(iSize[0], iSize[1], cTruth, CictT, TictT) 
print 'Truth image values: ', tImg.getValArr(0)

cntU = 0.5  # set counts uncertainty
tmpU = 0.08 # set temperature uncertainty 

# create error (structure)
mErr = Error(cntU, tmpU) 
 # generate measured image 
mImg = Mimage(tImg, mErr)
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