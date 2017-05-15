""" FIDUCEO FCDR ensemble generation 
    Author: Arta Dilo / NPL MM
    Date created: 07-03-2016
    Last update: 24-10-2016 

Modelling errors around the measured values, potential problems.
Old version, combined content from three files: 
    class declarations from classens.py
    visualisation function declarations from visens.py
    calls for ensemble generation from ensemble.py. """

import random
import numpy as np
import copy as cp
import matplotlib.pyplot as plt
#import matplotlib.image as mpimg 
import matplotlib.mlab as mlab


""" This class describes a satellite image band. In this test we use syntetic  
images, simplified version of satellite images (one-band instead of multi-band). 
The class has variables/attributes:
    - number of rows
    - number of columns
    - a 2D array of image (pixel) values.
Later may change to an image class from GDAL, and use the corresponding library 
to read, manipulate and display images. """
class ImgBand(object):
    
    def __init__(self, rows, columns):
        """ Initialises an image with the specified number of rows and columns. 
        Both arguments, rows and columns, are integers. Initialise with 0. """
        
        self.nor = rows
        self.noc = columns
        self.pxs = np.zeros((self.nor, self.noc),dtype=float)
        
    def getNoR(self): # get number of rows in the image
        return self.nor
        
    def getNoC(self): # get number of columns in the image
        return self.noc

    def getIVals(self): # get image values as a (numpy) 2D array
        """ Return a copy of pixel values so that the calling functions do not 
        change the values in the original image. """
        ivals = self.pxs.copy() # copy of numpy (2D) array of pixel values
        return ivals
        
    def setIVals(self, mat): # set image values to mat 
        self.pxs = mat # mat should be an np.array of dimension: nor x noc
        return self.pxs
        
    def getPxVal(self, x, y): # get the value at pixel (x,y)
        val = cp.copy(self.pxs[x, y]) # copy of pixel value
        return val
        
    def setPxVal(self, x, y, val): # set the value of pixel (x,y) to val
        self.pxs[x, y] = val

    def getValArr(self): # get values as a 1D array;
        vals = self.pxs.flatten() # a copy of the array collapsed to 1D
        return vals
        
    def imgMap(self, mapt, mapc=None, zran=None):
        """ The method diplays the image using matplotlib. Its arguments are
        - mapt: map title (string type)
        - mapc: color scheme - label for the library color schemes (string)
        - zran: value range for colors (array of numbers, length 2) """
        if mapc is None:
            mapc = 'gist_earth'
        if zran is None:
            zran = [min(self.getValArr()), max(self.getValArr())]
        
        plt.figure() 
        plt.title(mapt)
        iM = plt.imshow(self.pxs, cmap=mapc, vmin=zran[0], vmax=zran[1])
        plt.colorbar(iM)
        plt.show()
        return iM
        

""" The Error class contains attributes that describe the error structure. 
For now the attributes are:
    - relative error at pixel level, e.g. random error of Earth counts 
    - relative error at scanline level, e.g. combination of random errors of 
        PRTs' temperature and ICT counts. 
Methods of the class draw an error from error PDFs of counts and temperature 
and calculate pixel and scanline error using sensitivity coefficients. """
class Error(object):
    
    def __init__(self, pxErr, slnErr):
        self.rpxErr = pxErr # relative error at pixel level
        self.rlnErr = slnErr # relative error at scanline level

    def getPxErr(self): # get error at pixel level
        attr1 = cp.copy(self.rpxErr) # for now a silly name: attr1 
        return attr1
        
    def getSlnErr(self): # get error at scanline level
        attr2 = cp.copy(self.rlnErr) # pythonic confusion if attr2 <-> slnErr ??
        return attr2


""" This is an ensemble instance. It is computed from the measurements image 
and its errors, for the moment a random error at pixel and at scanline level. 
The class has attributes:
    - a 2D array of instance values calculated from measurements and err struct
    - a separate 2D array for errors (in each pixel). 
Currently it has two ways for calculating instance values: measurement - errors; 
measurement + errors. """
class eInstance(ImgBand):
    
    def __init__(self, measImg, isize, ierr):
        """ Initialises an ensemble instance image: sets the size of the image
        to measurements image measImg size, and calculates the instance values  
        and errors from measurements and error components from ierr attributes,  
        for now the relative errors at pixel and scanline levels. """
        
        self.nor = isize[0] # number of rows
        self.noc = isize[1] # number of columns
        self.pxs = measImg.getIVals() # initialise instance values with measured 
        self.errs = np.zeros((self.nor, self.noc), dtype=float) # init error
        
        # compute instance values and errors (both 2D np arrays)
        rpxErr = ierr.getPxErr()  # relative err at pixel level
        rlnErr = ierr.getSlnErr() # relative err at scanline
        for i in range(self.nor): # is it better to use local vars NoR, NoC ??
            lnStd = rlnErr * self.pxs[i,1] # scanline i standard deviation  
            lnErr = random.gauss(0, lnStd)
            
            for j in range(self.noc):
                pxStd = rpxErr * self.pxs[i,j] # standard dev. for pixel (i,j)
                pxErr = random.gauss(0, pxStd)
                self.errs[i,j] = lnErr + pxErr
                self.pxs[i,j] -= self.errs[i,j] # remove error from measured
      
    def getInstErrs(self): # get instance errors as a 2D array
        ierrs = self.errs.copy() # copy of instance errors
        return ierrs
        
    def getPxErr(self, x, y): # get intance error at pixel (x,y)
        return cp.copy(self.errs[x,y])
        
    """ Perform an alternate computation of ensemble instances (for testing), 
    add errors to the measurements image: instance value = meas + err. """        
    def compAPxs(self, measImg):
        apxs = measImg.getIVals() + self.errs
        return apxs


""" This is an ensemble of instances calculated from the measured image and
the estimated errors. The class has attributes:
    - size of the ensemble, i.e. number of instances
    - the measured image the ensemble was generated from
    - size of an ensemble instance, i.e. no of rows & columns in an instance
    - error structure of the ensemble
    - list of eInstance object, i.e. ensemble instances 
    - GUM standard uncertainty image. """
class Ensemble(object):
    
    def __init__(self, noInst, measImg, measErr):
        self.measurements = measImg
        noR = measImg.getNoR() # number of rows
        noC = measImg.getNoC() # number of columns
        self.insize = [noR, noC] # instance size, i.e. no rows & columns 
        self.errStr = measErr # error structure of the ensemble

        self.size = noInst
        self.instanceSet = [] # list of instances 
        # create and add noInst new instances 
        for k in range(noInst):
            newIns = eInstance(measImg, self.insize, self.errStr)
            self.instanceSet.append(newIns)
        
        # create GUM uncertainty image
        self.gumU = ImgBand(noR,noC) # initialise uncertainty image
        # ensemble's constant relative uncert -> will CHANGE, based on Error 
        rU = ((measErr.getPxErr())**2 + (measErr.getSlnErr())**2)**0.5 
        unc = rU*abs(measImg.getIVals()) # full array uncertainty
        self.gumU.setIVals(unc)
    
    def getESize(self): # get the size of the ensemble, i.e. no of instances
        return self.size
        
    def getEInsize(self): # get the ensemble instance size
        return self.insize

    def getEnErr(self): # get the ensemble error structure
        return self.errStr

    def getGumU(self): # get ensemble uncertainty (an image)
        return self.gumU

    def getInstance(self, k): # get instance k of the ensemble
        if k >= self.size:
            print k, ' is greater than the size of the ensemble'
            return 0 # do a proper error handling
        else:
            return self.instanceSet[k]

    def getESource(self): # get measurements (image) that built the ensemble
        return self.measurements
            
    def getPxEval(self, x, y): # get ensemble values and errors at pixel (x,y)
        noEI = self.size # number of instances
        evals = np.empty([2, noEI], dtype=float) 
        
        for k in range(noEI):
            instance = self.getInstance(k)
            evals[0,k] = instance.getPxVal(x,y)
            evals[1,k] = instance.getPxErr(x, y)
        return evals
                
""" This is the assumed truth image, which value we don't know in reality. 
Values will be simulated, here a constant value. It is subclass of ImgBand. """
class Timage(ImgBand):
    
    def __init__(self,rows, columns, value):
        # fill all pixels with the given value
        self.nor = rows
        self.noc = columns        
        self.pxs = np.full((self.nor, self.noc), value, dtype=float)
        
        
""" This image contains the measured values. For this test we simulate the 
measured values image from an assumed truth image and simple errors, 
a random error at pixel level, a random error at scanline. 
The class is subclass of ImgBand, thus has attributes:
    - number of rows
    - number of columns
    - a 2D array of pixel values """
class Mimage(ImgBand):
    
    def __init__(self, trueImg, measErr):
        """ Initialises the measured image: set the size, ie.e. rows and columns 
        to truth image trueImg size and simulate measured values for each pixel
        using truth and error structure in measErr. """
        
        self.nor = trueImg.getNoR()
        self.noc = trueImg.getNoC()
        self.pxs = trueImg.getIVals() # set pixel values to trueImg values 
        pixelErr = measErr.getPxErr() # get err at pixel level
        scanlineErr = measErr.getSlnErr() # get err at scanline
        
        nRow = self.nor     # number of rows
        nCol = self.noc     # number of columns
        for i in range(nRow):
            """ Standard deviation calculated from relative scanline error and
             the first element in the row (this last will change to ...) """
            lnStd = scanlineErr * self.pxs[i,1] # standard dev. at scanline i
            lnErr = random.gauss(0, lnStd)
            
            for j in range(nCol):
                pxStd = pixelErr * self.pxs[i,j] # standard dev. at pixel (i, j)
                pxErr = random.gauss(0, pxStd)
                self.pxs[i,j] += lnErr + pxErr 
        
## END OF classens


## START OF visens (functions declaration)
""" Plot value and error image for an ensemble instance. """
def mapInstance(ensemble, Idx, zran):
    # get ensemble instance
    inst = ensemble.getInstance(Idx)
    emeas = ensemble.getESource()

    # visualise image of instance values (meas - err)
    plt.figure() 
    plt.title('Ensemble instance %d' %Idx)
    imgV = plt.imshow(inst.getIVals(),cmap='BuGn',vmin=zran[0],vmax=zran[1]) 
    plt.colorbar(imgV)
    plt.show()
    # visualise image of alternate instance values (meas + err)
    plt.figure() 
    plt.title('Alternate values of instance %d' %Idx)
    imgAV = plt.imshow(inst.compAPxs(emeas),cmap='BuGn',vmin=zran[0],vmax=zran[1]) 
    plt.colorbar(imgAV)
    plt.show()
    
    # visualise image of instance errors
    plt.figure() 
    plt.title('Error in ensemble instance %d' %Idx)
    imgEr = plt.imshow(inst.getInstErrs(),cmap='OrRd') 
    plt.colorbar(imgEr)
    plt.show()

        
""" Plot ensemble distribution of values and errors for a given pixel """
def plotPxDist(ensemble, x, y):
    # get value and error of pixel (x,y) through all ensemble instances
    pxDist = ensemble.getPxEval(x, y)
    pxVals = pxDist[0,:]
    pxErrs = pxDist[1,:]
    # get the measured value and uncertainty at pixel (x,y)
    pxMeas = ensemble.getESource().getPxVal(x, y)
    pxUnc = ensemble.getGumU().getPxVal(x, y)
    pxEaM = pxMeas + pxErrs # errors around meared pixel value
    
    """ Plot histogram of ensemble values at pixel (x,y)"""
    no_bin = 15 # number of bins for the histograms
    plt.figure()
    plt.hist(pxVals, no_bin, normed=1, fc='green', alpha=0.5)
    title = 'Histogram of values at pixel ' + str((x,y))
    plt.title(title)
    plt.xlabel('Instance number')
    plt.ylabel('Percentage') 
    plt.show()
    
    """ Plot histogram of ensemble errors at pixel (x,y) """
    plt.figure()
    plt.hist(pxErrs, no_bin, normed=1, fc='red', alpha=0.5)
    title = 'Histogram of errors at pixel ' + str((x,y))
    plt.title(title)
    plt.xlabel('Instance number')
    plt.ylabel('Percentage')
    plt.show()
    
    """ Plot graph of ensemble values and errors at pixel (x,y) """
    plt.figure() 
    pv, = plt.plot(pxVals, 'go', label='Ens. values') # does not come in legend
    pav, = plt.plot(pxEaM, 'ro', alpha=0.4, label='Ens. alt. val') 
    # plot measured +/- standard uncertainty
    pm = plt.axhline(pxMeas, linewidth=3, color='c') # plot mean 
    pbnd = plt.axhspan(pxMeas - pxUnc, pxMeas + pxUnc, fc='c', alpha=0.4)
    # add legend to the graph
    plt.legend([pv,pav,pm,pbnd],['Ens. values','Alternate vals','Measured val','+/- Std uncert.'])
    title = 'Ensemble values & alternates at pixel ' + str((x,y))
    plt.title(title)
    plt.xlabel('Instance no.')
    plt.ylabel('Inst. (alternate) values')
    plt.show()
    

""" Plot images, scatterplot and histogram of measured and truth values. 
The graphs, scatterplot & histogram, are meaningful i.e. useful in some way, 
for the case when we start from a flat truth. Others should be thought of 
for a general case, i.e. satellite image of a non-homogenous region. """
def plotMeasured(truth, measured, error, zran=None):

    """ Collect values for plotting the graph and histogram """
    atruth = truth.getPxVal(0,0) # constant value used for truth
    # Constant uncertainty of measured from constant truth and relative errs
    mU = abs(atruth) * (error.getPxErr()**2 + error.getSlnErr()**2)**0.5

    # Summary statistics of the (simulated) measured image values
    mVals = measured.getValArr() # image values as 1D array
    meanM = np.mean(mVals) # mean of measured
    stdM = np.std(mVals) # standard deviation of measured 
    minM = min(mVals) # min and max measured value
    maxM = max(mVals)

    # Print summary data about measured image
    print '\nSample mean of measured values:', meanM 
    print 'Sample standard deviation of measured:', stdM
    print 'Standard uncertainty of measured image:', mU
    
    """ Plot the (simulated) truth and measured image """
    zran = [minM, maxM]
    truth.imgMap('Assumed truth image','Greens',zran)    
    measured.imgMap('Measured generated from flat truth','Greens',zran) 

    """ Plot graph of measured values and stats """
    plt.figure() 
    mv, = plt.plot(mVals, 'g+', label='Measured') # does not come in legend
    plt.title('Measured values generated from a flat truth')
    plt.xlabel('Pixel count')
    plt.ylabel('Pixel value')
    # plot mean of measured +/- standard deviation
    sm = plt.axhline(meanM, linewidth=3, color='c') # plot mean 
    sbnd = plt.axhspan(meanM - stdM, meanM + stdM, fc='c', alpha=0.4)
    # plot truth +/- constant (measured) uncertainty
    at = plt.axhline(atruth, linewidth=2, color='r') # plot truth value
    plt.axhline(atruth-mU, linewidth=3, color='b') 
    gumub = plt.axhline(atruth+mU,linewidth=3,color='b') # plot truth + uncert.
    # add legend to the graph
    plt.legend([mv,sm,sbnd,at,gumub],['Measured','Mean measured','+/-Standard dev.','Truth','+/- Std uncert.'])
    plt.show()
    
    """ Plot histogram of measured values """
    plt.figure()
    histB = max(10,np.ceil(3*len(mVals)**0.3)) # no. of bins for the histogram
    n, bins, patches = plt.hist(mVals, histB, normed=1, facecolor='blue', alpha=0.5)
    ym = mlab.normpdf(bins, meanM, stdM) # add a 'best fit' line
    plt.plot(bins, ym, 'm--')
    plt.title('Histogram of measured values')
    plt.xlabel('Pixel values')
    plt.ylabel('Percentage') # values don't look right, expected y in [0,1]
    plt.show()
    
        
""" Plot images, scatterplot and histograms for ensemble. The graphs have a  
meaning when we start from a flat truth. Others should be thought of for the 
general case, i.e. satellite image of a non-homogenous region. """
def plotEnsemble(ensemble, INo):
    
    """ Get measurements image and calculate summary stats """ 
    measured = ensemble.getESource() # get measurements image 
    mVals = measured.getValArr() # measurements as 1D array
    meanM = np.mean(mVals) # mean of measured values
    stdM = np.std(mVals) # standard deviation of measured

    histB = max(10,np.ceil(2*len(mVals)**0.3)) # no of bins for histograms
    # maybe will need length instead of number of bins
   
    """ Collect ensemble values and summary stats for graphs and histograms """
    # calculate the number of instances' values
    ilen = measured.getNoR() * measured.getNoC() # elements in an instance
    elen = len(INo) * ilen # elements in INo ensemble instances
    k = 0 # count instances in the ensemble
    
    # check for simple code (concatenate, append?) without waste of memory...
    insVals = np.zeros([elen]) # array of values of ensemble instances
    insErrs = np.zeros([elen]) # array of instances' errors
    insAVals = np.zeros([elen]) # array of instances' additional vals
    for iid in INo:
        instance = ensemble.getInstance(iid)

        vals = instance.getValArr() # get instance values
#        print 'Instance', iid, 'values:', vals
        insVals[k*ilen:(k+1)*ilen] = vals

        errs = instance.getInstErrs().flatten() # get instance errors
        insErrs[k*ilen:(k+1)*ilen] = errs
        avals = instance.compAPxs(measured).flatten() # get alternate values
        insAVals[k*ilen:(k+1)*ilen] = avals
        k += 1
#    ehistB = max(10,np.ceil(3*len(insVals)**0.3)) # no bins for histogram
    
    # calculate and print ensemble summary stats
    meanpE = np.mean(insVals) # mean & stdev of picked instances (calc: meas-err)
    stdpE = np.std(insVals) 
    minp = min(insVals) # min & max pixel value through all picked instances
    maxp = max(insVals)
    zran = [minp, maxp] # color-map range for instance values

    meanE2 = np.mean(insAVals) # mean & stdev of instances (calc: meas+err)
    stdE2 = np.std(insAVals) 
    
    # Get ensemble uncertainty (image), max & min uncertainty value
    uncert = ensemble.getGumU() # GUM uncertainty image
    eU = max(uncert.getValArr()) # max uncertainty value
    emU = np.mean(uncert.getValArr()) # average uncertainty value
    
    esize = ensemble.getESize() # ensemble size
    print '\nPicked no from the ensemble of', esize, 'instances:', INo
    print 'Sample mean of instances (meas-err, meas+err):', meanpE,',',meanE2
    print 'Sample st.dev of instances (meas-err, meas+err):', stdpE,',',stdE2
    print 'Maximum standard uncertainty of ensemble:', eU
    print 'Average standard uncertainty * sqrt(2):', emU * np.sqrt(2)

    """ Plot images: measured, GUM uncertainty, value & error of 2 instances """
    measured.imgMap('Measured image creating the ensemble','BuGn',zran) 
    uncert.imgMap('Ensemble\'s GUM uncertainty','Oranges') # display uncert.map
    mapInstance(ensemble, INo[0], zran) # map of instance value and error 
    
    """ Scatterplot of picked ensemble instances' values """
    plt.figure()
    ev, = plt.plot(insVals, 'g+', label='Ensemble values', zorder=1) 
    em = plt.axhline(meanpE, linewidth=3, color='c') # mean of picked inst. vals  
    sibnd = plt.axhspan(meanpE-stdpE, meanpE+stdpE, fc='c', alpha=0.4, zorder=2)
    mm = plt.axhline(meanM, linewidth=2, color='b') # mean of measured
    gumeub = plt.axhline(meanM+eU,linewidth=3,color='m') # mean meas.+GUM uncert
    plt.axhline(meanM-eU, linewidth=3, color='m') 
    plt.legend([ev,em,sibnd,mm,gumeub],['Ensemble values','Ensemble Mean','+/- Ensemble Std','Mean measured','+/- max std uncert.'])
    title = 'Pixel values of ensemble instances: ' + str(INo)
    plt.title(title)
    plt.xlabel('Pixel count')
    plt.ylabel('Pixel value')
    plt.show()
    
    """ Histograms of measured & ensemble instances (meas+errs, meas-errs) """
    # histogram to compare distribution of measured with ensemble values
    plt.figure()
    n, bins, patches = plt.hist(mVals,histB,normed=1,fc='blue',alpha=0.5)#,rwidth=histW)
    ym = mlab.normpdf(bins, meanM, stdM) # gaussian 'best fit' for measured
    plt.plot(bins, ym, 'm--')
    en, ebin, epat = plt.hist(insVals,histB,normed=1,fc='green',alpha=0.3)#,rwidth=histW)
    ye = mlab.normpdf(ebin, meanpE, stdpE) # gaussian 'best fit' for ensemble
    plt.plot(ebin, ye, 'r--')
    title = 'Histogram of measured and instance ' + str(INo) + ' values'
    plt.title(title)
    plt.xlabel('Pixel values')
    plt.ylabel('Percentage')
    plt.show()
    
    # hist. comparing the 2 ways of calculating ensemble from measured & errors
    plt.figure()
    data = [mVals, insVals, insAVals] 
    datalbl = ['Measured values','Ensemble values','Ens.alternate vals']
    plt.hist(data, histB, alpha=0.5, label=datalbl)
    plt.title('Histogram of measured and ensemble values')
    plt.xlabel('Pixel values')
    plt.ylabel('Counts')
    plt.legend()
    plt.show()
            
## END OF funens


""" Generate measured image from an assumed truth, and plot values. """
ctruth = 1 # set a constant value for truth
iSize = [100,100] # set image size: no of rows & columns
tImg = Timage(iSize[0], iSize[1], ctruth) # create truth image
print 'Truth image values: ', tImg.getValArr()

rpxe = 0.02 # set relative error at pixel level 
rsle = 0.03 # set relative error at scanline level 
mErr = Error(rpxe, rsle) # create error (structure)
print 'Error structure (pixel err, scanline err):', mErr.getPxErr(), mErr.getSlnErr()

mImg = Mimage(tImg, mErr) # generate measured image 
print 'Measured image values: ', mImg.getValArr()
plotMeasured(tImg, mImg, mErr) # visualise measured image and graphs of values


""" Generate ensemble from the measured image with a given ensemble size, 
pick a number of instances and pixels for displaying images and graphs. """
eSize = 80 # set ensemble size, i.e. number of instances to generate
if eSize < 10: # ens. size will be an input var, limited to range [10, 100]
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