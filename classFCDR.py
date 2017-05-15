""" FIDUCEO test FCDR generation 
    Author: Arta Dilo, NPL M&M
    Date created: 07-03-2016
    Last update: 31-10-2016 

CLASS declarations for working with FCDRs, full, easy, ensemble formats: 
test FCDR generation with different error structure... 
Classes are data structures that shape the information to be stored in 
full-FCDR format, and functions to extract and compile information components
to be stored in easy-FCDR and ensemble-FCDR formats. Classes include
    - raw satellite image, i.e. counts (level L0 data), here a simplified 
    format/content and only one image band; with subclasses 
        * truth image, assumed, in reality not known
        * measured image, generated from the assumed truth
    - error structure <- replicated from raw image structure (?)
    - harmonisation data
    - full FCDR 
    - easy FCDR 
    - ensemble instance
    - ensemble FCDR
    - satellite instrument describing a meta level of satellite sensor 
    information, e.g. number of bands, telemetry data, bands' SRFs, etc.
In this test images are one-band, used for seeing results of GUM application  
around the measured values (assuming a homegenous Earth region is observed). 
...Will probably be modified to multi-band images so that we account for 
accross channel correlation. """

import random
import copy as cp
import matplotlib.pyplot as plt
#import matplotlib.image as mpimg 
from funFCDR import *


""" This class describes a satellite image band, DN values and calibration data. 
Syntetic images are used here for testing purpose, calibration data is processed 
from raw measurements so that it is (simpler and) close to harmonisation vars. 
The class has variables/attributes:
    - number of rows
    - number of columns
    - a 2D array of Earth counts (space clamped), i.e. image DN values,
    - array of (space clamped) internal calibration target (ICT) counts
    - array of ICT temperatures, i.e. average of PRTs' temperatures.
This may change to an image class from GDAL, if there is a class for L0 data,
and use GDAL library to read, manipulate and display satellite images. """
class ImgBand(object):
    
    def __init__(self, rows, columns):
        """ Initialises an image with the specified number of rows and columns. 
        Both arguments, rows and columns, are integers. Initialise image values 
        and calibration data arrays with 0 (type float). Image rows are assumed
        to be sensor scanlines. """
        
        self.nor = rows
        self.noc = columns
        self.pxs = np.zeros((self.nor, self.noc),dtype=float)
        self.cnts = np.zeros(self.nor,dtype=float)
        self.tmps = np.zeros(self.nor,dtype=float)
        
    def getNoR(self): # get number of rows in the image
        return self.nor
        
    def getNoC(self): # get number of columns in the image
        return self.noc

    def getDNs(self): # get image DNs, i.e. Earth counts as numpy 2D array
        """ Return a copy of pixel values so that the calling functions do not 
        change the values in the original image. """
        ivals = self.pxs.copy() # copy of numpy (2D) array of pixel values
        return ivals
        
    def getCalDt(self): # get image calibration data
        cvals = (self.cnts.copy(), self.tmps.copy())
        return cvals
        
    def evaLE(self): # create radiance image from the calibration model
        Limg = np.zeros((self.nor, self.noc),dtype=float)
        """ fill values in nested loop for image rows, i.e. scanlines: 
        assuming LEarth function works with arrays """
        coef = a # calibration coeff.; global var from funFCDR
        ch = ch11 # channel parameters; global var from funFCDR
        for i in range(self.nor): # loop through rows
            CEvals = cp.copy(self.pxs[i, :]) # Earth counts for scanline i
            Cict = cp.copy(self.cnts[i])    # ICT count for scanline i
            Tict = cp.copy(self.tmps[i])    # ICT temperature at scanline i
            Limg[i,:] = LEarth(CEvals, Cict, Tict, coef, ch)
        return Limg
        
    def setIVals(self, mat): # set image DN and calibration values
        """ mat is an np.array of dimension: nor x (noc + 2) - NOT USED """
        self.pxs = mat # CHANGE this to image DNs and calibration arrays
        return self.pxs
        
    def getPxVal(self, x, y): # get the DN value at pixel (x,y)
        val = cp.copy(self.pxs[x, y]) # copy of pixel value
        return val
#    def setPxVal(self, x, y, val): # set the DN of pixel (x,y) to val
#        self.pxs[x, y] = val
#
    def getValArr(self, level): # get values as a 1D array;
        if level == 1:
            img = self.evaLE()
        else:
            img = self.pxs
        vals = img.flatten() # a copy of the image collapsed to 1D
        return vals
        
    def imgMap(self, level, mapt, mapc=None, zran=None):
        """ The method diplays the image using matplotlib. Its arguments are
        - level: plot radiance imge if 1, plot counts if 0 and otherwise
        - mapt: map title (string type)
        - mapc: color scheme; label for the library color schemes (string)
        - zran: value range for colors (array of numbers, length 2) """
        if mapc is None:
            mapc = 'gist_earth'
        if zran is None:
            zran = [min(self.getValArr()), max(self.getValArr())]
        if level == 1:
            img = self.evaLE()
        else:
            img = self.pxs
        
        plt.figure() 
        plt.title(mapt)
        iM = plt.imshow(img, cmap=mapc, vmin=zran[0], vmax=zran[1])
        plt.colorbar(iM)
        plt.show()
        return iM
        
        
""" This is an assumed truth image (truth is never known). A homogenous image 
is simulated for the test, thus a constant value for all image pixels and 
constant calibration values, i.e. constant ICT count and ICT temperature  
for all scanlines. Timage is a subclass of ImgBand. """
class Timage(ImgBand):
    
    def __init__(self,rows, columns, CEval, ICnt, ICTmp):
        # fill image pixels and calibration arrays with the given constants
        self.nor = rows
        self.noc = columns        
        self.pxs = np.full((self.nor, self.noc), CEval, dtype=float)
        self.cnts = np.full(self.nor, ICnt, dtype=float)
        self.tmps = np.full(self.nor, ICTmp, dtype=float)
        

""" The Error class contains attributes that describe the error structure, i.e. 
standard uncertainty for each error type of each measured variable. Here 
the variables 
For now the attributes are:
    - counts standard uncertainty, random component
    - (PRTs) temperature standard uncertainty, random component
    - temperature standard uncertainty, systematic component. 
Methods of the class store the standard uncertainty for each of the measured 
variables stored in full-FCDR-full,  from PDFs of counts and temperature 
and calculate pixel and scanline error using sensitivity coefficients. """
class Error(object):
    
    def __init__(self, uCnt, uPRT):
        self.uCounts = uCnt # counts uncertainty
        self.uTemp = uPRT # ICT temperature uncertainty

    def getPxErr(self): # get Earth count error: pixel level
        Cstd = cp.copy(self.uCounts) * np.sqrt(2) # standard dev for the distrib.
        pxErr = random.gauss(0, Cstd) # draw error from the normal PDF
        return pxErr
        
    def getSlnErr(self): # get error at scanline level
        Tstd = cp.copy(self.uTemp) # normal distrib. st.dev for ICT temperature
        Cstd = cp.copy(self.uCounts) * np.sqrt(2) # st.dev for ICT count
        return attr2


""" This image contains the measured values. For this test we simulate the 
measured values image from an assumed truth image and an error structure, (now) 
random error at pixel level, random error at scanline. It is subclass of ImgBand.
This will be image data from an instrument we work with in Fiduceo. Probably a 
specialisation of a GDAL class for multi-band image/raster (adding 
satellite instrument...). """
class Mimage(ImgBand):
    
    def __init__(self, trueImg, measErr):
        """ Initialises the measured image: set the size (rows & columns) to
        truth image size (trueImg) and simulates measured values for each pixel
        using truth and error structure in measErr. """
        
        self.nor = trueImg.getNoR()
        self.noc = trueImg.getNoC()
        self.pxs = trueImg.getIVals() # set pixel values to trueImg values 
        scanlineErr = measErr.getSlnErr() # get err at scanline
        
        nRow = self.nor     # number of rows
        nCol = self.noc     # number of columns
        for i in range(nRow):
            """ Standard deviation calculated from relative scanline error and
             the first element in the row (this last will change to ...) """
            lnStd = scanlineErr * self.pxs[i,1] # standard dev. at scanline i
            lnErr = random.gauss(0, lnStd)
            
            for j in range(nCol):
                pxErr = measErr.getPxErr() # get err at pixel level
                self.pxs[i,j] += lnErr + pxErr 


class hInput(object):
    # TO CHANGE
    def __init__(self, CE, Cict, Tict):
        self.ECount = CE # Earth counts 
        self.Cict = Cict # ICT counts 
        self.Tict = Tict # ICT temperature 


""" This is an ensemble instance. It is computed from the measurements image 
and its errors, (for now) a random error at pixel and at scanline level. 
The class has attributes:
    - a 2D array of instance values (calculated from measurements and err struct)
    - a separate 2D array for errors (in each pixel). 
Currently it has two ways for calculating instance values: measurement - errors; 
measurement + errors. """
class eInstance(ImgBand):
    
    def __init__(self, measImg, isize, ierr):
        """ Initialises an ensemble instance image: sets the size of the image
        to measImg size (measurements), and calculates the instance values and 
        errors from measurements and error components from ierr (attributes),  
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
            # Use Weibull distribution instead of Gaussian
#            lnErr = random.weibullvariate(1, lnStd)
            
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


""" This is an ensemble (of instances) calculated from the measured image and
the estimated errors. The class has attributes:
    - the measured image the ensemble was generated from
    - size of the ensemble, i.e. number of instances
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


""" Class instrument contains meta-information about a satellite instrument. 
The class defines the structure of data to be stored in full-FCDR format. Not 
sure yet how to use it for setting automatically the structure of full-FCDR.
Its attributes define the content of full-FCDR. They are:
    - number of bands
    - spectral response function of each band
    - structure of telemetry, i.e. raw measurement variables
    - independent variables for the harmonisation derived from telemetry vars
        * Earth and calibration counts per band
        * internal calibration target (ICT/IWCT) temperature
    - error structure inferred from the structure of harmonisation vars 
        * error in counts (different per count type, i.e. Earth, ICT, space;
                            and band ??) 
        * error in ICT temperature
    - ... """
class instrument(object):
    
    def __init__(self, noB, srfs, telemetry, harVars, errVars):
        self.NoB = noB # number of bands; integer
        self.SRFs = srfs # array of arrays(?) with SRF values
        self.RawMeas = telemetry # list of raw measurement vars
        self.HData = harVars # list of harmonisation variables
        self.HUnc = errVars # list of error variables for harmonisation
        
                
# digitisation noise: err = 0.5*(2**b), with b the number of bits (i.e. length)