""" FIDUCEO FCDR generation 
    Author: Arta Dilo, NPL M&M
    Date created: 07-03-2016
    Last update: 28-10-2016 

FUNCTION declaration for visualising
    - ensemble instances' value and error images
    - distribution of ensemble values and errors at a pixel
    -  ...
    - measured & truth image and graphs of their values
    - ensemble instances (images and) graphs of values
Currently we consider images to have only one band. """


import matplotlib.mlab as mlab
from classFCDR import *


""" Plot value and error image for an ensemble instance. """
def mapInstance(ensemble, Idx, zran):
    # get ensemble instance
    inst = ensemble.getInstance(Idx)
    emeas = ensemble.getESource()

    # visualise image of instance values (meas - err)
    plt.figure() 
    plt.title('Ensemble instance %d' %Idx)
    imgV = plt.imshow(inst.getDNs(),cmap='BuGn',vmin=zran[0],vmax=zran[1]) 
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
    pxVals = pxDist[0,:] # ensemble values at (x,y)
    pxErrs = pxDist[1,:] # ensemble errors
    # get the measured value and uncertainty at pixel (x,y)
    pxMeas = ensemble.getESource().getPxVal(x, y)
    pxUnc = ensemble.getGumU().getPxVal(x, y) # GUM std. uncertainty at (x,y)
    pxEaM = pxMeas + pxErrs # ensemble alternate values at (x,y)
    # get the spread of ensemble values at pixel (x,y)
    pxVStd = np.std(pxVals) # st.dev. ensemble values at (x,y)
    pxDU = pxUnc * np.sqrt(2) # doubled (GUM) uncertainty at (x,y)
    
    """ Plot histogram of ensemble values at pixel (x,y)"""
    no_bin = 20 # number of bins for the histograms
    plt.figure()
    plt.hist(pxVals, no_bin, normed=1, fc='green', alpha=0.5)
    title = 'Histogram of values at pixel ' + str((x,y))
    plt.title(title)
    plt.xlabel('Ensemble value')
    plt.ylabel('Percentage') 
    plt.show()
    
    """ Plot histogram of ensemble errors at pixel (x,y) """
    plt.figure()
    plt.hist(pxErrs, no_bin, normed=1, fc='red', alpha=0.5)
    title = 'Histogram of errors at pixel ' + str((x,y))
    plt.title(title)
    plt.xlabel('Ensemble error')
    plt.ylabel('Percentage')
    plt.show()
    
    print '\nStandard deviation of ensemble values at pixel (',x,',',y,'):',pxVStd
    print 'GUM standard uncertainty at pixel (',x,',',y,'):', pxUnc
    
    """ Plot graph of ensemble values and errors at pixel (x,y) """
    plt.figure() 
    pv, = plt.plot(pxVals, 'go', label='Ens. values') # does not come in legend
    pav, = plt.plot(pxEaM, 'ro', alpha=0.4, label='Ens. alt. val') 
    # plot measured value and spread of ensemble values at (x,y)
    pm = plt.axhline(pxMeas, linewidth=3, color='c') # plot measured value
    pbnd = plt.axhspan(pxMeas - pxVStd, pxMeas + pxVStd, fc='c', alpha=0.4)
    pUnc = plt.axhline(pxMeas - pxUnc, linewidth=2, color='b') 
    plt.axhline(pxMeas + pxUnc, linewidth=2, color='b') 
    peUnc = plt.axhline(pxMeas - pxDU, linewidth=1, color='m') 
    plt.axhline(pxMeas + pxDU, linewidth=1, color='m') 
    # add legend to the graph
    pvLbl = ['Ensemble vals','Alternate vals','Measured val','+/- Ens. St.Dev.','St.uncertainty']
    plt.legend([pv, pav, pm, pbnd,pUnc], pvLbl)
    title = 'Ensemble values & alternates at pixel ' + str((x,y))
    plt.title(title)
    plt.xlabel('Instance number')
    plt.ylabel('Instance (alternate) values')
    plt.show()
    

""" Plot images, scatterplot and histogram of measured (and truth) values. 
The graphs, scatterplot & histogram, are meaningful (i.e. useful in some way) 
only for the case when we start from a flat truth. Others should be thought of 
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
    print '\n(Sample) mean of measured values:', meanM 
    print '(Sample) standard deviation of measured:', stdM
    print 'Standard uncertainty of measured image:', mU
    
    """ Plot the (simulated) truth and measured image """
    zran = [minM, maxM]
    truth.imgMap('Assumed truth image','Greens',zran)    
    measured.imgMap('Measured generated from flat truth','Greens',zran) 

    """ Plot graph of measured values and stats """
    plt.figure() 
    mv, = plt.plot(mVals, 'g+')#, label='Measured')
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
    mvLbl = ['Measured','Mean measured','+/-Standard dev.','Truth','+/- Std uncert.']
    plt.legend([mv,sm,sbnd,at,gumub], mvLbl)
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
    histW = (maxp-minp)/max(10,np.ceil(3*len(insVals)**0.3))
    hRng = np.arange(minp, maxp, histW)

    meanE2 = np.mean(insAVals) # mean & stdev of instances (calc: meas+err)
    stdE2 = np.std(insAVals) 
    
    # Get ensemble uncertainty (image), max & min uncertainty value
    uncert = ensemble.getGumU() # GUM uncertainty image
    eU = max(uncert.getValArr()) # max uncertainty value
    emU = np.mean(uncert.getValArr()) # average uncertainty value
    
    esize = ensemble.getESize() # ensemble size
    print '\nPicked from the ensemble of', esize, 'instances:', INo
    print '(Sample) mean of picked instances (meas-err, meas+err):', meanpE,',',meanE2
    print '(Sample) st.dev of instances (meas-err, meas+err):', stdpE,',',stdE2
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
    spLbl = ['Ensemble values','Ensemble Mean','+/- Ensemble Std','Mean measured','+/- max std uncert.']
    plt.legend([ev,em,sibnd,mm,gumeub],spLbl)
    title = 'Pixel values of ensemble instances: ' + str(INo)
    plt.title(title)
    plt.xlabel('Pixel count')
    plt.ylabel('Pixel value')
    plt.show()
    
    """ Histograms of measured & ensemble instances (meas+errs, meas-errs) """
    # histogram to compare distribution of measured with ensemble values
    plt.figure()
    n, bins, patches = plt.hist(mVals,hRng,normed=1,fc='blue',alpha=0.5)
    ym = mlab.normpdf(bins, meanM, stdM) # gaussian 'best fit' for measured
    plt.plot(bins, ym, 'm--')
    en, ebin, epat = plt.hist(insVals,hRng,normed=1,fc='green',alpha=0.5)#,rwidth=histW)
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