""" FIDUCEO FCDR generation 
    Author: Arta Dilo, NPL M&M
    Date created: 28-10-2016
    Last update: 31-10-2016 

FUNCTION declaration for evaluation of
    - Earth radiance from counts and calibration data i.e. calibration model
    - radiance (L) from temperature: Planck's law
    - temperature from radiance conversion, Planck's inverse function
    - radiance uncertainty via GUM (partial derivatives and vars uncertainty)
    - ...
Some physics constants e.g. 1st and 2nd radiantion contant, and global constants
/variables are declare here. """

import numpy as np

""" Declare global constants and variables """
c1 = 1.191062e-5    # 1st radiation constant
c2 = 1.4387863     # 2nd radiation constant

a = (1.16, 0, 8e-6) # calibration coefficients; will be output of harmonisation
ch11 = (918.956, 0.613784, 0.998183) # channel 11 central freq, aval and val

uCnt = 0.5      # counts standard uncertainty ?
uTmp = 0.008    # ICT temperature std. uncertainty


""" Convert radiance from temperature and central band frequency: 
Planck's law"""
def radPL(temp, freq, aval, bval):
    radFt = (c1 * np.power(freq,3)) 
    radFt = radFt / (np.exp(c2 * freq/(aval + bval * temp)) - 1)
    return(radFt)

""" Convert temperature from radiance: inverse Planck's law """
def tempIPL(rad, freq, aval, bval):
    et = c2 * freq / np.log(c1 * np.power(freq,3)/rad + 1) # effective temp.
    temp = (et - aval) / bval                 # temperature
    return(temp)


""" Evaluate radiance uncertainty using GUM law:
currently calculated only with temperature uncertainty from the partial 
derivative of Planck's law on temperature, i.e. assumes uncertainties from 
band frequency and aval & bval variables are 0. Check this again! """
def Lunc(temp, freq, aval, bval, utemp):
    et = aval + bval * temp
    pdev = (c1*c2 * np.power(freq,4) * bval) 
    pdev = pdev / (np.power(et,2) * (np.exp(c2 * freq /et) - 1))
    unc = utemp * pdev
    return(unc)


""" Evaluate Earth radiance from Earth counts and calibration data, this is a  
simplified calibration model. For this test define some values for model 
coefficients. The coeficients will be estimated from the harmonisation. """
def LEarth(CE, Cict, Tict, acoef, chp):
    Lict = radPL(Tict, chp[0], chp[1], chp[2])
    LE = acoef[0] + CE * (Lict / Cict + acoef[1]) + acoef[2] * np.power(CE, 2)
    return LE
    
#CE = 350
#Tict = 295
#Cict = 580
#LE = LEarth(CE, Cict, Tict, a, ch11)
#print 'Earth radiance ', LE