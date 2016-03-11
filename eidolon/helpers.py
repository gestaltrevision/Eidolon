#==============================================================================
# Imports
#==============================================================================
import numpy as np

#==============================================================================
# Eidolon imports
#==============================================================================
from noise import *

#==============================================================================
# Helper functions
#==============================================================================
#// This replaces the integral over scale by a finite sum
#// as a complication, the samples are not equispaced over the scale domain
#// Here the trapezoid rule is used
#// In order to obtain the integral one needs to add the lowest resolution image ("rockbottomPlane")
def StackIntegrate(aScaleSpace, numScaleLevels, scaleLevels, picSize):
    tmp = np.zeros(picSize)

    first = aScaleSpace.next()
    for k in range (numScaleLevels-1):
        interval = 0.5 * (scaleLevels[k] - scaleLevels[k+1])
        second = aScaleSpace.next()
        tmp = tmp + (first + second) * (interval * scaleLevels[k])
        first = second
        
    return tmp

    
# Same as StackIntegrate, only for multiple scalespaces at once
def StackIntegrateMultiple(aScaleSpace, numScaleLevels, scaleLevels, picSize):
    first = list()
    second = list()
    tmp = list()

    first = aScaleSpace.next()
    numberOfScaleSpaces = len(first)
    
    for i in range (numberOfScaleSpaces):
        tmp.append(np.zeros(picSize))
        
    for k in range (numScaleLevels-1):
        interval = 0.5 * (scaleLevels[k] - scaleLevels[k+1])      
        second = aScaleSpace.next()
        for i in range (numberOfScaleSpaces):
            tmp[i] = tmp[i] + (first[i] + second[i]) * (interval * scaleLevels[k])
        first = second        

    return tmp
    
    
#// Superficial pixel-shifting: no scalespaces needed! Any dataplane can be used
# xDisplacements and yDisplacements are matrices same size as image, reach is a number
def DataPlaneDisarray(dataPlane, xDisplacements, yDisplacements, reach):
    h, w = dataPlane.shape[0], dataPlane.shape[1]
    # grid contains x coordinates in layer 0, y coordinates in layer 1
    grid = np.indices((h,w))
    # shuffles the dataPlane     
    xNew = (np.clip((xDisplacements * reach + grid[0]), 0, h-1)).astype(int)
    yNew = (np.clip((yDisplacements * reach + grid[1]), 0, w-1)).astype(int)
    
    return dataPlane[xNew, yNew]


def QuartilesAndPercentLevels(data, theLevel = 0.5): 
    if theLevel < 0 or theLevel >= 1:
        raise ValueError('Quantile has to be between 0 (inclusive) and 1 (exclusive)!')
        
    tmp = np.sort(data.reshape(data.size))
    q = {
            'q0.01' : tmp[np.round(0.01*tmp.size)],
            'q0.05' : tmp[np.round(0.05*tmp.size)],
            'q0.25' : tmp[np.round(0.25*tmp.size)],
            'q0.50' : tmp[np.round(0.50*tmp.size)],
            'q0.75' : tmp[np.round(0.75*tmp.size)],
            'q0.95' : tmp[np.round(0.95*tmp.size)],
            'q0.99' : tmp[np.round(0.93*tmp.size)],
            'qLevel' : tmp[np.round(theLevel*tmp.size)]
        }
    return q

    
def DataToImage(data):
    quantiles = QuartilesAndPercentLevels(data)
    mn = (quantiles['q0.95'] + quantiles['q0.05']) / 2.0
    amp = (quantiles['q0.95'] - quantiles['q0.05']) / 2.0
    data = 127.5 + 127.5 * (data-mn) / amp
    data[data < 0] = 0
    data[data > 255] = 255
    return data.astype('uint8')


#// All scalespace layers are individually dislocated by independent but statistically equal displacement fields
def LotzeDisarray(aDOGScaleSpace, reach, grain, numScaleLevels, w, h):
    tmp = np.zeros((h,w))
    i1 = IncoherentGaussianDataStack(numScaleLevels, w, h, grain) # xDisplacements 
    i2 = IncoherentGaussianDataStack(numScaleLevels, w, h, grain) # yDisplacements
    #// synthesizes with disarrayed DOG scale space layers
    #// this essentially yields an exact integral over scale 
    #// because the DOG samples are slices, not poit samples
    for i in range (numScaleLevels):
        tmp += DataPlaneDisarray(aDOGScaleSpace.next(), i1.next(), i2.next(), reach)    
    return tmp


#// All scalespace layers are individually dislocated by independent displacement fields, scaled with the local resolution
def HelmholtzDisarray(aDOGScaleSpace, reach, numScaleLevels, w, h, MAX_SIGMA, scaleLevels):
    tmp = np.zeros((h,w))
    i1 = IncoherentScaledGaussianDataStack(numScaleLevels, w, h, MAX_SIGMA, scaleLevels)
    i2 = IncoherentScaledGaussianDataStack(numScaleLevels, w, h, MAX_SIGMA, scaleLevels)

    for i in range (numScaleLevels):
        tmp += DataPlaneDisarray(aDOGScaleSpace.next(), i1.next(), i2.next(), MAX_SIGMA * reach)    
    return tmp
    
    
#// All scalespace layers are individually dislocated by mutually dependent displacement fields, scaled with the local resolution
#// Coarse resolution layers contribute to fine resolution layers, thus large RFs drag small RFs along with them
def CoherentDisarray(aDOGScaleSpace, reach, w, h, MAX_SIGMA, numScaleLevels, scaleLevels):
    tmp = np.zeros((h,w))
    c1 = CoherentRandomGaussianDataStack(numScaleLevels, w, h, MAX_SIGMA, scaleLevels)
    c2 = CoherentRandomGaussianDataStack(numScaleLevels, w, h, MAX_SIGMA, scaleLevels)

    for i in range (numScaleLevels):
        tmp += DataPlaneDisarray(aDOGScaleSpace.next(), c1.next(), c2.next(), reach)    
    return tmp    
    

#// this replaces the integral over scales by a finite sum
#// as a complication, the samples are not equispaced over the scale domain
#// Here the trapezoid rule is used
def StackDisarrayDiscrete(aScaleSpace, xP, yP, xQ, yQ, xR, yR, reach, numScaleLevels, scaleLevels, picSize):  
    sdd = StackDisarrayDiscreteGenerator(aScaleSpace, xP, yP, xQ, yQ, xR, yR, reach)  
    return StackIntegrateMultiple(sdd, numScaleLevels, scaleLevels, picSize)  


# helper for StackDisarrayDiscrete  
def StackDisarrayDiscreteGenerator(aScaleSpace, xP, yP, xQ, yQ, xR, yR, reach):
    for p, q, r in aScaleSpace:
        pD = DataPlaneDisarray(p, xP.next(), yP.next(), reach)
        qD = DataPlaneDisarray(q, xQ.next(), yQ.next(), reach)
        rD = DataPlaneDisarray(r, xR.next(), yR.next(), reach)
        yield (pD, qD, rD)


def SuppressSmallActivity(dataPlane, fraction):
    if fraction < 0 or fraction >= 1:
        raise ValueError('Fraction has to be between 0 (inclusive) and 1 (exclusive)!')
    
    for data in dataPlane:
        tmp = np.sort(data.reshape(data.size))   
        threshold = tmp[np.round(fraction*tmp.size)]        
        data[data < threshold] = 0
        yield data


#// Planar vector image - hue drom direction, saturation signifies magnitude
def VectorImage(x, y):   
    eidolonDataPlane = np.ones((x.shape[0], x.shape[1], 3)) * 255    
    rho = x**2 + y**2
    phi = np.arctan2(y,x)    
    maxRho = np.max(rho)

    limit = 0.9 #// the maximum saturation
    cnst = maxRho * (1.0-limit)/limit

    R, G, B = DirectionColor(phi)
    f = rho/(cnst + rho)

    eidolonDataPlane[:,:,0] = LerpColor(eidolonDataPlane[:,:,0], R, f)
    eidolonDataPlane[:,:,1] = LerpColor(eidolonDataPlane[:,:,1], G, f)
    eidolonDataPlane[:,:,2] = LerpColor(eidolonDataPlane[:,:,2], B, f)
    
    return eidolonDataPlane


def DirectionColor(direction):
    direction = (1 + direction/np.pi) * 3.0  # // From (-PI,+PI) to (0,6)
    i = np.floor(direction);
    f = direction - i;
    resultR = np.ones(direction.shape) * 127
    resultG = np.ones(direction.shape) * 127
    resultB = np.ones(direction.shape) * 127

    # YELLOW to GREEN
    resultR[i==0] = LerpColor(255, 0, f[i==0])
    resultG[i==0] = 255
    resultB[i==0] = 0
    # GREEN to CYAN
    resultR[i==1] = 0
    resultG[i==1] = 255
    resultB[i==1] = LerpColor(0, 255, f[i==1])
    # CYAN to BLUE
    resultR[i==2] = 0
    resultG[i==2] = LerpColor(255, 0, f[i==2])
    resultB[i==2] = 255
    # BLUE to MAGENTA
    resultR[i==3] = LerpColor(0, 255, f[i==3])
    resultG[i==3] = 0
    resultB[i==3] = 255
    # MAGENTA to RED
    resultR[i==4] = 255
    resultG[i==4] = 0
    resultB[i==4] = LerpColor(255, 0, f[i==4])
    # RED to YELLOW
    resultR[i==5] = 255
    resultG[i==5] = LerpColor(0, 255, f[i==5])
    resultB[i==5] = 0
    
    return resultR, resultG, resultB

def LerpColor(a, b, f):
    # this basically returns a matrix with values between the values of 
    # a and b, depending on the value of f (f=0 returns a, f=1 returns b )
    return np.floor((b-a)*f + a)


#// Shows the 3D line finder vector basis in RGB
def LineFinderBasisImage(secondOrderP, secondOrderQ, secondOrderR):
    eidolonDataPlane = np.zeros((secondOrderP.shape[0], secondOrderP.shape[1], 3))
    
    qpp = QuartilesAndPercentLevels(secondOrderP, 0.75)
    qpq = QuartilesAndPercentLevels(secondOrderQ, 0.75)
    qpr = QuartilesAndPercentLevels(secondOrderR, 0.75)
    factor = 127.0 / max(qpp['qLevel'], qpq['qLevel'], qpr['qLevel'])

    eidolonDataPlane[:,:,0] = 127 + factor * secondOrderP
    eidolonDataPlane[:,:,1] = 127 + factor * secondOrderQ
    eidolonDataPlane[:,:,2] = 127 + factor * secondOrderR
    
    return eidolonDataPlane


def PartiallyCoherentDisarray(aDOGScaleSpace, reach, degree, sigma, w, h, MAX_SIGMA, numScaleLevels, scaleLevels):
    xDisplacements = PartiallyCoherentScaledGaussianDataStack(numScaleLevels, w, h, sigma, MAX_SIGMA, scaleLevels, degree)
    yDisplacements = PartiallyCoherentScaledGaussianDataStack(numScaleLevels, w, h, sigma, MAX_SIGMA, scaleLevels, degree)
    return StackDisarray(aDOGScaleSpace, xDisplacements, yDisplacements, reach)


#// synthesizes with disarrayed DOG scale space layers
#// this essentially yields an exact integral over scale 
#// because the DOG samples are slices, not poit samples
def StackDisarray(aDOGScaleSpace, xDisplacements, yDisplacements, reach):
    sd = StackDisarrayGenerator(aDOGScaleSpace, xDisplacements, yDisplacements, reach)
    tmp = None  
    for plane in sd:
        if tmp == None:
            tmp = np.zeros(plane.shape)
        tmp += plane
    return tmp


# helper for StackDisarray
def StackDisarrayGenerator(aScaleSpace, x, y, reach):
    for scaleSpace in aScaleSpace:
        yield DataPlaneDisarray(scaleSpace, x.next(), y.next(), reach)


def ImageToOpponentRepresentation(dataRed, dataGreen, dataBlue):
    # be sure that matrices are of type to float otherwise calculations go horribly wrong!
    kw = dataRed + dataGreen + dataBlue
    rg = dataRed - dataGreen 
    yb = dataRed + dataGreen - 2.0 * dataBlue
    return kw, rg, yb

def OpponentRepresentationToImage(kw, rg, yb):
    # be sure that matrices are of type to float otherwise calculations go horribly wrong!
    r = np.round((2.0 * kw + 3.0 * rg + yb) / 6.0)
    g = np.round((2.0 * kw - 3.0 * rg + yb) / 6.0)
    b = np.round((kw - yb) / 3.0) 
    return r, g, b


def ITORGeneratorKW(scaleSpaceR, scaleSpaceG, scaleSpaceB):
    i = 0
    for sp in scaleSpaceR:
        kw, rg, yb = ImageToOpponentRepresentation(sp, scaleSpaceG[i], scaleSpaceB[i])
        i += 1
        yield kw   
        
def ITORGeneratorRG(scaleSpaceR, scaleSpaceG, scaleSpaceB):
    i = 0
    for sp in scaleSpaceR:
        kw, rg, yb = ImageToOpponentRepresentation(sp, scaleSpaceG[i], scaleSpaceB[i])
        i += 1
        yield rg   
        
def ITORGeneratorYB(scaleSpaceR, scaleSpaceG, scaleSpaceB):
    i = 0
    for sp in scaleSpaceR:
        kw, rg, yb = ImageToOpponentRepresentation(sp, scaleSpaceG[i], scaleSpaceB[i])
        i += 1
        yield yb

