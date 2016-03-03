#/* *************************************************************************** */
#/* *************************************************************************** */
#/* THESE ARE FEW EXAMPLES OF TYPICAL APPLICATIONS OF THE TOOLKIT               */
#/*                                                                             */
#/* They have been designed mainly so as to ilustrate to you how to roll your   */
#/* own applications.                                                           */
#/* Select an example that seems closest to your goal and edit it.              */
#/* In most cases you will find the tools you need used in the examples.        */
#/* In some cases you will need to write some code similar to what you find in  */
#/* the implementation of the functions under the tab "Disarray", whereas       */
#/* various low level functions are available under the tab "Utilities".        */
#/* In rare cases you may need to write raw code.                               */
#/* There are examples of any kind.                                             */
#/* *************************************************************************** */


#==============================================================================
# Imports
#==============================================================================
from PIL import Image
import numpy as np
import time

#==============================================================================
# Eidolon imports
#==============================================================================
from picture import *
from scalespaces import *
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

def StackIntegrateMultiple(aScaleSpace, numScaleLevels, scaleLevels, picSize):
    first = list()
    second = list()
    tmp = list()

    first = aScaleSpace.next()
    numberOfScaleSpaces = len(first)
    
    for i in range (numberOfScaleSpaces - 1):
        tmp.append(np.zeros(picSize))
        
    for k in range (numScaleLevels-1):
        interval = 0.5 * (scaleLevels[k] - scaleLevels[k+1])      
        second = aScaleSpace.next()
        for i in range (numberOfScaleSpaces - 1):
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


    
#==============================================================================
# Functions
#==============================================================================

#* *************************************************************************** */
#/* *************************************************************************** */
#/* SYNTHESIS FROM D.O.G. ACTIVITY: A SANITY CHECK                              */
#/*                                                                             */
#/* This may seem like an unremarkable example: the output just the image       */
#/* again!), but it actually did decompose the image into "edges" and from      */
#/* these put it together again. Some feat! In producing eidolons you will      */
#/* simply dislocate the parts produced by the decomposition before attempting  */
#/* the synthesis. This is the generic way to produce eidolons.                 */
#/* *************************************************************************** */
def synthesisFromDOGActivity(pic):
    print "Embarking on synthesis from DOG activity"
    start_time =  time.clock()

    #// B. SET UP NECESSARY DATA STRUCTURES
    fiducialDOGScaleSpace = DOGScaleSpace(pic)     #// this is the datastructure that is really needed
  
    print "Computing eidolon..."
      
    #// C. CONSTRUCT THE EIDOLON
    eidolonDataPlane = np.zeros(pic.fatFiducialDataPlane.shape)
    #// here the eidolon is constructed
    for i in range(pic.numScaleLevels):
        eidolonDataPlane += fiducialDOGScaleSpace.next() # this is a generator!
    eidolon = pic.DisembedDataPlane(eidolonDataPlane) #// convert dataplane to image
    
    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    Image.fromarray(eidolon, 'L').show()


#/* *************************************************************************** */
#/* *************************************************************************** */
#/* SYNTHESIS FROM THE LAPLACIAN: A SANITY CHECK                         */
#/*                                                                             */
#/* This may seem like an unremarkable example: the output just the image       */
#/* again!).         */
#/* *************************************************************************** */
def synthesisFromLaplacian(pic, INTEGRATION_FUDGE_FACTOR):
    print "Embarking on synthesis from Laplacians"
    start_time =  time.clock()

    #// B. SET UP NECESSARY DATA STRUCTURES
    rockBottomPlaneGenerator = RockBottomPlane(pic)       # // needed because of the lowest resolution plane
    fiducialLaplacianScaleSpace = FiducialLaplacian(pic)          #// sets up the Laplacian activity
    numScaleLevels = pic.numScaleLevels
    scaleLevels = pic.scaleLevels
    picSize = pic.fatFiducialDataPlane.shape
    
    print "Computing eidolon..."
      
    #// C. CONSTRUCT THE EIDOLON
    eidolonDataPlane = np.zeros(pic.fatFiducialDataPlane.shape)
    #// here the eidolon is constructed
    eidolonDataPlane = rockBottomPlaneGenerator.next() + StackIntegrate(fiducialLaplacianScaleSpace, numScaleLevels, scaleLevels, picSize) * INTEGRATION_FUDGE_FACTOR
    
    eidolon = pic.DisembedDataPlane(eidolonDataPlane) #// convert dataplane to image
    
    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    Image.fromarray(eidolon, 'L').show()
    

#/* *************************************************************************** */
#/* *************************************************************************** */
#/* SYNTHESIS FROM SIMPLE CELL ACTIVITY: A SANITY CHECK                         */
#/*                                                                             */
#/* This may seem like an unremarkable example: the output just the image       */
#/* again!), but it actually did decompose the image into "edges" of all        */
#/* possible resolution and from these put it together again. Some feat! In     */
#/* producing eidolons you will simply dislocate the parts produced by the      */
#/* decomposition before attempting the synthesis. This is the generic way to   */
#/* produce eidolons. It can be simplified by dislocating various coarser       */
#/* decompositions like the DOG-representation (no orientations) or even the    */
#/* pixel representation (no multiple scales). These simpler methods may        */
#/* suffice for most of your needs, so do not ignore them.                      */
#/* *************************************************************************** */
def synthesisFromSimpleCellActivity(pic, INTEGRATION_FUDGE_FACTOR):
    THREE_FOURTH = 3.0/4
    print "Embarking on synthesis from simple cell activity"
    start_time =  time.clock()

    #// B. SET UP NECESSARY DATA STRUCTURES
    rockBottomPlaneGenerator = RockBottomPlane(pic)       # // needed because of the lowest resolution plane
    fiducialSecondOrderGenerator = FiducialSecondOrder(pic)          # // sets up the simple cell activity
    numScaleLevels = pic.numScaleLevels
    scaleLevels = pic.scaleLevels
    picSize = pic.fatFiducialDataPlane.shape

    print "Computing eidolon..."
      
    #// C. CONSTRUCT THE EIDOLON
    integratedList = StackIntegrateMultiple(fiducialSecondOrderGenerator, numScaleLevels, scaleLevels, picSize)
    eidolonDataPlane = np.zeros(pic.fatFiducialDataPlane.shape)
    # make sum of all the stacks
    for il in integratedList:
        eidolonDataPlane += il

    eidolonDataPlane = rockBottomPlaneGenerator.next() + eidolonDataPlane * INTEGRATION_FUDGE_FACTOR * THREE_FOURTH
 
    eidolon = pic.DisembedDataPlane(eidolonDataPlane) #// convert dataplane to image
    
    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    Image.fromarray(eidolon, 'L').show()


#* *************************************************************************** */
#/* *************************************************************************** */
#/* SUPERFICIAL DISARRAY ("PIXEL SHIFTING")                                     */
#/*                                                                             */
#/* This implements "CAPTCHA"-style rendering.                                  */
#/* *************************************************************************** */
def SuperficialDisarray(pic, reach,  grain):                    #// parameters defined in function definition
    print "Embarking on superficial disarray"
    start_time =  time.clock()

    #// B. SET UP NECESSARY DATA STRUCTURES  - NOT NEEDED HERE

    print "Computing eidolon..."
      
    #// C. CONSTRUCT THE EIDOLON  
    (h, w) = pic.fatFiducialDataPlane.shape
    b1 = BlurredRandomGaussianDataPlane(w, h, grain)
    xDisplacements = b1.blurredRandomGaussianDataPlane 
#    # normalize and *5 (the 5 is a magic number, to match Jan's version)
#    xDisplacements = (xDisplacements / np.max(xDisplacements)) * 5
    b2 = BlurredRandomGaussianDataPlane(w, h, grain)
    yDisplacements = b2.blurredRandomGaussianDataPlane 
#    # normalize and *5 (the 5 is a magic number, to match Jan's version)
#    yDisplacements = (yDisplacements / np.max(yDisplacements)) * 5
      
    eidolonDataPlane = DataPlaneDisarray(pic.fatFiducialDataPlane, xDisplacements, yDisplacements, reach)
  
    eidolon = pic.DisembedDataPlane(eidolonDataPlane) #// convert dataplane to image
    
    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    Image.fromarray(eidolon, 'L').show()


#/* *************************************************************************** */
#/* LOTZE TYPE DISARRAY                                                         */
#/*                                                                             */
#/* This implements "Lotze-type" disarray of the edges (DOG) representation.    */
#/* All scale levels are disarryed by statistically identical (though           */
#/* independent) perturbation.                                                  */
#/* *************************************************************************** */
def LotzeTypeDisarray(pic, reach, grain):                   #  // parameter defined in function definition
    print "Embarking on Lotze-type disarray"
    start_time =  time.clock()
    
    (h,w) = pic.fatFiducialDataPlane.shape
    numScaleLevels = pic.numScaleLevels
    
    #// B. SET UP NECESSARY DATA STRUCTURES
    fiducialDOGScaleSpace = DOGScaleSpace(pic)     #// this is the datastructure that is really needed
  
    print "Computing eidolon..."
      
    #// C. CONSTRUCT THE EIDOLON
    eidolonDataPlane = np.zeros(pic.fatFiducialDataPlane.shape)
    
    #// here the eidolon is constructed  
    eidolonDataPlane = lotzeDisarray(fiducialDOGScaleSpace, reach, grain, numScaleLevels, w, h)

    eidolon = pic.DisembedDataPlane(eidolonDataPlane) #// convert dataplane to image
    
    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    Image.fromarray(eidolon, 'L').show()


#// All scalespace layers are individually dislocated by independent but statistically equal displacement fields
def lotzeDisarray(aDOGScaleSpace, reach, grain, numScaleLevels, w, h):
    tmp = np.zeros((h,w))
    i1 = IncoherentGaussianDataStack(numScaleLevels, w, h, grain) # xDisplacements 
    i2 = IncoherentGaussianDataStack(numScaleLevels, w, h, grain) # yDisplacements
    #// synthesizes with disarrayed DOG scale space layers
    #// this essentially yields an exact integral over scale 
    #// because the DOG samples are slices, not poit samples
    for i in range (numScaleLevels):
        tmp += DataPlaneDisarray(aDOGScaleSpace.next(), i1.next(), i2.next(), reach)    
    return tmp








#==============================================================================
# 
# Program
# 
#==============================================================================
def testFunction(): 
    SZ = 512
    MIN_SIGMA = 1/sqrt(2) #// 1f/sqrt(2f)
    MAX_SIGMA = SZ/4.0
    SIGMA_FACTOR = sqrt(2)  #// sqrt(2f)

    INTEGRATION_FUDGE_FACTOR = 1.25; #// Serves to improves the stack integration by the trapezoidal rule 
                                     #// The scale intervals are unequal in width and depend systematically on scale
                                     #// In a future update we may update the numerical integration routine
                                     #// Expected difference are slight, present method works quite well


    # parameters
    REACH = 2.0 #// These are just starting values - adjust as you see fit (try factors of 2 from here for a start,
                #// in some cases you need REALLY large changes before noticing a difference).
    GRAIN = 8.0  #// Not all of these parameters are necessarily relevant for a given method (see above).
    LEVEL = 6  #// These numbers are not of any essence:
    FRACTION = 0.90  #// Once you get the hang of it, you won't need the precooked stuff anyway!    
    INCOHERENCE = 0.10  #// It is imperative to obtain an intuitive feel for the effects of these parameters (and their interactions) though.


    # make a picture class to be used in the examples
#    pic = Picture('Baby_face.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)
#    pic = Picture('berries.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)
    pic = Picture('DenisHopper.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)
#    pic = Picture('Hanna.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)
#    pic = Picture('Hanna2.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)
#    pic = Picture('Hanna_Salience.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)
#    pic = Picture('Mona_Lisa.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)
#    pic = Picture('Rudy_Dekeerschieter.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)
#    pic = Picture('Rudy_Dekeerschieter_ZW.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)

#    synthesisFromDOGActivity(pic)    
#    synthesisFromLaplacian(pic, INTEGRATION_FUDGE_FACTOR)    
#    synthesisFromSimpleCellActivity(pic, INTEGRATION_FUDGE_FACTOR)    
#    SuperficialDisarray(pic, REACH, GRAIN)     
    LotzeTypeDisarray(pic, REACH, GRAIN)
    
    
    
    
    
    
    
#==============================================================================
# main
#==============================================================================
if __name__ == "__main__":
    testFunction()
    print "Done!"