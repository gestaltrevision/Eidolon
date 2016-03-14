#/* *************************************************************************** */
#/* *************************************************************************** */
#/* *                                                                         * */
#/* *                           EIDOLON TOOLKIT                               * */
#/* *                                                                         * */
#/* *             de Clootcrans - Jan Koenderink 28 nov 2015                  * */
#/* *                                                                         * */
#/* *                             version 0.0                                 * */
#/* *                                                                         * */
#/* *                                                                         * */
#/* *  please report any bugs to me at KoenderinkJan@gmail.com mentioning     * */
#/* * "Eidolon toolkit" in the subject slot of your message.                  * */
#/* *                                                                         * */
#/* * I can accept no responsibility for your disasters, sorry about that!    * */
#/* *                                                                         * */
#/* *************************************************************************** */
#/* *************************************************************************** */
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
from eidolon.picture import *
from eidolon.scalespaces import *
from eidolon.noise import *
from eidolon.helpers import *

#==============================================================================
# Functions
#==============================================================================

#* **************************************************************************** */
#/* *************************************************************************** */
#/* SYNTHESIS FROM D.O.G. ACTIVITY: A SANITY CHECK                              */
#/*                                                                             */
#/* This may seem like an unremarkable example: the output just the image       */
#/* again!), but it actually did decompose the image into "edges" and from      */
#/* these put it together again. Some feat! In producing eidolons you will      */
#/* simply dislocate the parts produced by the decomposition before attempting  */
#/* the synthesis. This is the generic way to produce eidolons.                 */
#/* *************************************************************************** */
def SynthesisFromDOGActivity(pic):
    print "Embarking on synthesis from DOG activity"
    start_time = time.clock()

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
    
    return Image.fromarray(eidolon, 'L')


#/* *************************************************************************** */
#/* *************************************************************************** */
#/* SYNTHESIS FROM THE LAPLACIAN: A SANITY CHECK                                */
#/*                                                                             */
#/* This may seem like an unremarkable example: the output just the image       */
#/* again!).                                                                    */
#/* *************************************************************************** */
def SynthesisFromLaplacian(pic, INTEGRATION_FUDGE_FACTOR):
    print "Embarking on synthesis from Laplacians"
    start_time = time.clock()

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
    
    return Image.fromarray(eidolon, 'L')
    

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
def SynthesisFromSimpleCellActivity(pic, INTEGRATION_FUDGE_FACTOR):
    THREE_FOURTH = 3.0/4
    print "Embarking on synthesis from simple cell activity"
    start_time = time.clock()

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
    
    return Image.fromarray(eidolon, 'L')


#/* *************************************************************************** */
#/* *************************************************************************** */
#/* SUPERFICIAL DISARRAY ("PIXEL SHIFTING")                                     */
#/*                                                                             */
#/* This implements "CAPTCHA"-style rendering.                                  */
#/* *************************************************************************** */
def SuperficialDisarray(pic, reach,  grain):                                  #// parameters defined in function definition
    print "Embarking on superficial disarray"
    start_time = time.clock()

    #// B. SET UP NECESSARY DATA STRUCTURES  - NOT NEEDED HERE

    print "Computing eidolon..."
      
    #// C. CONSTRUCT THE EIDOLON  
    (h, w) = pic.fatFiducialDataPlane.shape
    b1 = BlurredRandomGaussianDataPlane(w, h, grain)
    xDisplacements = b1.blurredRandomGaussianDataPlane 
    b2 = BlurredRandomGaussianDataPlane(w, h, grain)
    yDisplacements = b2.blurredRandomGaussianDataPlane 
      
    eidolonDataPlane = DataPlaneDisarray(pic.fatFiducialDataPlane, xDisplacements, yDisplacements, reach)
  
    eidolon = pic.DisembedDataPlane(eidolonDataPlane) #// convert dataplane to image
    
    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    return Image.fromarray(eidolon, 'L')


#/* *************************************************************************** */
#/* *************************************************************************** */
#/* LOTZE TYPE DISARRAY                                                         */
#/*                                                                             */
#/* This implements "Lotze-type" disarray of the edges (DOG) representation.    */
#/* All scale levels are disarryed by statistically identical (though           */
#/* independent) perturbation.                                                  */
#/* *************************************************************************** */
def LotzeTypeDisarray(pic, reach, grain):                                     # // parameter defined in function definition
    print "Embarking on Lotze-type disarray"
    start_time = time.clock()
    
    (h,w) = pic.fatFiducialDataPlane.shape
    numScaleLevels = pic.numScaleLevels
    
    #// B. SET UP NECESSARY DATA STRUCTURES
    fiducialDOGScaleSpace = DOGScaleSpace(pic)     #// this is the datastructure that is really needed
  
    print "Computing eidolon..."
      
    #// C. CONSTRUCT THE EIDOLON
    eidolonDataPlane = np.zeros(pic.fatFiducialDataPlane.shape)
    
    #// here the eidolon is constructed  
    eidolonDataPlane = LotzeDisarray(fiducialDOGScaleSpace, reach, grain, numScaleLevels, w, h)

#    eidolon = pic.DisembedDataPlane(eidolonDataPlane) #// convert dataplane to image
    eidolon = DataToImage(pic.DisembedDataPlane(eidolonDataPlane, clip=False)) #// convert dataplane to image
    
    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    return Image.fromarray(eidolon.astype('uint8'), 'L')


#/* *************************************************************************** */
#/* *************************************************************************** */
#/* HELMHOLTZ TYPE DISARRAY                                                     */
#/*                                                                             */
#/* This implements "Helmholtz-type" disarray of the edges (DOG)                */
#/* representation.  All scale levels are disarryed by statistically identical  */
#/* (though independent) perturbation, scaled by the local resolution.          */
#/* *************************************************************************** */
def HelmholtzTypeDisarray(pic, reach):  
    MAX_SIGMA = pic.MAX_SIGMA                           
    print "Embarking on Helmholtz-type disarray"
    start_time = time.clock()
    
    (h,w) = pic.fatFiducialDataPlane.shape
    numScaleLevels = pic.numScaleLevels
    scaleLevels = pic.scaleLevels
    
    #// B. SET UP NECESSARY DATA STRUCTURES
    fiducialDOGScaleSpace = DOGScaleSpace(pic)     #// this is the datastructure that is really needed
  
    print "Computing eidolon..."
      
    #// C. CONSTRUCT THE EIDOLON
    eidolonDataPlane = np.zeros(pic.fatFiducialDataPlane.shape)
    
    #// here the eidolon is constructed  
    eidolonDataPlane = HelmholtzDisarray(fiducialDOGScaleSpace, reach, numScaleLevels, w, h, MAX_SIGMA, scaleLevels)
    
    eidolon = pic.DisembedDataPlane(eidolonDataPlane) #// convert dataplane to image
    
    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    return Image.fromarray(eidolon, 'L')


#/* *************************************************************************** */
#/* *************************************************************************** */
#/* COHERENT DISARRAY OF EDGES                                                  */
#/*                                                                             */
#/* Here the edge representation (DOG activity is disarrayed on all scales      */
#/* levels. The reach at a level is proprtional to the resolution. Moreover,    */
#/* large RFs drag smaller ones along with them.                                */
#/* *************************************************************************** */
def CoherentDisarrayOfEdges(pic, reach):                         # // parameter "reach" in function definition
    MAX_SIGMA = pic.MAX_SIGMA       
    print "Embarking on coherent disarray of edges"  
    start_time = time.clock()
    
    (h,w) = pic.fatFiducialDataPlane.shape
    numScaleLevels = pic.numScaleLevels
    scaleLevels = pic.scaleLevels
    
    #// B. SET UP NECESSARY DATA STRUCTURES
    fiducialDOGScaleSpace = DOGScaleSpace(pic)     #// this is the datastructure that is really needed
  
    print "Computing eidolon..."
      
    #// C. CONSTRUCT THE EIDOLON
    eidolonDataPlane = np.zeros(pic.fatFiducialDataPlane.shape)
    
    eidolonDataPlane = CoherentDisarray(fiducialDOGScaleSpace, reach, w, h, MAX_SIGMA, numScaleLevels, scaleLevels) # // here the eidolon is constructed
  
    eidolon = pic.DisembedDataPlane(eidolonDataPlane) #// convert dataplane to image
    
    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    return Image.fromarray(eidolon, 'L')


#/* *************************************************************************** */
#/* *************************************************************************** */
#/* SYNTHESIS FROM DISARRAYED SIMPLE CELL ACTIVITY                              */
#/*                                                                             */
#/* This is a fairly complicated example. However, it uses only "straight out   */
#/* the box" toolkit methods.                                                   */
#/* *************************************************************************** */
def OrientedEdgeDisarray(pic, reach, INTEGRATION_FUDGE_FACTOR):
    MAX_SIGMA = pic.MAX_SIGMA       
    THREE_FOURTH = 3.0/4
    
    print "Embarking on synthesis from simple cell activity"  
    start_time = time.clock()
    
    (h,w) = pic.fatFiducialDataPlane.shape
    picSize = pic.fatFiducialDataPlane.shape
    numScaleLevels = pic.numScaleLevels
    scaleLevels = pic.scaleLevels
    
    #// B. SET UP NECESSARY DATA STRUCTURES
    rockBottomPlaneGenerator = RockBottomPlane(pic)       # // needed because of the lowest resolution plane
    print("rockBottomPlaneGenerator --- %s seconds ---" % ( time.clock() - start_time))
    start_time = time.clock()
    
    fiducialSecondOrder = FiducialSecondOrder(pic)        # // sets up the simple cell activity
    print("fiducialSecondOrder --- %s seconds ---" % ( time.clock() - start_time))
    start_time = time.clock()

    #  surface.setTitle("Computing scale spaces P ...");               
    noiseStackPx = CoherentRandomGaussianDataStack(numScaleLevels, w, h, MAX_SIGMA, scaleLevels)
    noiseStackPy = CoherentRandomGaussianDataStack(numScaleLevels, w, h, MAX_SIGMA, scaleLevels)
    #  surface.setTitle("Computing scale spaces Q ...");
    noiseStackQx = CoherentRandomGaussianDataStack(numScaleLevels, w, h, MAX_SIGMA, scaleLevels)
    noiseStackQy = CoherentRandomGaussianDataStack(numScaleLevels, w, h, MAX_SIGMA, scaleLevels)
    #  surface.setTitle("Computing scale spaces R ...");
    noiseStackRx = CoherentRandomGaussianDataStack(numScaleLevels, w, h, MAX_SIGMA, scaleLevels)
    noiseStackRy = CoherentRandomGaussianDataStack(numScaleLevels, w, h, MAX_SIGMA, scaleLevels)
    print("Building noise stacks --- %s seconds ---" % ( time.clock() - start_time))
    start_time = time.clock()
  
    print "Computing eidolon..."
      
    #// C. CONSTRUCT THE EIDOLON
    eidolonDataPlane = np.zeros(pic.fatFiducialDataPlane.shape)

    planeP, planeQ, planeR = StackDisarrayDiscrete(
                                fiducialSecondOrder, 
                                noiseStackPx, 
                                noiseStackPy, 
                                noiseStackQx, 
                                noiseStackQy, 
                                noiseStackRx, 
                                noiseStackRy, 
                                reach, 
                                numScaleLevels, 
                                scaleLevels, 
                                picSize)

    eidolonDataPlane = rockBottomPlaneGenerator.next() + (planeP + planeQ + planeR) * THREE_FOURTH * INTEGRATION_FUDGE_FACTOR
               
    eidolon = pic.DisembedDataPlane(eidolonDataPlane) #// convert dataplane to image
    
    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    return Image.fromarray(eidolon, 'L')


#/* *************************************************************************** */
#/* *************************************************************************** */
#/* DISPLAY A SCALESPACE LEVEL                                                  */
#/* This example shows a pattern that is slightly different from the usual.     */
#/* You can use this technique to view any dataplane you want.                  */
#/* *************************************************************************** */
def DisplayScalespaceLevel(pic, theLevel):   
    if theLevel >= pic.numScaleLevels:
        raise ValueError('Not so much levels as ' + str(theLevel) + ', maximum = ' + str(pic.numScaleLevels - 1) + ' (0 - ' + str(pic.numScaleLevels - 1) + ')!')
    
    print "Embarking on display scalespace level"
    start_time = time.clock()
    
    #// B. SET UP NECESSARY DATA STRUCTURES
    fiducialScaleSpace = ScaleSpace(pic)     
  
    print "Computing eidolon..."
      
    #// C. CONSTRUCT THE EIDOLON
    eidolonDataPlane = np.zeros(pic.fatFiducialDataPlane.shape)
  
    i = 0
    while i <= theLevel:
        eidolonDataPlane = fiducialScaleSpace.next()
        i += 1
 
    eidolon = pic.DisembedDataPlane(eidolonDataPlane) #// convert dataplane to image
    
    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    return Image.fromarray(eidolon, 'L')


#/* *************************************************************************** */
#/* *************************************************************************** */
#/* DISPLAY A DISARRAY FIELD                                                    */
#/* This example shows a fractal disarray field.                                */
#/* this is useful! But you can do MUCH better by constructing both Cartesian   */
#/* components and displaying them by way of a nice colored code (hue for       */
#/* direction, saturation for magnitude, or whatever you fancy).                */
#/* *************************************************************************** */
def DisplayFractalNoise(pic, theLevel):         
    MAX_SIGMA = pic.MAX_SIGMA       
    print "Embarking on display fractal noise"  
    start_time = time.clock()
    
    (h,w) = pic.fatFiducialDataPlane.shape
    numScaleLevels = pic.numScaleLevels
    scaleLevels = pic.scaleLevels
    
    #// B. SET UP NECESSARY DATA STRUCTURES
    fiducialScaleSpace = ScaleSpace(pic)     
  
    print "Computing eidolon..."
      
    #// C. CONSTRUCT THE EIDOLON
    eidolonDataPlane = np.zeros(pic.fatFiducialDataPlane.shape)

    coherentRandomGaussianDataStack = CoherentRandomGaussianDataStack(numScaleLevels, w, h, MAX_SIGMA, scaleLevels)

    i = 0
    while i <= theLevel:
        eidolonDataPlane = coherentRandomGaussianDataStack.next()
        i += 1
        
    eidolon = DataToImage(pic.DisembedDataPlane(eidolonDataPlane, clip=False)) #// convert dataplane to image
    
    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    return Image.fromarray(eidolon, 'L')


#/* *************************************************************************** */
#/* *************************************************************************** */
#/* COHERENT DISARRAY OF EDGES IN RGB CHANNELS                                  */
#/* This example shows a pattern that is slightly different from the usual.     */
#/* The image channels are independently disarrayed, resulting in annoying      */
#/* "chromatic aberration".                                                     */
#/* Doing the disarray in opponent channels yields far nicer results, you can   */
#/* implement it on the same pattern.                                           */
#/* *************************************************************************** */
def CoherentDisarrayOfEdgesRGBChannels(pic, reach):
    MAX_SIGMA = pic.MAX_SIGMA       
    print "Embarking on coherent disarray of edges RGB channels" 
    start_time = time.clock()
    
    (h,w) = pic.fatFiducialDataPlane.shape
    numScaleLevels = pic.numScaleLevels
    scaleLevels = pic.scaleLevels
    
    #// B. SET UP NECESSARY DATA STRUCTURES
    fiducialScaleSpace = ScaleSpace(pic)     
  
    print "Computing eidolon..."                         
  
    #// B & C COMBINED AND REPEATED THREE TIMES
 
    print "Computing RED scale spaces..."
    pic.color = 'red' # setting the color variable to ensure the right data plane is used
    fiducialDOGScaleSpace = DOGScaleSpace(pic) 
    print "Computing RED eidolon dataplane..."
    eidolonRedDataPlane = CoherentDisarray(fiducialDOGScaleSpace, reach, w, h, MAX_SIGMA, numScaleLevels, scaleLevels)
      
    print "Computing GREEN scale spaces..."
    pic.color = 'green' # setting the color variable to ensure the right data plane is used
    fiducialDOGScaleSpace = DOGScaleSpace(pic) 
    print "Computing GREEN eidolon dataplane..."
    eidolonGreenDataPlane = CoherentDisarray(fiducialDOGScaleSpace, reach, w, h, MAX_SIGMA, numScaleLevels, scaleLevels)
      
    print "Computing BLUE scale spaces..."
    pic.color = 'blue' # setting the color variable to ensure the right data plane is used
    fiducialDOGScaleSpace = DOGScaleSpace(pic) 
    print "Computing BLUE eidolon dataplane..."
    eidolonBlueDataPlane = CoherentDisarray(fiducialDOGScaleSpace, reach, w, h, MAX_SIGMA, numScaleLevels, scaleLevels)
    
    pic.color = None # restore to b&w, just in case
    
    eidolon = np.zeros((pic.embeddingData['h'], pic.embeddingData['w'], 3)).astype('uint8')
    
    eidolon[:,:,0] = DataToImage(pic.DisembedDataPlane(eidolonRedDataPlane, clip=False))
    eidolon[:,:,1] = DataToImage(pic.DisembedDataPlane(eidolonGreenDataPlane, clip=False))
    eidolon[:,:,2] = DataToImage(pic.DisembedDataPlane(eidolonBlueDataPlane, clip=False))
    
#    eidolon[:,:,0] = pic.DisembedDataPlane(eidolonRedDataPlane)
#    eidolon[:,:,1] = pic.DisembedDataPlane(eidolonGreenDataPlane)
#    eidolon[:,:,2] = pic.DisembedDataPlane(eidolonBlueDataPlane)
        
    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    return Image.fromarray(eidolon)


#/* *************************************************************************** */
#/* *************************************************************************** */
#/* DISPLAY THE FIRST ORDER GRADIENT DIRECTION                                  */
#/* This example shows a pattern that is slightly different from the usual.     */
#/* The gradient direction is displayed in graylevel, 0-360 degrees mapped on   */
#/* the black-white scale. Notice that angles are periodic, thus black==white!  */
#/* The apparent singularity is just the branch cut.                            */
#/* It would have been better to use the (periodic!) hue circle here. But this  */
#/* is just an instructive example.                                             */
#/* A good exercise (not too hard) is to dosplay the line finder orientation    */
#/* map. Such displays are often convenient when you need to "debug" a novel    */
#/* method.                                                                     */
#/* *************************************************************************** */
def DisplayGradientDirection(pic, theLevel):  
    if theLevel >= pic.numScaleLevels:
        raise ValueError('Not so much levels as ' + str(theLevel) + ', maximum = ' + str(pic.numScaleLevels - 1) + ' (0 - ' + str(pic.numScaleLevels - 1) + ')!')

    print "Embarking on display gradient direction"
    start_time = time.clock()
    
    #// B. SET UP NECESSARY DATA STRUCTURES
    # from Jan -> constructFiducialFirstOrder();
    fiducialFirstOrderXScaleSpace = DifferentialScaleSpace(pic,1,0)
    fiducialFirstOrderYScaleSpace = DifferentialScaleSpace(pic,0,1)
  
    print "Computing eidolon..."
      
    #// C. CONSTRUCT THE EIDOLON
    eidolonDataPlane = np.zeros(pic.fatFiducialDataPlane.shape)
  
    i = 0
    while i <= theLevel:
        gx = fiducialFirstOrderXScaleSpace.next()
        gy = fiducialFirstOrderYScaleSpace.next()
        i += 1

    eidolonDataPlane = gx / np.sqrt(gx**2 + gy**2)
 
    eidolon = DataToImage(pic.DisembedDataPlane(eidolonDataPlane, clip=False))
#    eidolon = pic.DisembedDataPlane(eidolonDataPlane) #// convert dataplane to image
    
    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    return Image.fromarray(eidolon, 'L')


#/* *************************************************************************** */
#/* *************************************************************************** */
#/* SYNTHESIS FROM LARGEST D.O.G. ACTIVITY ONLY                                 */
#/* Here we discard fraction (say 0.9 for 90%) of the DOG activity. Thus only   */
#/* the strongest edges (at any scale) contribute,                              */
#/* This is a straightforward example, except that we change the DOF scalespace */
#/* thus the essential representation. However, this is not a problem if the    */
#/* is one-shot: it doesn't have to clean up the mess it leaves behind.         */
#/* *************************************************************************** */
def SynthesisFromLargestDOGActivity(pic, fraction):
    print "Embarking on synthesis from largest DOG activity"
    start_time = time.clock()

    #// B. SET UP NECESSARY DATA STRUCTURES
    fiducialDOGScaleSpace = DOGScaleSpace(pic)     #// this is the datastructure that is really needed
  
    print "Computing eidolon..."
      
    #// C. CONSTRUCT THE EIDOLON
    eidolonDataPlane = np.zeros(pic.fatFiducialDataPlane.shape)
    #// here the eidolon is constructed

    ssa = SuppressSmallActivity(fiducialDOGScaleSpace, fraction)
    # stack addition
    for data in ssa:
        eidolonDataPlane += data

 
    eidolon = DataToImage(pic.DisembedDataPlane(eidolonDataPlane, clip=False))
#    eidolon = pic.DisembedDataPlane(eidolonDataPlane) #// convert dataplane to image
    
    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    return Image.fromarray(eidolon, 'L')


#/* *************************************************************************** */
#/* *************************************************************************** */
#/* GRADIENT FIELD                                                              */
#/* The first order (gradient) field for a given scale level.                   */
#/* Hue codes direction, saturation codes magnitude                             */
#/* *************************************************************************** */
def GradientField(pic, theLevel):
    if theLevel >= pic.numScaleLevels:
        raise ValueError('Not so much levels as ' + str(theLevel) + ', maximum = ' + str(pic.numScaleLevels - 1) + ' (0 - ' + str(pic.numScaleLevels - 1) + ')!')

    print "Embarking on gradient field"
    start_time = time.clock()
    
    #// B. SET UP NECESSARY DATA STRUCTURES
    # from Jan -> constructFiducialFirstOrder();
    fiducialFirstOrderXScaleSpace = DifferentialScaleSpace(pic,1,0)
    fiducialFirstOrderYScaleSpace = DifferentialScaleSpace(pic,0,1)
  
    print "Computing eidolon..."
      
    #// C. CONSTRUCT THE EIDOLON
    eidolon = np.zeros((pic.embeddingData['h'], pic.embeddingData['w'], 3)).astype('uint8')

    i = 0
    while i <= theLevel:
        gx = fiducialFirstOrderXScaleSpace.next()
        gy = fiducialFirstOrderYScaleSpace.next()
        i += 1  

    eidolonDataPlane = VectorImage(gx, gy)

    #// convert dataplane to image
    eidolon[:,:,0] = pic.DisembedDataPlane(eidolonDataPlane[:,:,0])
    eidolon[:,:,1] = pic.DisembedDataPlane(eidolonDataPlane[:,:,1])
    eidolon[:,:,2] = pic.DisembedDataPlane(eidolonDataPlane[:,:,2])
    
    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    return Image.fromarray(eidolon, 'RGB')


#/* *************************************************************************** */
#/* *************************************************************************** */
#/* LINE FINDER ACTIVITY                                                        */
#/* The second order (gradient) field for a given scale level.                  */
#/* The P, Q, R orietationsare mapped on RGB                                    */
#/* *************************************************************************** */
def LineFinderField(pic, theLevel):
    if theLevel >= pic.numScaleLevels:
        raise ValueError('Not so much levels as ' + str(theLevel) + ', maximum = ' + str(pic.numScaleLevels - 1) + ' (0 - ' + str(pic.numScaleLevels - 1) + ')!')

    print "Embarking on line finder field"
    start_time = time.clock()
    
    #// B. SET UP NECESSARY DATA STRUCTURES
    fiducialSecondOrder = FiducialSecondOrder(pic)  
    
    print "Computing eidolon..."
      
    #// C. CONSTRUCT THE EIDOLON
    eidolon = np.zeros((pic.embeddingData['h'], pic.embeddingData['w'], 3)).astype('uint8')
  
    i = 0
    while i <= theLevel:
        fiducialSecondOrderPScaleSpace, fiducialSecondOrderQScaleSpace, fiducialSecondOrderRScaleSpace = fiducialSecondOrder.next()
        i += 1    

    eidolonDataPlane = LineFinderBasisImage(fiducialSecondOrderPScaleSpace, fiducialSecondOrderQScaleSpace, fiducialSecondOrderRScaleSpace)

    #// convert dataplane to image
    eidolon[:,:,0] = pic.DisembedDataPlane(eidolonDataPlane[:,:,0])
    eidolon[:,:,1] = pic.DisembedDataPlane(eidolonDataPlane[:,:,1])
    eidolon[:,:,2] = pic.DisembedDataPlane(eidolonDataPlane[:,:,2])
    
    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    return Image.fromarray(eidolon, 'RGB')


#/* *************************************************************************** */
#/* *************************************************************************** */
#/* PARTIALLY COHERENT DISARRAY                                                 */
#/*                                                                             */
#/* This implements partially coherent disarray                                 */
#/* *************************************************************************** */
def PartiallyCoherentDisarrayExample(pic, theReach, theCoherence, theGrain):                # // parameter defined in function definition
    MAX_SIGMA = pic.MAX_SIGMA       
    print "Embarking on partially coherent disarray"
    start_time = time.clock()

    #// B. SET UP NECESSARY DATA STRUCTURES
    fiducialDOGScaleSpace = DOGScaleSpace(pic)     #// this is the datastructure that is really needed
  
    (h,w) = pic.fatFiducialDataPlane.shape  
  
    print "Computing eidolon..."
      
    #// C. CONSTRUCT THE EIDOLON
    eidolonDataPlane = np.zeros(pic.fatFiducialDataPlane.shape)

    eidolonDataPlane = PartiallyCoherentDisarray(fiducialDOGScaleSpace, theReach, theCoherence, theGrain, w, h, MAX_SIGMA, pic.numScaleLevels, pic.scaleLevels)

    eidolon = pic.DisembedDataPlane(eidolonDataPlane) #// convert dataplane to image
    
    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    return Image.fromarray(eidolon, 'L')
            

#/* *************************************************************************** */
#/* *************************************************************************** */
#/* COHERENT DISARRAY OF EDGES IN OPPONENT CHANNELS                             */
#/* This example shows a pattern that is slightly different from the usual.     */
#/* This method still follows the standard scheme in principal, yet there are   */
#/* many (necessary) deviations, so it may take you some effort to grok.        */
#/* This example uses superficial disarray, but - of course - that is easy      */
#/* enough to change.                                                           */
#/* *************************************************************************** */
def CoherentDisarrayOfEdgesOpponentChannels(pic, theReach, theGrain):
    print "Embarking on coherent disarray of edges opponent channels"
    start_time = time.clock()
    
    (h,w) = pic.fatFiducialDataPlane.shape
    numScaleLevels = pic.numScaleLevels
    scaleLevels = pic.scaleLevels     
  
    print "Computing eidolon..."   
                      
    #// B & C COMBINED AND REPEATED THREE TIMES
    pic.color = 'red' # have to set a color to make the pic produce color data planes
    # need to cast to float otherwise calculations go horribly wrong!
    r = pic.colorFatFiducialDataPlane['red'].astype('float') # since type = uint8 => cast to type float
    g = pic.colorFatFiducialDataPlane['green'].astype('float') # since type = uint8 => cast to type float
    b = pic.colorFatFiducialDataPlane['blue'].astype('float') # since type = uint8 => cast to type float
    kw, rg, yb = ImageToOpponentRepresentation(r, g, b)
    pic.color = None # reset color, just to be on the safe side
      
    print "BLACK-WHITE disarray going on..."
    xDisplacements = BlurredRandomGaussianDataPlane(w, h, theGrain).blurredRandomGaussianDataPlane
    yDisplacements = BlurredRandomGaussianDataPlane(w, h, theGrain).blurredRandomGaussianDataPlane
    kwPlane = DataPlaneDisarray(kw, xDisplacements, yDisplacements, theReach)
    print "RED_GREEN disarray going on..."
    xDisplacements = BlurredRandomGaussianDataPlane(w, h, theGrain).blurredRandomGaussianDataPlane
    yDisplacements = BlurredRandomGaussianDataPlane(w, h, theGrain).blurredRandomGaussianDataPlane
    rgPlane = DataPlaneDisarray(rg, xDisplacements, yDisplacements, theReach)
    print "YELLOW_BLUE disarray going on..."
    xDisplacements = BlurredRandomGaussianDataPlane(w, h, theGrain).blurredRandomGaussianDataPlane
    yDisplacements = BlurredRandomGaussianDataPlane(w, h, theGrain).blurredRandomGaussianDataPlane
    ybPlane = DataPlaneDisarray(yb, xDisplacements, yDisplacements, theReach)
     
    print "Opponent representation to RGB-image..."
    r, g, b = OpponentRepresentationToImage(kwPlane, rgPlane, ybPlane)
    
    eidolon = np.zeros((pic.embeddingData['h'], pic.embeddingData['w'], 3)).astype('uint8')
    
    eidolon[:,:,0] = pic.DisembedDataPlane(r)
    eidolon[:,:,1] = pic.DisembedDataPlane(g)
    eidolon[:,:,2] = pic.DisembedDataPlane(b)

    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    return Image.fromarray(eidolon, 'RGB')


#/* *************************************************************************** */
#/* *************************************************************************** */
#/* LOTZE DISARRAY OF EDGES IN OPPONENT CHANNELS                                */
#/* Notice how things soon get complicated, although standard tools suffice.    */
#/* *************************************************************************** */
def LotzeDisarrayOfEdgesOpponentChannels(pic, theReach, theGrain):
    print "Embarking on Lotze disarray of edges opponent channels"
    start_time = time.clock()
    
    (h,w) = pic.fatFiducialDataPlane.shape
    numScaleLevels = pic.numScaleLevels  
  
    print "Computing eidolon..."   
                      
    #// B & C COMBINED AND REPEATED THREE TIMES
    print "Computing RED scale spaces..."
    pic.color = 'red' # setting the color variable to ensure the right data plane is used
    fiducialDOGScaleSpaceR = list(DOGScaleSpace(pic)) # dump generator into list for the next step
    print "Computing GREEN scale spaces..."
    pic.color = 'green' # setting the color variable to ensure the right data plane is used
    fiducialDOGScaleSpaceG = list(DOGScaleSpace(pic)) # dump generator into list for the next step
    print "Computing BLUE scale spaces..."
    pic.color = 'blue' # setting the color variable to ensure the right data plane is used
    fiducialDOGScaleSpaceB = list(DOGScaleSpace(pic)) # dump generator into list for the next step

    pic.color = None # reset color just in case 
    
    fiducialDOGScaleSpaceKW = ITORGeneratorKW(fiducialDOGScaleSpaceR, fiducialDOGScaleSpaceG, fiducialDOGScaleSpaceB)
    fiducialDOGScaleSpaceRG = ITORGeneratorRG(fiducialDOGScaleSpaceR, fiducialDOGScaleSpaceG, fiducialDOGScaleSpaceB)
    fiducialDOGScaleSpaceYB = ITORGeneratorYB(fiducialDOGScaleSpaceR, fiducialDOGScaleSpaceG, fiducialDOGScaleSpaceB)

    # get rid of lists, they take too much memory
    del (fiducialDOGScaleSpaceR)
    del (fiducialDOGScaleSpaceG)
    del (fiducialDOGScaleSpaceB)

    print "BLACK-WHITE disarray going on..."
    kwPlane = LotzeDisarray(fiducialDOGScaleSpaceKW,theReach,theGrain, numScaleLevels, w, h);
    print "RED_GREEN disarray going on..."
    rgPlane = LotzeDisarray(fiducialDOGScaleSpaceRG,theReach,theGrain, numScaleLevels, w, h);
    print "YELLOW_BLUE disarray going on..."
    ybPlane = LotzeDisarray(fiducialDOGScaleSpaceYB,theReach,theGrain, numScaleLevels, w, h);
    
    print "Opponent representation to RGB-image..."

    r, g, b = OpponentRepresentationToImage(kwPlane, rgPlane, ybPlane)
    
    eidolon = np.zeros((pic.embeddingData['h'], pic.embeddingData['w'], 3)).astype('uint8')
    
    eidolon[:,:,0] = pic.DisembedDataPlane(r)
    eidolon[:,:,1] = pic.DisembedDataPlane(g)
    eidolon[:,:,2] = pic.DisembedDataPlane(b)

    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    return Image.fromarray(eidolon, 'RGB')


#/* *************************************************************************** */
#/* *************************************************************************** */
#/* COHERENT DISARRAY OF EDGES IN OPPONENT CHANNELS                             */
#/* Notice how things soon get complicated, although standard tools suffice.    */
#/* *************************************************************************** */
def CoherentDisarrayOfOpponentChannels(pic, theReach):
    MAX_SIGMA = pic.MAX_SIGMA  
    print "Embarking on Lotze disarray of edges opponent channels"
    start_time = time.clock()
    
    (h,w) = pic.fatFiducialDataPlane.shape
    numScaleLevels = pic.numScaleLevels     
    scaleLevels = pic.scaleLevels      
  
    print "Computing eidolon..."   
                      
    #// B & C COMBINED AND REPEATED THREE TIMES
    print "Computing RED scale spaces..."
    pic.color = 'red' # setting the color variable to ensure the right data plane is used
    fiducialDOGScaleSpaceR = list(DOGScaleSpace(pic)) # dump generator into list for the next step
    print "Computing GREEN scale spaces..."
    pic.color = 'green' # setting the color variable to ensure the right data plane is used
    fiducialDOGScaleSpaceG = list(DOGScaleSpace(pic)) # dump generator into list for the next step
    print "Computing BLUE scale spaces..."
    pic.color = 'blue' # setting the color variable to ensure the right data plane is used
    fiducialDOGScaleSpaceB = list(DOGScaleSpace(pic)) # dump generator into list for the next step

    pic.color = None # reset color just in case 
    
    fiducialDOGScaleSpaceKW = ITORGeneratorKW(fiducialDOGScaleSpaceR, fiducialDOGScaleSpaceG, fiducialDOGScaleSpaceB)
    fiducialDOGScaleSpaceRG = ITORGeneratorRG(fiducialDOGScaleSpaceR, fiducialDOGScaleSpaceG, fiducialDOGScaleSpaceB)
    fiducialDOGScaleSpaceYB = ITORGeneratorYB(fiducialDOGScaleSpaceR, fiducialDOGScaleSpaceG, fiducialDOGScaleSpaceB)

    # get rid of lists, they take too much memory
    del (fiducialDOGScaleSpaceR)
    del (fiducialDOGScaleSpaceG)
    del (fiducialDOGScaleSpaceB)

    print "BLACK-WHITE disarray going on..."
    kwPlane = CoherentDisarray(fiducialDOGScaleSpaceKW, theReach, w, h, MAX_SIGMA, numScaleLevels, scaleLevels)
    print "RED_GREEN disarray going on..."
    rgPlane = CoherentDisarray(fiducialDOGScaleSpaceRG, theReach, w, h, MAX_SIGMA, numScaleLevels, scaleLevels)
    print "YELLOW_BLUE disarray going on..."
    ybPlane = CoherentDisarray(fiducialDOGScaleSpaceYB, theReach, w, h, MAX_SIGMA, numScaleLevels, scaleLevels)
   
    print "Opponent representation to RGB-image..."

    r, g, b = OpponentRepresentationToImage(kwPlane, rgPlane, ybPlane)
    
    eidolon = np.zeros((pic.embeddingData['h'], pic.embeddingData['w'], 3)).astype('uint8')
    
    eidolon[:,:,0] = pic.DisembedDataPlane(r)
    eidolon[:,:,1] = pic.DisembedDataPlane(g)
    eidolon[:,:,2] = pic.DisembedDataPlane(b)

    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    return Image.fromarray(eidolon, 'RGB')    
    
    
#/* *************************************************************************** */
#/* *************************************************************************** */
#/* HELMHOLTZ DISARRAY OF EDGES IN OPPONENT CHANNELS                            */
#/* Notice how things soon get complicated, although standard tools suffice.    */
#/* *************************************************************************** */
def HelmholtzDisarrayOfEdgesOpponentChannels(pic, theReach):
    MAX_SIGMA = pic.MAX_SIGMA  
    print "Embarking on Helmholtz disarray of edges opponent channels"
    start_time = time.clock()
    
    (h,w) = pic.fatFiducialDataPlane.shape
    numScaleLevels = pic.numScaleLevels     
    scaleLevels = pic.scaleLevels      
  
    print "Computing eidolon..."   
                      
    #// B & C COMBINED AND REPEATED THREE TIMES
    print "Computing RED scale spaces..."
    pic.color = 'red' # setting the color variable to ensure the right data plane is used
    fiducialDOGScaleSpaceR = list(DOGScaleSpace(pic)) # dump generator into list for the next step
    print "Computing GREEN scale spaces..."
    pic.color = 'green' # setting the color variable to ensure the right data plane is used
    fiducialDOGScaleSpaceG = list(DOGScaleSpace(pic)) # dump generator into list for the next step
    print "Computing BLUE scale spaces..."
    pic.color = 'blue' # setting the color variable to ensure the right data plane is used
    fiducialDOGScaleSpaceB = list(DOGScaleSpace(pic)) # dump generator into list for the next step

    pic.color = None # reset color just in case 
    
    fiducialDOGScaleSpaceKW = ITORGeneratorKW(fiducialDOGScaleSpaceR, fiducialDOGScaleSpaceG, fiducialDOGScaleSpaceB)
    fiducialDOGScaleSpaceRG = ITORGeneratorRG(fiducialDOGScaleSpaceR, fiducialDOGScaleSpaceG, fiducialDOGScaleSpaceB)
    fiducialDOGScaleSpaceYB = ITORGeneratorYB(fiducialDOGScaleSpaceR, fiducialDOGScaleSpaceG, fiducialDOGScaleSpaceB)

    # get rid of lists, they take too much memory
    del (fiducialDOGScaleSpaceR)
    del (fiducialDOGScaleSpaceG)
    del (fiducialDOGScaleSpaceB)

    print "BLACK-WHITE disarray going on..."
    kwPlane = HelmholtzDisarray(fiducialDOGScaleSpaceKW, theReach, numScaleLevels, w, h, MAX_SIGMA, scaleLevels)
    print "RED_GREEN disarray going on..."
    rgPlane = HelmholtzDisarray(fiducialDOGScaleSpaceRG, theReach, numScaleLevels, w, h, MAX_SIGMA, scaleLevels)
    print "YELLOW_BLUE disarray going on..."
    ybPlane = HelmholtzDisarray(fiducialDOGScaleSpaceYB, theReach, numScaleLevels, w, h, MAX_SIGMA, scaleLevels)
   
    print "Opponent representation to RGB-image..."

    r, g, b = OpponentRepresentationToImage(kwPlane, rgPlane, ybPlane)
    
    eidolon = np.zeros((pic.embeddingData['h'], pic.embeddingData['w'], 3)).astype('uint8')
    
    eidolon[:,:,0] = pic.DisembedDataPlane(r)
    eidolon[:,:,1] = pic.DisembedDataPlane(g)
    eidolon[:,:,2] = pic.DisembedDataPlane(b)

    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    return Image.fromarray(eidolon, 'RGB')   
    

#/* *************************************************************************** */
#/* *************************************************************************** */
#/* MISSING FROM SCALESPACEREPRESENTATIO                                        */
#/* This example shows what the scalespace representation discards.             */
#/* This is a special method that hopefully shows up as average gray!           */
#/* *************************************************************************** */
def MissingFromScaleSpaceRepresentation(pic):
    MIN_SIGMA = pic.MIN_SIGMA  
    MAX_SIGMA = pic.MAX_SIGMA  
    print"Embarking on missing from SS representation"
    start_time = time.clock()
    
    (h,w) = pic.fatFiducialDataPlane.shape
    numScaleLevels = pic.numScaleLevels     
    scaleLevels = pic.scaleLevels   

    print "Finding missing high resolution detail"
    sp = ScaleSpace(pic)
    fdMinSigma = sp.FuzzyDerivative(pic.fatFiducialDataPlane, 0, 0, MIN_SIGMA)
    fdMaxSigma = sp.FuzzyDerivative(pic.fatFiducialDataPlane, 0, 0, MAX_SIGMA)
    missingHiResData = pic.fatFiducialDataPlane.astype('float') - fdMinSigma
    
    mean = np.mean(pic.fatFiducialDataPlane)
    averageGray = np.ones((h,w)) * mean
       
    missingLoResData = fdMaxSigma - averageGray

    averageGray = np.ones((h,w)) * 127.5          #  // Use uniform average gray for zero level
    dataR = averageGray + missingHiResData        #  // Positive deviations RED, negative GREEN
    dataG = averageGray - missingHiResData
    dataB = averageGray

    dataR += missingLoResData                     #  // Positive deviations YELLOW, negative BLUE
    dataG += missingLoResData
    dataB -= missingLoResData

    eidolon = np.zeros((pic.embeddingData['h'], pic.embeddingData['w'], 3)).astype('uint8')
    
    eidolon[:,:,0] = pic.DisembedDataPlane(dataR)
    eidolon[:,:,1] = pic.DisembedDataPlane(dataG)
    eidolon[:,:,2] = pic.DisembedDataPlane(dataB)

    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    return Image.fromarray(eidolon, 'RGB')   
  

#/* *************************************************************************** */
#/* *************************************************************************** */
#/* DISCRETIZATION ERROR                                                        */
#/* This example shows the discretization error.                                */
#/* *************************************************************************** */
def DiscretizationError(pic, INTEGRATION_FUDGE_FACTOR):
    MIN_SIGMA = pic.MIN_SIGMA  
    MAX_SIGMA = pic.MAX_SIGMA  
    print "Embarking on discretization error"
    start_time = time.clock()

    (h,w) = pic.fatFiducialDataPlane.shape  
    numScaleLevels = pic.numScaleLevels
    scaleLevels = pic.scaleLevels
    picSize = pic.fatFiducialDataPlane.shape
    
    #// B. SET UP NECESSARY DATA STRUCTURES
    fiducialScaleSpace = ScaleSpace(pic);
    fiducialDOGScaleSpace = DOGScaleSpace(pic);
    fiducialLaplacianScaleSpace = FiducialLaplacian(pic);  
  
    print "Computing eidolon..."
      
    #// C. CONSTRUCT THE EIDOLON
    print "Synthesizing and finding difference"

    sp = ScaleSpace(pic)
    fdMinSigma = sp.FuzzyDerivative(pic.fatFiducialDataPlane, 0, 0, MIN_SIGMA)
    fdMaxSigma = sp.FuzzyDerivative(pic.fatFiducialDataPlane, 0, 0, MAX_SIGMA)
    
    error = (fdMinSigma - fdMaxSigma) - StackIntegrate(fiducialLaplacianScaleSpace, numScaleLevels, scaleLevels, picSize) * INTEGRATION_FUDGE_FACTOR
    errorData = QuartilesAndPercentLevels(error)   

    print "Error median is: " + str(errorData['q0.50']) + ", inter quartile range is " + str(errorData['q0.75'] - errorData['q0.25']) + "."

    averageGray = np.ones((h,w)) * 127.5            # // Use uniform average gray for zero level
    dataR = averageGray + error                     # // Positive deviations YELLOW, negative BLUE
    dataG = averageGray + error
    dataB = averageGray - error

    eidolon = np.zeros((pic.embeddingData['h'], pic.embeddingData['w'], 3)).astype('uint8')
    
    eidolon[:,:,0] = pic.DisembedDataPlane(dataR)
    eidolon[:,:,1] = pic.DisembedDataPlane(dataG)
    eidolon[:,:,2] = pic.DisembedDataPlane(dataB)

    print("--- %s seconds ---" % ( time.clock() - start_time))
    
    return Image.fromarray(eidolon, 'RGB')   
  


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
#    pic = Picture('Hanna_Salience.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)
#    pic = Picture('Mona_Lisa.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)


#    im = SynthesisFromDOGActivity(pic) # sanity check
#    im = SynthesisFromLaplacian(pic, INTEGRATION_FUDGE_FACTOR) # sanity check
#    im = SynthesisFromSimpleCellActivity(pic, INTEGRATION_FUDGE_FACTOR) # sanity check  
#    im = SuperficialDisarray(pic, REACH, GRAIN)     
#    im = LotzeTypeDisarray(pic, REACH, GRAIN)
#    im = HelmholtzTypeDisarray(pic, REACH)
#    im = CoherentDisarrayOfEdges(pic, REACH)
#    im = OrientedEdgeDisarray(pic, REACH, INTEGRATION_FUDGE_FACTOR)

#    for theLevel in range(pic.numScaleLevels):
#        print "level =", theLevel
#        im = DisplayScalespaceLevel(pic, theLevel)
#        im = DisplayFractalNoise(pic, theLevel)
#        im = DisplayGradientDirection(pic, theLevel)
#        im = GradientField(pic, theLevel) 
#        im = LineFinderField(pic, theLevel)
#        im.show()

#    im = CoherentDisarrayOfEdgesRGBChannels(pic, REACH) # is dit wel de bedoeling, werkt niet zo bij Jan...        
#    im = SynthesisFromLargestDOGActivity(pic, FRACTION)
    im = PartiallyCoherentDisarrayExample(pic, REACH, INCOHERENCE, GRAIN)
#    im = CoherentDisarrayOfEdgesOpponentChannels(pic, REACH, GRAIN)   
#    im = LotzeDisarrayOfEdgesOpponentChannels(pic, REACH, GRAIN) 
#    im = CoherentDisarrayOfOpponentChannels(pic, REACH)
#    im = HelmholtzDisarrayOfEdgesOpponentChannels(pic, REACH)
#    im = MissingFromScaleSpaceRepresentation(pic)
#    im = DiscretizationError(pic, INTEGRATION_FUDGE_FACTOR)
    
    im.show()
    
   
#==============================================================================
# main
#==============================================================================
if __name__ == "__main__":
    testFunction()

    print "Examples Done!"