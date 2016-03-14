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

#==============================================================================
# Eidolon imports
#==============================================================================
from examples import *

#==============================================================================
# 
# CONSTANTS - just to make life easier? 
# 
#==============================================================================

#// Ignore the numbers, they're just labels.
#// See tab "Examples" for the implementations.
#// These are supposed to help you get on to steam.
#// MONOCHROME EIDOLON GENERATION
SUPERFICIAL_DISARRAY = 1  #// random pixel shifting - parameters reach and grain
LOTZE_DISARRAY = 2  #// disarray of boundariness over scales - parameters reach and grain
HELMHOLTZ_DISARRAY = 3  #// disarray of boundariness over scales, scaled by resolution - parameter reach
COHERENT_DISARRAY = 4  #// disarray of boundariness over scales, large RFs drag small RFs along - parameter reach
PARTIAL_COHERENCE = 5  #// mixture of above, parameter reach, incoherence & grain
ORIENTED_EDGE_DISARRAY = 6  #// simple cells of all locations and orientations disarrayed - parameter reach
#// MONOCHROME EIDOLON GENERATION BY SELECTION
LARGEST_DOG_ACTIVITY = 7  #// synthesis using only part of the boundariness - parameter fraction
#// CHROMATIC EIDOLON GENERATION
CHROMATIC_ABERRATION = 8  #// superficial disrray of RGB-channels, yields color fringing - parameters reach and grain
SUPERFICIAL_OPPONENT = 9  #// superficial disrray of opponent-channels - parameters reach and grain
LOTZE_OPPONENT = 10  #// Lotze-type disarray of opponent channels - parameters reach and grain
HELMHOLTZ_OPPONENT = 11  #// Helmholtz-type disarray of opponent channels - parameter reach
COHERENT_OPPONENT = 12  #// coherent-type disarray of opponent channels - parameter reach
#// SANITY CHECKS
SANITY_CHECK_I = 13  #// synthesis from "boundariness" ("D.O.G. activity")
SANITY_CHECK_II = 14  #// synthesis from "boundariness" discrete (Laplacian)
SANITY_CHECK_III = 15  #// synthesis from oriented 2nd order derivatives ("simple cells", "line dtectors")
#// ERROR ESTIMATES
MISSING_FROM_REPRESENTATION = 16  #// Shows what is missing from the scalespace representation
DISCRETIZATION_ERROR = 17  #// discretization error for approximate stack integration
#// DISPLAY OF INTERNAL DATA ATRUCTURES
DISPLAY_SCALESPACE_LEVEL = 18  #// scale space level of fixed resolution - parameter level
GRADIENT_FIELD = 19  #// 1st order directional derivative ("edge detector") field - parameter level
EDGE_DIRECTION_MAP = 20  #// 1st order directional derivative ("edge detector") direction field - parameter level
LINE_FINDER_FIELD = 21  #// 2nd order orientational derivative ("line detector") field - parameter level
#// DISPLAY OF NOISE FIELDS
DISPLAY_FRACTAL_NOISE_PLANE = 22  #// displays a fractal noise plane (one Cartesian coordinate of coherent disarray)



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

    # make a picture class to be used in the examples
#    pic = Picture('Baby_face.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)
#    pic = Picture('berries.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)
#    pic = Picture('DenisHopper.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)
    pic = Picture('Hanna.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)
#    pic = Picture('Hanna_Salience.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)
#    pic = Picture('Mona_Lisa.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)




    #==============================================================================
    # What are we running?
    #==============================================================================
    PRECOOKED_EXAMPLE = DISCRETIZATION_ERROR  #// PICK ANY OF THE CONSTANTS ABOVE & DEFINE PARAMETERS BELOW
    
    # parameters
    REACH = 2.0 #// These are just starting values - adjust as you see fit (try factors of 2 from here for a start,
                #// in some cases you need REALLY large changes before noticing a difference).
    GRAIN = 8.0  #// Not all of these parameters are necessarily relevant for a given method (see above).
    LEVEL = 6  #// These numbers are not of any essence:
    FRACTION = 0.90  #// Once you get the hang of it, you won't need the precooked stuff anyway!    
    INCOHERENCE = 0.10  #// It is imperative to obtain an intuitive feel for the effects of these parameters (and their interactions) though.


    #==============================================================================
    # What are we running?
    #==============================================================================
    if   PRECOOKED_EXAMPLE == SUPERFICIAL_DISARRAY: 
        im = SuperficialDisarray(pic,REACH,GRAIN)
    elif PRECOOKED_EXAMPLE == LOTZE_DISARRAY: 
        im = LotzeTypeDisarray(pic,REACH, GRAIN) 
    elif PRECOOKED_EXAMPLE == HELMHOLTZ_DISARRAY: 
        im = HelmholtzTypeDisarray(pic, REACH) 
    elif PRECOOKED_EXAMPLE == COHERENT_DISARRAY: 
        im = CoherentDisarrayOfEdges(pic, REACH) 
    elif PRECOOKED_EXAMPLE == ORIENTED_EDGE_DISARRAY: 
        im = OrientedEdgeDisarray(pic, REACH, INTEGRATION_FUDGE_FACTOR)
    elif PRECOOKED_EXAMPLE == DISPLAY_SCALESPACE_LEVEL: 
        im = DisplayScalespaceLevel(pic, LEVEL)        
    elif PRECOOKED_EXAMPLE == DISPLAY_FRACTAL_NOISE_PLANE: 
        im = DisplayFractalNoise(pic, LEVEL)        
    elif PRECOOKED_EXAMPLE == CHROMATIC_ABERRATION: 
        im = CoherentDisarrayOfEdgesRGBChannels(pic, REACH)  # different from Jan's version...  
    elif PRECOOKED_EXAMPLE == EDGE_DIRECTION_MAP: 
        im = DisplayGradientDirection(pic, LEVEL) 
    elif PRECOOKED_EXAMPLE == SANITY_CHECK_I: 
        im = SynthesisFromDOGActivity(pic)       
    elif PRECOOKED_EXAMPLE == SANITY_CHECK_III: 
        im = SynthesisFromSimpleCellActivity(pic, INTEGRATION_FUDGE_FACTOR)
    elif PRECOOKED_EXAMPLE == SANITY_CHECK_II: 
        im = SynthesisFromLaplacian(pic, INTEGRATION_FUDGE_FACTOR)
    elif PRECOOKED_EXAMPLE == LARGEST_DOG_ACTIVITY: 
        im = SynthesisFromLargestDOGActivity(pic, FRACTION)
    elif PRECOOKED_EXAMPLE == GRADIENT_FIELD: 
        im = GradientField(pic, LEVEL) 
    elif PRECOOKED_EXAMPLE == LINE_FINDER_FIELD: 
        im = LineFinderField(pic, LEVEL) 
    elif PRECOOKED_EXAMPLE == PARTIAL_COHERENCE: 
        im = PartiallyCoherentDisarrayExample(pic, REACH, INCOHERENCE, GRAIN) 
    elif PRECOOKED_EXAMPLE == SUPERFICIAL_OPPONENT: 
        im = CoherentDisarrayOfEdgesOpponentChannels(pic, REACH,GRAIN) 
    elif PRECOOKED_EXAMPLE == LOTZE_OPPONENT: 
        im = LotzeDisarrayOfEdgesOpponentChannels(pic, REACH, GRAIN) 
    elif PRECOOKED_EXAMPLE == COHERENT_OPPONENT: 
        im = CoherentDisarrayOfOpponentChannels(pic, REACH)
    elif PRECOOKED_EXAMPLE == HELMHOLTZ_OPPONENT: 
        im = HelmholtzDisarrayOfEdgesOpponentChannels(pic, REACH)
    elif PRECOOKED_EXAMPLE == MISSING_FROM_REPRESENTATION: 
        im = MissingFromScaleSpaceRepresentation(pic)
    elif PRECOOKED_EXAMPLE == DISCRETIZATION_ERROR: 
        im = DiscretizationError(pic, INTEGRATION_FUDGE_FACTOR)
    else: 
        im = pic.resizedOriginalImage # just show the original image resized if no valid choices are made
    
    # show the calculated image
    im.show()
    
#==============================================================================
# main
#==============================================================================
if __name__ == "__main__":
    testFunction()

    print "Done!"