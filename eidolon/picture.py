#==============================================================================
# Imports
#==============================================================================
from PIL import Image
import numpy as np
from math import log, ceil, sqrt

#==============================================================================
# class Picture
#==============================================================================
class Picture(object): 
    """
    Picture calculates a number of things based on the constants 
    SZ, MIN_SIGMA, MAX_SIGMA and SIGMA_FACTOR:
        - new width and height based on SZ.
        - fiducialImage (basically a grayscale version of the original 
          image (with new width and height)).
        - fatFiducialDataPlane (enlarged fiducialImage with original 
          fiducialImage in the middle, borders are mirrored pieces of 
          fiducialImage). Common trick used to prevent bordereffects.
        - numScaleLevels and scaleLevels based on MIN_SIGMA, 
          MAX_SIGMA and SIGMA_FACTOR.
    It also is able to return the original image, the resized original image 
    and the fiducial image.
    I also added a color variable to make a red, green or blue fatFiducialDataPlane.
    Since all computations are done on the fatFiducialDataPlane, the class also
    has a method for getting the endresult out of the bigger data plane called
    "DisembedDataPlane".
    """
    def __init__(self, theImagePath, SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR):
        self.theImagePath = theImagePath
        self.size = SZ      
        self.embeddingData = self.GetEmbeddingData()    
        self.numScaleLevels, self.scaleLevels = self.BuildScaleLevels(MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)
        # setting the color variable ensures the right data plane is used
        self._color = None
        self.fatFiducialDataPlane = self.EmbedImage(self.FiducialDataPlane(self.fiducialImage))
        self.colorFatFiducialDataPlane = None
        self.MIN_SIGMA = MIN_SIGMA
        self.MAX_SIGMA = MAX_SIGMA
        self.SIGMA_FACTOR = SIGMA_FACTOR

    def ResetFatFiducialDataPlane(self):
        # check if we already have a colorFatFiducialDataPlane, if not, make one
        # we only make the color data plane if we need one, conserves memory
        if self.colorFatFiducialDataPlane == None:
            self.colorFatFiducialDataPlane = self.ColorFatFiducialDataPlane()
            
        if self._color == 'red':
            self.fatFiducialDataPlane = self.colorFatFiducialDataPlane['red']
            # print 'using red dp'
        elif self._color == 'green':
            self.fatFiducialDataPlane = self.colorFatFiducialDataPlane['green']
            # print 'using green dp'
        elif self._color == 'blue':
            self.fatFiducialDataPlane = self.colorFatFiducialDataPlane['blue']
            # print 'using blue dp'
        else:
            self.fatFiducialDataPlane = self.EmbedImage(self.FiducialDataPlane(self.fiducialImage))
            # print 'using b&w dp'
        
    @property
    def color(self):
        return self._color
    
    @color.setter # setting the color variable ensures the right data plane is used
    def color(self, value):
        if value in ('red', 'green', 'blue', None):
            self._color = value
            self.ResetFatFiducialDataPlane()
        else:
            raise ValueError("Wrong color code, acceptable values are: 'red', 'green', 'blue' and None.")

    @property # this doesn't store the image, thus conserves memory
    def originalImage(self):
        return Image.open(self.theImagePath)
        
    @property # this doesn't store the image, thus conserves memory
    def resizedOriginalImage(self):
        width = self.embeddingData['w']
        height = self.embeddingData['h']
        return (Image.open(self.theImagePath)).resize((width, height), resample = Image.LANCZOS).convert('RGB')
        
    @property # this doesn't store the image, thus conserves memory
    def fiducialImage(self):
        # fiducialImage is the resized version of the grayscale image, Jan uses the green channel for this
        return ((self.resizedOriginalImage).convert(mode='L'))    
                 
    @property # this doesn't store the image, thus conserves memory
    def colorFiducialDataPlane(self):
        return self.FiducialDataPlane(self.resizedOriginalImage)

    def ColorFatFiducialDataPlane(self):
        # fatColorFiducialDataPlane is a dictionary of the red, green and blue channels)
        R = self.EmbedImage(self.colorFiducialDataPlane[:,:,0])
        G = self.EmbedImage(self.colorFiducialDataPlane[:,:,1])
        B = self.EmbedImage(self.colorFiducialDataPlane[:,:,2])        
        return {'red':R, 'green':G, 'blue':B}

    #// This sets up the levels of resolution that are of relevance to the fiducial image
    def BuildScaleLevels(self, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR):
        numScaleLevels = int(ceil(2 * log(MAX_SIGMA) / log(2)))
        scaleLevels = np.zeros(numScaleLevels)
        scaleLevels[0] = MIN_SIGMA
        for i in range(1, numScaleLevels):
            scaleLevels[i] = SIGMA_FACTOR * scaleLevels[i-1]
        return numScaleLevels, scaleLevels

    def FiducialDataPlane(self, im):
        return np.array(im) 

    def EmbedImage(self, fiducialDataPlane):
        w  = self.embeddingData['w']
        h  = self.embeddingData['h']
        dw = self.embeddingData['dw']
        dh = self.embeddingData['dh']
        ww = self.embeddingData['ww']
        hh = self.embeddingData['hh']

        fatFiducialDataPlane = np.zeros((self.embeddingData['hh'], self.embeddingData['ww']), dtype=np.uint8)
        
        # fill in middle of fatFiducialDataPlane (0+dh -> dh+h)
        # copies flipped (lr) left part (0->dw) of original picture to left and middle of fatFiducialDataPlane
        fatFiducialDataPlane[dh:h+dh, 0:dw] = np.fliplr(fiducialDataPlane[:, 0:dw])
        # copies original picture in the middle of fatFiducialDataPlane
        fatFiducialDataPlane[dh:h+dh, dw:dw+w] = fiducialDataPlane
        # copies flipped (lr) right part (dw->w) of original picture to right and middle of fatFiducialDataPlane
        fatFiducialDataPlane[dh:h+dh, ww-dw:ww] = np.fliplr(fiducialDataPlane[:, w-dw:w])

        # fill up top and bottom of fatFiducialDataPlane
        # copies flipped (ud) top part with height dh of newly created picture 
        # in middle fatFiducialDataPlane to top of fatFiducialDataPlane
        fatFiducialDataPlane[0:dh, :] = np.flipud(fatFiducialDataPlane[dh:2*dh, :])
        # copies flipped (ud) bottom part with height dh of newly created picture 
        # in middle fatFiducialDataPlane to bottom of fatFiducialDataPlane
        fatFiducialDataPlane[hh-dh:hh, :] = np.flipud(fatFiducialDataPlane[hh-2*dh:hh-dh, :])
        
        return fatFiducialDataPlane

    def DisembedDataPlane(self, dataPlane, clip = True):
        w  = self.embeddingData['w']
        h  = self.embeddingData['h']
        dw = self.embeddingData['dw']
        dh = self.embeddingData['dh']
        if clip:
            dp = (dataPlane[dh:h+dh, dw:w+dw])
            dp[dp < 0] = 0 
            dp[dp > 255] = 255
            return dp.astype('uint8')
        else:
            return dataPlane[dh:h+dh, dw:w+dw]
        
    def GetEmbeddingData(self):
        # this calculates the new size of the picture and the size of the fatFiducialDataPlane
        # it resizes the image and makes it smaller (if size > SZ) or bigger (if size < SZ)
        w, h = self.originalImage.size
        aspectRatio = float(h)/w
        # get new size for picture
        if (aspectRatio<=1):
            w = self.size 
            h = int(round(aspectRatio * self.size))
        else:
            h = self.size
            w = int(round(self.size/aspectRatio))
        
        # get size for fatFiducialDataPlane
        ww = self.PowerOfTwoCeiling((w*3)/2.0)        
        hh = self.PowerOfTwoCeiling((h*3)/2.0)        
        # left, right, top and bottom marge of fatFiducialDataPlane
        dw = int((ww-w)/2.0)
        dh = int((hh-h)/2.0)
        
        return {"w":w, "h":h, "ww":ww, "hh":hh, "dw":dw, "dh":dh}

    # returns number that is a power of 2, bigger than or equal to the given number - from utilities
    def PowerOfTwoCeiling(self, n):
        # don't do numbers smaller than 1
        if n < 1: return 1      
        
        logBase2 = log(n)/log(2)
        # if n is a power of 2, return n, else return the next bigger power of 2 number
        if logBase2 == int(logBase2):
            return n
        else:
            return int(2**(ceil(logBase2)))

            


#==============================================================================
# 
# Program
# 
#==============================================================================
def testFunction(): 
    SZ = 512
#    SZ = 16
    MIN_SIGMA = 1/sqrt(2) #// 1f/sqrt(2f)
    MAX_SIGMA = SZ/4.0
    SIGMA_FACTOR = sqrt(2)  #// sqrt(2f)
#    pic = Picture('test.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)
    pic = Picture('Hanna.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)
#    pic = Picture('Rudy_Dekeerschieter.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)
#    Image.fromarray(pic.fatFiducialDataPlane, 'L').show()
#    pic.originalImage.show()
#    pic.resizedOriginalImage.show()
#    pic.fiducialImage.show()
#    Image.fromarray(pic.DisembedDataPlane(pic.fatFiducialDataPlane), 'L').show()
#    print pic.numScaleLevels, pic.scaleLevels 
    
#    (Image.fromarray(pic.colorFatFiducialDataPlane['red'], 'L')).show()
#    (Image.fromarray(pic.colorFatFiducialDataPlane['green'], 'L')).show()
#    (Image.fromarray(pic.colorFatFiducialDataPlane['blue'], 'L')).show()


#==============================================================================
# main
#==============================================================================
if __name__ == "__main__":
    testFunction()

    print("Picture Done!")