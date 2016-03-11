#==============================================================================
# Imports
#==============================================================================
import numpy as np

from scalespaces import Convolution

#==============================================================================

#==============================================================================
#/* *************************************************************************** */
#/*                                                                             */
#/*                          GENERATE FRESH DATAPLANES                          */
#/*                                                                             */
#/* *************************************************************************** */
#==============================================================================
# class RandomDataPlane
#==============================================================================
class RandomDataPlane(object):            
    """
    Base class for random data planes.
    """
#    def __init__(self, w, h, loc=0.0, scale=127, low=0.0, high=255):     
    def __init__(self, w, h, loc=0.0, scale=1.0, low=-127, high=128):     
#    def __init__(self, w, h, loc=0.0, scale=1.0, low=0.0, high=1):     
#    def __init__(self, w, h, loc=0.0, scale=255, low=0.0, high=1):     
        self.w = w
        self.h = h
        self.loc = loc # for gaussian = mean
        self.scale = scale # for gaussian = standard deviation
        self.low = low # for uniform = lower bound - inclusive
        self.high = high # for uniform = upper bound - exclusive

    @property # returns a random gaussion data plane - values between 0 and 1
    def randomGaussianDataPlane(self):
        data = np.random.normal(loc=self.loc, scale=self.scale, size=(self.h, self.w))
        data = 127.5 * data/np.max(np.abs(data))
        data[data < -127] = -127
        data[data > 127] = 127
        return data
    
    @property # returns a random uniform data plane - values between 0 and 1 exclusive
    def randomUniformDataPlane(self):
        data = np.random.uniform(low=self.low, high=self.high, size=(self.h, self.w))
        data[data > 127] = 127
        return data


#==============================================================================
# class BlurredRandomGaussianDataPlane
#==============================================================================
class BlurredRandomGaussianDataPlane(RandomDataPlane): 
    """
    This builds a Blurred Random Gaussian Data Plane.
    """
    def __init__(self, w, h, sigma):
        super(BlurredRandomGaussianDataPlane, self).__init__(w, h)
        self.sigma = sigma

    @property
    def blurredRandomGaussianDataPlane(self):
        convolution = Convolution(self.randomGaussianDataPlane, 0, 0, self.w, self.h, self.sigma)
        c = (convolution.convolved)
        m = np.mean(c)
        v = np.sqrt(np.var(c))  
        return c - m/v

    
#==============================================================================
# class ScaledBlurredRandomGaussianDataPlane
#==============================================================================
class ScaledBlurredRandomGaussianDataPlane(BlurredRandomGaussianDataPlane): 
    """
    This builds a Scaled Blurred Random Gaussian Data Plane.
    """
    def __init__(self, w, h, sigma, MAX_SIGMA):
        super(ScaledBlurredRandomGaussianDataPlane, self).__init__(w, h, sigma)
        self.MAX_SIGMA = MAX_SIGMA

    @property
    def scaledBlurredRandomGaussianDataPlane(self):
        return self.blurredRandomGaussianDataPlane * (self.sigma / self.MAX_SIGMA)


#==============================================================================
# class DataStack
#==============================================================================
class DataStack(object): 
    """
    This is a data stack generator, mother class.
    """
    def __init__(self, numScaleLevels, w, h):     
        self.numScaleLevels = numScaleLevels
        self.w = w
        self.h = h
        self.current = 0
        
    def __iter__(self):
        return self
        
        
#==============================================================================
# class IncoherentGaussianDataStack
#==============================================================================
class IncoherentGaussianDataStack(DataStack): 
    """
    This is a data stack generator, generates an Incoherent Gaussian Data Stack.
    """
    def __init__(self, numScaleLevels, w, h, sigma):     
        super(IncoherentGaussianDataStack, self).__init__(numScaleLevels, w, h)
        self.sigma = sigma
        
    def next(self):
        if self.current == self.numScaleLevels:
            raise StopIteration("Out of bounds! The number of scale levels is " + str(self.numScaleLevels) + "!")
        else:
            self.current += 1
            b = BlurredRandomGaussianDataPlane(self.w, self.h, self.sigma)
            return b.blurredRandomGaussianDataPlane


#==============================================================================
# class IncoherentScaledGaussianDataStack
#==============================================================================
class IncoherentScaledGaussianDataStack(DataStack): 
    """
    This is a data stack generator, generates an Incoherent Gaussian Data Stack.
    """
    def __init__(self, numScaleLevels, w, h, MAX_SIGMA, scaleLevels):     
        super(IncoherentScaledGaussianDataStack, self).__init__(numScaleLevels, w, h)
        self.MAX_SIGMA = MAX_SIGMA
        self.scaleLevels = scaleLevels

    def next(self):
        if self.current == self.numScaleLevels:
            raise StopIteration("Out of bounds! The number of scale levels is " + str(self.numScaleLevels) + "!")
        else:
            k = self.current
            self.current += 1
            s = ScaledBlurredRandomGaussianDataPlane(self.w, self.h, self.scaleLevels[k], self.MAX_SIGMA)
            return s.scaledBlurredRandomGaussianDataPlane

               
#==============================================================================
# class CoherentRandomGaussianDataStack
#==============================================================================
class CoherentRandomGaussianDataStack(DataStack): 
    """
    // Generate a stack of independent scaled noise planes
    This is a data stack generator, generates an Coherent Random Gaussian Data Stack.
    What this does is that it returns the sum of an Incoherent Scaled Gaussian Data 
    Stack (which consists of a stack of Scaled Blurred Random Gaussian Data Planes)
    as first element, the next element is that sum minus the first matrix in the 
    Incoherent Scaled Gaussian Data Stack, and so on.
    So it returns:
        element 0 = listOfDataPlanes[0] + listOfDataPlanes[1] + ... + listOfDataPlanes[n-1] + listOfDataPlanes[n]
        element 1 = listOfDataPlanes[1] + ... + listOfDataPlanes[n-1] + listOfDataPlanes[n] (= previous - listOfDataPlanes[0])
        element 2 = listOfDataPlanes[2] + ... + listOfDataPlanes[n-1] + listOfDataPlanes[n] (= previous - listOfDataPlanes[1])
        ...
        element n = listOfDataPlanes[n]   
    """
    def __init__(self, numScaleLevels, w, h, MAX_SIGMA, scaleLevels):     
        super(CoherentRandomGaussianDataStack, self).__init__(numScaleLevels, w, h)
        self.MAX_SIGMA = MAX_SIGMA
        self.scaleLevels = scaleLevels

        # dump all dataplanes in a list
        self.listOfDataPlanes = list()        
        for k in range(self.numScaleLevels):
            s = ScaledBlurredRandomGaussianDataPlane(self.w, self.h, self.scaleLevels[k], self.MAX_SIGMA)
            self.listOfDataPlanes.append(s.scaledBlurredRandomGaussianDataPlane)
            if k == 0:
                self.sum = np.zeros(self.listOfDataPlanes[k].shape) # initialize sum matrix
            self.sum += self.listOfDataPlanes[k]
        # reverse list because we'll pop the elements off 
        self.listOfDataPlanes.reverse()
       
    def next(self):
        if self.current == self.numScaleLevels:
            raise StopIteration("Out of bounds! The number of scale levels is " + str(self.numScaleLevels) + "!")
        else:
            k = self.current
            self.current += 1
            if k == 0:
                return self.sum
            else:
                self.sum -= self.listOfDataPlanes.pop() 
                # delete the list of matrices if at end
                if self.current == self.numScaleLevels:
                    del self.listOfDataPlanes 
                return self.sum        


#==============================================================================
# class PartiallyCoherentScaledGaussianDataStack
#==============================================================================
class PartiallyCoherentScaledGaussianDataStack(DataStack): 
    """
    This is a data stack generator, generates an Coherent Random Gaussian Data Stack.
    """
    def __init__(self, numScaleLevels, w, h, sigma, MAX_SIGMA, scaleLevels, degree):     
        super(PartiallyCoherentScaledGaussianDataStack, self).__init__(numScaleLevels, w, h)
        self.degree = degree
        self.incoh = IncoherentGaussianDataStack(numScaleLevels, w, h, sigma)
        self.coh = CoherentRandomGaussianDataStack(numScaleLevels, w, h, MAX_SIGMA, scaleLevels)
       
    def next(self):
        if self.current == self.numScaleLevels:
            raise StopIteration("Out of bounds! The number of scale levels is " + str(self.numScaleLevels) + "!")
        else:
            self.current += 1
            return self.incoh.next() * self.degree + self.coh.next()          





    
#==============================================================================
# 
# Program
# 
#==============================================================================
import time

from math import sqrt
def testFunction(): 
    numScaleLevels = 5
    w = 4
    h = 3
    sigma = sqrt(2)
    MAX_SIGMA = 1
    scaleLevels = (1/sqrt(2), 1, sqrt(2), 2, sqrt(5))

    w = 1024
    h = 512    
    scaleLevels = (0.70710678, 1., 1.41421356, 2., 2.82842712, 4., 5.65685425, 8., 11.3137085, 16., 22.627417, 32., 45.254834, 64.)
    numScaleLevels = len(scaleLevels)
    
    degree = 0.666

    r = RandomDataPlane(w, h)      
    print "randomGaussianDataPlane \n", r.randomGaussianDataPlane    
    print "randomUniformDataPlane \n", r.randomUniformDataPlane
    
    r = RandomDataPlane(w, h, loc=10.0, scale=20)      
    print "randomGaussianDataPlane \n", r.randomGaussianDataPlane    
    print "randomUniformDataPlane \n", r.randomUniformDataPlane
    
    r = RandomDataPlane(w, h, low=10.0, high=100.0)      
    print "randomGaussianDataPlane \n", r.randomGaussianDataPlane    
    print "randomUniformDataPlane \n", r.randomUniformDataPlane
    
    b = BlurredRandomGaussianDataPlane(w, h, sigma)
    print "blurredRandomGaussianDataPlane \n", b.blurredRandomGaussianDataPlane
    
    s = ScaledBlurredRandomGaussianDataPlane(w, h, sigma, MAX_SIGMA)
    print "scaledBlurredRandomGaussianDataPlane \n", s.scaledBlurredRandomGaussianDataPlane

    i1 = IncoherentGaussianDataStack(numScaleLevels, w, h, sigma)
    teller = 1
    for x in i1:
        print teller, "\n"
        teller += 1
        print x
    
    i2 = IncoherentScaledGaussianDataStack(numScaleLevels, w, h, MAX_SIGMA, scaleLevels)
    teller = 1
    for x in i2:
        print teller, "\n"
        teller += 1
        print x
    
    start_time =  time.clock()    
    c = CoherentRandomGaussianDataStack(numScaleLevels, w, h, MAX_SIGMA, scaleLevels)
    teller = 1
    for x in c:
#        print teller, "\n"
#        teller += 1
#        print x
#        x += 1
        pass
    print("===========> --- %s seconds ---" % ( time.clock() - start_time))
    
    start_time =  time.clock()       
    p = PartiallyCoherentScaledGaussianDataStack(numScaleLevels, w, h, sigma, MAX_SIGMA, scaleLevels, degree)
    teller = 1
    for x in p:
#        print teller, "\n"
#        teller += 1
#        print x
        pass
    print("===========> --- %s seconds ---" % ( time.clock() - start_time))




#==============================================================================
# main
#==============================================================================
if __name__ == "__main__":
    testFunction()

    print "Noice Done!"    