#==============================================================================
# Imports
#==============================================================================
import numpy as np

#==============================================================================
# 
# Classes 
# 
#==============================================================================
#==============================================================================
# class Convolution
#==============================================================================
class Convolution(object):            
    """
    This allows you to get a kernel and a convolution of that kernel with a
    given data plane. Both kernel and convolved are accessible as class
    variables, respectively kernel and convolved.
    """
    def __init__(self, dataPlane, xOrder, yOrder, w, h, sigma):     
        self.dataPlane = dataPlane
        self.xOrder = xOrder
        self.yOrder = yOrder
        self.w = w
        self.h = h
        self.sigma = sigma

    @property # this calculates the kernel and is presented as a class variable
    def kernel(self):
        return self.Kernel(self.xOrder, self.yOrder, self.w, self.h, self.sigma)
    
    # this calculates the convolution between the kernel and a given dataplane 
    # and is presented as a class variable 
    @property 
    def convolved(self):
        return self.Convolve(self.dataPlane, self.kernel)
           
#    # convolve - from convolution
#    def Convolve(self, dataPlane, kernel):
#        dataPlaneFFT = np.fft.rfft2(dataPlane)
#        kernelFFT = np.fft.rfft2(kernel) 
#        return np.fft.irfft2((dataPlaneFFT * kernelFFT)) 
              
    # convolve - from convolution - slightly faster with flattened matrices
    def Convolve(self, dataPlane, kernel):
        (h,w) = dataPlane.shape
        return (np.fft.irfft( np.fft.rfft(dataPlane.flatten()) * np.fft.rfft(kernel.flatten()) )).reshape(h,w)
    
    # make kernel to convolve with
    def Kernel(self, xOrder, yOrder, w, h, sigma): 
        v = np.array(range(h/2 + 1) + range(h/2-h + 1, 0))
        column = ((-1.0/(sigma * sqrt(2)))**yOrder) * self.HermitePolynomial(yOrder, v/(sigma * sqrt(2))) * self.Gaussian(v, sigma)  
        u = np.array(range(w/2 + 1) + range(w/2-w + 1, 0))
        row = ((-1.0/(sigma * sqrt(2)))**xOrder) * self.HermitePolynomial(xOrder, u/(sigma * sqrt(2))) * self.Gaussian(u, sigma) 
        return row * column[:, np.newaxis]

    #// Hermite polynomials, via the recursion relation - from utilities
    def HermitePolynomial(self, n, x): 
        if int(n) != n:
            raise ValueError('The value of n must be an integer!')        
        if n == 0:
            return 1.0
        elif n == 1:
            return 2.0 * x 
        else:
            return 2.0 * (x * self.HermitePolynomial(n-1, x) - (n-1) * self.HermitePolynomial(n-2, x))
    
    #// Gaussian - from utilities
    def Gaussian(self, x, sigma):
        return np.exp(-x**2 / (2.0 * sigma**2)) / (np.sqrt(2 * np.pi) * sigma)

            
#==============================================================================
# class DifferentialScaleSpace
#==============================================================================
class DifferentialScaleSpace(object): 
    """
    This allows you to implement any differential order.
    It's the same thing as the Scalespace, but with xOrder and yOrder != 0.
    """
    def __init__(self, picture, xOrder, yOrder):     
        self.pic = picture
        self.current = 0
        self.xOrder = xOrder
        self.yOrder = yOrder
        
    def __iter__(self):
        return self

    def next(self):
        if self.current == self.pic.numScaleLevels:
            raise StopIteration("Out of bounds! The number of scale levels is " + str(self.pic.numScaleLevels) + "!")
        else:
            k = self.current
            self.current += 1
            return self.FuzzyDerivative(self.pic.fatFiducialDataPlane, self.xOrder, self.yOrder, self.pic.scaleLevels[k])

    # make a fuzzy pic
    def FuzzyDerivative(self, dataPlane, orderH, orderV, sigma):
        (h,w) = dataPlane.shape
        c = Convolution(dataPlane, orderH, orderV, w, h, sigma)
        return c.convolved


#==============================================================================
# class ScaleSpace
#==============================================================================
class ScaleSpace(DifferentialScaleSpace): 
    """
    This builds the scalespace.
    It does not have much immediate uses, but you may enjoy watching what the various layers look like.
    It's the same thing as the DifferentialScaleSpace, but with xOrder and yOrder = 0.
    """
    def __init__(self, picture):
        super(ScaleSpace, self).__init__(picture, 0, 0)
        
                
#==============================================================================
# class RockBottomPlane
#==============================================================================
class RockBottomPlane(DifferentialScaleSpace): 
    """
    This builds the rockBottomPlane which is the last value of the scalespace.
    I still use a generator for this, since that conserves memory and also for
    consistency, since the usage is the same as the rest of the scale spaces.
    """
    def __init__(self, picture):
        super(RockBottomPlane, self).__init__(picture, 0, 0)

    def next(self):
        if self.current == 1:
            raise StopIteration("Out of bounds! Only 1 RockBottomPlane!")
        else:
            k = self.current
            self.current += 1
            return self.FuzzyDerivative(self.pic.fatFiducialDataPlane, self.xOrder, self.yOrder, self.pic.scaleLevels[-1])


#==============================================================================
# class DOGScaleSpace
#==============================================================================
class DOGScaleSpace(object): 
    """
    These are actually "scalespace slices", if you add them you regain the image.
    The nature of this decomposition is that it is by "edges" in the phenomenonological
    sense an "edge" is a boundary of two regions and partakes qualities from both
    this is very different from the "edginess" as computable in the first order
    "edginess" simply traces high gradient magnitude.
    """
    def __init__(self, picture):
        self.numScaleLevels = picture.numScaleLevels
        self.scaleStackGenerator = ScaleSpace(picture)
        self.current = 0        
        self.first = self.scaleStackGenerator.next()
       
    def __iter__(self):
        return self

    def next(self):        
        if self.current == self.numScaleLevels:
            raise StopIteration("DOGScaleSpace out of bounds! The number of scale levels is " + str(self.numScaleLevels) + "!")
        else:
            first = np.copy(self.first)
            self.current += 1
            if self.current < self.numScaleLevels:
                self.second = self.scaleStackGenerator.next()
                self.first = self.second    
                second = np.copy(self.second)
                return first - second
            else:   
                return self.second


#==============================================================================
# class FiducialSecondOrder
#==============================================================================
class FiducialSecondOrder(object): 
    """
    This constructs the "line finder" activity, a simple formal model for the representation in primary visual cortex.
    Here I transform from the Hessian matrix representation to a nice isotropic basis of three line finders 
    at 60 degrees orientation differences (fiducialSecondOrderPScaleSpace , fiducialSecondOrderQScaleSpace 
    and fiducialSecondOrderRScaleSpace).
    The coefficients are easy enough to find with a little algebra. 
    """
    def __init__(self, picture):
        self.numScaleLevels = picture.numScaleLevels
        self.hessianXX = DifferentialScaleSpace(picture, 2, 0) #// Compute the Hessian matrix coefficients
        self.hessianXY = DifferentialScaleSpace(picture, 1, 1) #// here the Hessian itself is discarded
        self.hessianYY = DifferentialScaleSpace(picture, 0, 2)
        self.current = 0        
       
    def __iter__(self):
        return self

    def next(self):
        if self.current == self.numScaleLevels:
            raise StopIteration("FiducialSecondOrder out of bounds! The number of scale levels is " + str(self.numScaleLevels) + "!")
        else:
            self.current += 1
            fiducialSecondOrderPScaleSpace = self.hessianXX.next()
            # Here I transform from the Hessian matrix representation to a nice
            # isotropic basis of three line finders at 60 degrees orientation differences.
            # The coefficients are easy enough to find with a little algebra. 
            tmp1 = fiducialSecondOrderPScaleSpace * (1.0/8) + self.hessianYY.next() * (3.0/8)
            tmp2 = self.hessianXY.next() * (-np.sqrt(3)/4.0)
                       
            # fiducialSecondOrderPScaleSpace = hessianXX #// these are the simple cell ("line finder") activities
            fiducialSecondOrderQScaleSpace = tmp1 + tmp2 #// the basis consists of three line finders at 120 degrees orientation increments
            fiducialSecondOrderRScaleSpace = tmp1 - tmp2 #// this suffices to easily find the activity for ANY orientation   

            return fiducialSecondOrderPScaleSpace, fiducialSecondOrderQScaleSpace, fiducialSecondOrderRScaleSpace


#==============================================================================
# class FiducialLaplacian
#==============================================================================                
class FiducialLaplacian(object): 
    """
    This is the true Laplacian defined as a second order differential invariant.
    It is almost identical to the DOG representation, so there is actually little need for it.
    You may want to include it in some "sanity checks" though.
    You will find it easy to construct the Gaussian curvature scalespace in a very similar way.
    """
    def __init__(self, picture):
        self.numScaleLevels = picture.numScaleLevels
        self.hessianXX = DifferentialScaleSpace(picture, 2, 0) #// Compute the Hessian matrix coefficients
        self.hessianYY = DifferentialScaleSpace(picture, 0, 2)
        self.current = 0        
       
    def __iter__(self):
        return self

    def next(self):
        if self.current == self.numScaleLevels:
            raise StopIteration("FiducialLaplacian out of bounds! The number of scale levels is " + str(self.numScaleLevels) + "!")
        else:               
            self.current += 1           
            return self.hessianXX.next() + self.hessianYY.next() 



#==============================================================================
# 
# Program
# 
#==============================================================================
from picture import *

def testFunction(): 
#    SZ = 16
    SZ = 512
    MIN_SIGMA = 1/sqrt(2) #// 1f/sqrt(2f)
    MAX_SIGMA = SZ/4.0
    SIGMA_FACTOR = sqrt(2)  #// sqrt(2f)

    pic = Picture('Hanna.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)
#    pic = Picture('test.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)
    
    dataPlane = pic.fatFiducialDataPlane
    numScaleLevels = pic.numScaleLevels

    scaleStackGenerator = ScaleSpace(pic)
    p = np.zeros(dataPlane.shape)

#==============================================================================
# test scaleStackGenerator
#==============================================================================
    scaleStackGenerator = ScaleSpace(pic)
#    for i in range(numScaleLevels):
#        x = scaleStackGenerator.next()
#        Image.fromarray(x.astype('uint8'), 'L').show()
   
#==============================================================================
# test RockBottomPlane
#==============================================================================
    rockBottomPlaneGenerator = RockBottomPlane(pic)
#    x = rockBottomPlaneGenerator.next()
#    Image.fromarray(x.astype('uint8'), 'L').show()
        
#==============================================================================
# test DOGScaleSpace
#==============================================================================
    DOGScaleStackGenerator = DOGScaleSpace(pic)
#    for i in range(numScaleLevels):
#        ds = DOGScaleStackGenerator.next()
##        Image.fromarray(ds.astype('uint8'), 'L').show()
#        p += ds

#==============================================================================
# test fiducialSecondOrderGenerator
#==============================================================================
    fiducialSecondOrderGenerator = FiducialSecondOrder(pic)
#    for i in range(numScaleLevels):    
#        P, Q, R = fiducialSecondOrderGenerator.next()
#        Image.fromarray(P.astype('uint8'), 'L').show()
#        Image.fromarray(Q.astype('uint8'), 'L').show()
#        Image.fromarray(R.astype('uint8'), 'L').show()


#==============================================================================
# test FiducialLaplacian
#==============================================================================
    fiducialLaplacianGenerator = FiducialLaplacian(pic)
    for i in range(numScaleLevels):    
        x = fiducialLaplacianGenerator.next()
#        print i, "\n", x
        p += x
#        Image.fromarray(x.astype('uint8'), 'L').show()

#==============================================================================
# test DifferentialScaleSpace
#==============================================================================
    fiducialFirstOrderXScaleSpace = DifferentialScaleSpace(pic,1,0)
    fiducialFirstOrderYScaleSpace = DifferentialScaleSpace(pic,0,1)
#    for i in range(numScaleLevels):    
#        x = fiducialFirstOrderXScaleSpace.next()
#        p += x
#        Image.fromarray(x.astype('uint8'), 'L').show()
#        y = fiducialFirstOrderYScaleSpace.next()
#        Image.fromarray(y.astype('uint8'), 'L').show()
#        p += y


# show combination picture
    Image.fromarray(p.astype('uint8'), 'L').show()




#==============================================================================
# main
#==============================================================================
if __name__ == "__main__":
    testFunction()

    print "Scalespaces Done!"