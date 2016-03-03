# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 18:07:51 2016

@author: u0076465
"""


from utilities import *        
def BuildScaleSpace(dataPlane, numScaleLevels, scaleLevels):
    scaleStack = list()
    for k in range(numScaleLevels): 
        scaleStack.append(FuzzyDerivative(dataPlane, 0, 0, scaleLevels[k]))
#    rockBottomPlane = np.copy(scaleStack[numScaleLevels-1]) # make copy
    return scaleStack        
        
#==============================================================================
# form convolution 
#==============================================================================
# no clue as to what this is for - from convolution 
#float[][] fuzzyDerivative(float[][] dataPlane, int orderH, int orderV, float sigma) {
def FuzzyDerivative(dataPlane, orderH, orderV, sigma):
    w = dataPlane.shape[1]
    h = dataPlane.shape[0]
    kernel = Kernel(orderV, orderH, w, h, sigma)    
    return Convolve(dataPlane, kernel)   

# convolve - from convolution
def Convolve(dataPlane, kernel):
    dataPlaneFFT = np.fft.rfft2(dataPlane)
    kernelFFT = np.fft.rfft2(kernel)
    return (np.fft.irfft2((dataPlaneFFT * kernelFFT))).astype(np.uint8)



## convolve - from convolution
#def Convolve(dataPlane, kernel):
#    dataPlaneComplex = dataPlane + 0j
#    kernelComplex = kernel + 0j
#
#    dataPlaneFFT = np.fft.fft2(dataPlaneComplex)
#    kernelFFT = np.fft.fft2(kernelComplex)
#
#    return (np.fft.ifft2((dataPlaneFFT * kernelFFT))).real



# no clue as to what this is for - from convolution
def Kernel(xOrder, yOrder, w, h, sigma):  
    v = np.array(range(h/2 + 1) + range(h/2-h + 1, 0))
    column = ((-1.0/(sigma * sqrt(2)))**yOrder) * HermitePolynomial(yOrder, v/(sigma * sqrt(2))) * Gaussian(v, sigma)

    u = np.array(range(w/2 + 1) + range(w/2-w + 1, 0))    
    row = ((-1.0/(sigma * sqrt(2)))**xOrder) * HermitePolynomial(xOrder, u/(sigma * sqrt(2))) * Gaussian(u, sigma)

    return row * column[:, np.newaxis]







#==============================================================================
# 
# General functions
# 
#==============================================================================
## make a fuzzy pic
#def FuzzyDerivative(dataPlane, orderH, orderV, sigma):
#    w = dataPlane.shape[1]
#    h = dataPlane.shape[0]
#    kernel = Kernel(orderV, orderH, w, h, sigma)    
#    return Convolve(dataPlane, kernel)   
#
## convolve - from convolution
#def Convolve(dataPlane, kernel):
#    dataPlaneFFT = np.fft.rfft2(dataPlane)
#    kernelFFT = np.fft.rfft2(kernel)
#    return (np.fft.irfft2((dataPlaneFFT * kernelFFT))).astype(np.uint8)   
#
## make kernel to convolve with
#def Kernel(xOrder, yOrder, w, h, sigma):  
#    v = np.array(range(h/2 + 1) + range(h/2-h + 1, 0))
#    column = ((-1.0/(sigma * sqrt(2)))**yOrder) * HermitePolynomial(yOrder, v/(sigma * sqrt(2))) * Gaussian(v, sigma)
#    u = np.array(range(w/2 + 1) + range(w/2-w + 1, 0))    
#    row = ((-1.0/(sigma * sqrt(2)))**xOrder) * HermitePolynomial(xOrder, u/(sigma * sqrt(2))) * Gaussian(u, sigma)
#    return row * column[:, np.newaxis]
#
##// Hermite polynomials, via the recursion relation - from utilities
#def HermitePolynomial(n, x): 
#    if int(n) != n:
#        raise ValueError('The value of n must be an integer!')        
#    if n == 0:
#        return 1.0
#    elif n == 1:
#        return 2.0 * x
#    else:
#        return 2.0 * (x * HermitePolynomial(n-1, x) - (n-1) * HermitePolynomial(n-2, x))
#
##// Gaussian - from utilities
#def Gaussian(x, sigma):
#    return np.exp(-x**2 / (2.0 * sigma**2)) / (np.sqrt(2 * np.pi) * sigma)        