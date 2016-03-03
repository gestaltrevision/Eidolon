import numpy as np
from math import sqrt
from scipy import signal

from picture import *

# make kernel to convolve with
def Kernel1(xOrder, yOrder, w, h, sigma):  
    v = np.array(range(h/2 + 1) + range(h/2-h + 1, 0))
    column = ((-1.0/(sigma * sqrt(2)))**yOrder) * HermitePolynomial(yOrder, v/(sigma * sqrt(2))) * Gaussian(v, sigma)
    u = np.array(range(w/2 + 1) + range(w/2-w + 1, 0))    
    row = ((-1.0/(sigma * sqrt(2)))**xOrder) * HermitePolynomial(xOrder, u/(sigma * sqrt(2))) * Gaussian(u, sigma)
    return row * column[:, np.newaxis]  

# make kernel to convolve with
def Kernel2(xOrder, yOrder, w, h, sigma):  
    theKernel = np.zeros((w,h))
    row = np.zeros(w)
    column = np.zeros(h)
    for j in range(h):
        v = j
        if (j > h/2):
            v = v - h
        column[j] = ((-1.0/(sigma*np.sqrt(2)))**yOrder) * HermitePolynomial(yOrder,(v/(sigma*np.sqrt(2)))) * Gaussian(v,sigma)

    for i in range(w):
        u = i
        if (i > w/2):
            u = u - w
        row[i] = ((-1.0/(sigma*np.sqrt(2)))**xOrder) * HermitePolynomial(xOrder,(u/(sigma*np.sqrt(2)))) * Gaussian(u,sigma)             
        
    for i in range(w):
        for j in range(h):
            theKernel[i][j] = row[i]*column[j]
    theKernel = np.transpose(theKernel)

    return theKernel
 
#// Hermite polynomials, via the recursion relation - from utilities
def HermitePolynomial(n, x): 
    if int(n) != n:
        raise ValueError('The value of n must be an integer!')        
    if n == 0:
        return 1.0
    elif n == 1:
        return 2.0 * x 
    else:
        return 2.0 * (x * HermitePolynomial(n-1, x) - (n-1) * HermitePolynomial(n-2, x))

#// Gaussian - from utilities
def Gaussian(x, sigma):
    return np.exp(-x**2 / (2.0 * sigma**2)) / (np.sqrt(2 * np.pi) * sigma)    

# convolve - from convolution
def Convolve(dataPlane, kernel):
    dataPlaneFFT = np.fft.rfft2(dataPlane)
    kernelFFT = np.fft.rfft2(kernel) 
    print "dataPlaneFFT \n", dataPlaneFFT
    print "kernelFFT \n", kernelFFT

    return np.fft.irfft2((dataPlaneFFT * kernelFFT))  

# convolve - from convolution
def Convolve2(dataPlane, kernel):
    dataPlaneComplex = dataPlane + 0j
    kernelComplex = kernel + 0j

    #  fft.complexForward( dataPlaneComplex);
    #  fft.complexForward(kernelComplex);
    dataPlaneFFT = np.fft.fft2(dataPlaneComplex)
    kernelFFT = np.fft.fft2(kernelComplex)

    #  complexMultiply(dataPlaneComplex, kernelComplex);  
    #  fft.complexInverse(dataPlaneComplex, true);
    #  extractRealDataPlane(dataPlaneComplex,dataPlane);
    return (np.fft.ifft2((dataPlaneFFT * kernelFFT))).real

# convolve - from convolution
def Convolve3(dataPlane, kernel):
    print 'in convolve'
    return signal.convolve2d(dataPlane, kernel, mode='same', boundary='wrap')

def intPow(val, exp):
    print "val =", val        
    print "exp+1 =", exp+1
    if (exp>0):
        return(val*intPow(val,exp-1))
    elif (exp<0):
        return(intPow(val,exp+1)/val)
    else:
        print "return 1"
        return 1




SZ = 512
MIN_SIGMA = 1/sqrt(2) #// 1f/sqrt(2f)
MAX_SIGMA = SZ/4.0
SIGMA_FACTOR = sqrt(2)  #// sqrt(2f)

pic = Picture('DenisHopper.jpg', SZ, MIN_SIGMA, MAX_SIGMA, SIGMA_FACTOR)              

fiducialImage = pic.fiducialImage
npFiducialImage = np.array(fiducialImage)

(w,h) = fiducialImage.size

xOrder = 0
yOrder = 0

sigma = sqrt(2)

k = Kernel1(xOrder, yOrder, w, h, sigma)
c = Convolve3(npFiducialImage, k)

Image.fromarray(c.astype('uint8'), 'L').show()


#
#print intPow(2.5,0)


#dataPlane = ((np.arange(1,33)).reshape(4,8).astype('float64') * np.random.normal(loc=128, scale=64, size=(4,8))).astype('uint8')
#print dataPlane
#
#
#xOrder = 2
#yOrder = 0
#(h, w) = dataPlane.shape 
#sigma = sqrt(2)
#
#a = Kernel1(xOrder, yOrder, w, h, sigma)
#b = Kernel2(xOrder, yOrder, w, h, sigma)
#
#print "a \n", a
#print
#print "b \n", b
#print
#print "a-b \n", a-b
#
#c1 = Convolve(dataPlane, a)
#c2 = Convolve2(dataPlane, a)
#
#print "c1 \n", c1
#print "c2 \n", c2
#print "c1 - c2 \n", c1 - c2
#
#
#randomGaussianDataPlane = np.random.normal(loc=0, scale=1, size=(h, w)) * 100
#print "randomGaussianDataPlane \n", randomGaussianDataPlane
#
#c3 = Convolve(randomGaussianDataPlane, a)
#print "c3 \n", c3


#m = np.mean(c)
#v = np.sqrt(np.var(c))#, ddof=-1))
#print m
#print v
#print c - m/v
#
#w    = 3
#h    = 3
#x  = 0
#xx = 0
#for i in range(w):
#    columnSum = 0
#    columnSumSquared = 0
#    for j in range(h): 
#        columnSum        += c[i][j]
#        columnSumSquared += (c[i][j])**2 
#
#    x  += columnSum
#    xx += columnSumSquared
#
#x  =  x/(w*h);
#xx = xx/(w*h);
#
#m = x;
#v = sqrt(xx-x**2);
#
#print "m", m
#print "v", v
#print c - m/v


print "done!"