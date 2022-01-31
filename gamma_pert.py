#!/usr/bin/env python3

#"""
#Build grib files with smooth noise and time correlation
#"""

import os
import scipy
import numpy
import math
import metview as mv


from scipy import signal, ndimage 

from matplotlib import pyplot


def gamma_scaled(err,size):
    """
    Sample from a gamma distribution with mean=1 and stdev=err
    low errors (below 0.3) approx gaussian shape, large errors (above 1.0) approx lognormal shape
    """

    s=err**2
    a=1/s
    GAMMA=numpy.random.gamma(a,s,size)
    return GAMMA

def BuildRandomSample(err,nlon,nlat,width) :
    """
    Build the random sample at the approximate resolution of the kernel 
    width 2*sigmax (otherwise spread will be smoothed too much to an average value)
    """
    #
    # Build the sample 
    #

    rnlat=int(nlat/width)
    rnlon=int(nlon/width)

    ndim = rnlat * rnlon
    
    m = gamma_scaled(err, ndim).T 

    m = numpy.reshape(m, (rnlat,rnlon))    

    x=numpy.arange(0,rnlat,float(rnlat)/nlat)[0:nlat]
    rx=numpy.arange(rnlat)
    y=numpy.arange(0,rnlon,float(rnlon)/nlon)[0:nlon]
    ry=numpy.arange(rnlon) 

    xx,yy = numpy.meshgrid(rx,ry)  
    f = scipy.interpolate.RectBivariateSpline(rx,ry,m,kx=1,ky=1)
    m = f(x,y)

    #
    # Return 
    #
    return m

def Kernel(n, sigma, nd=1):
  '''
  Generates n-dimensional Gaussian kernel
  n : Number of points in the output window.
  sigma: The standard deviation,
  nd: dimension of the kernel
  '''
  #
  # Creation of the kernel
  #
  #k = scipy.signal.gaussian(n,sigma);
  k = scipy.signal.gaussian(n,sigma)
  for i in range(nd-1): k = k*k[:,None]
  k=k**0.5
  #normalize by the integral to keep the average close to 1
  k=k/k.sum()

  #
  # Return
  #
  return k

def Smooth(x, sigmax):
  '''
  Applies smoothing on n-dimensional data x
  '''
  #
  # Initialisation of the size of the window (3x2sigma)
  #
  sigma = sigmax
  wlength = int(numpy.ceil(sigma * 6))
  #
  # Creation of the kernel
  #
  knl = Kernel(n=wlength,sigma=sigma,nd=x.ndim)
  #
  # Apply kernel
  #
  y = ndimage.convolve(x,knl, mode="wrap")
  #
  # Return
  #
  return y


def BuildMemberPert (imember,err, sigmax, nlon, nlat) : 
    """
    """
    #
    # Build the time correlated smaple
    #
    sample = BuildRandomSample(err,nlon, nlat, 2 * sigmax) 
    #
    data = numpy.reshape(sample, (nlat,nlon))
    out = Smooth(data, sigmax)
    #
    # Retun
    #  
    return out

#
#-------------------------------------------------
# MAIN
#-------------------------------------------------
#
if __name__ == "__main__":
    #
    # User input
    #
    # 
    gribtable=216
    spec="CO"
    ltyp=      ['agr' , 'res' , 'eif', 'shp' , 'swd' , 'tra' , 'bio' , 'fire']
    lgribno=   [ 101  ,  102  ,  103 ,  104  ,  105  ,  106  ,  107  ,  108  ] 
    lerr =     [  1.5 , 0.75  ,  0.5 ,  0.75 , 0.4   , 1.0   ,  0.3  , 1.0   ] #relative error in the random noise
    lsigmax  = [   5  ,     5,      5,    15 ,   5   ,   5   ,  15   ,  2    ] # 2D grid points correlation 50km per cell
    members = 50 
    #
    # Read dummy field
    #
    dummy = mv.read("dummy.grib")
    nlon = int(dummy[0].grib_get_long("Ni"))
    nlat = int(dummy[0].grib_get_long("Nj"))
    #
    # Loop over the members and the sectors
    #

    for typ,gribno,err,sigmax in zip(ltyp,lgribno,lerr,lsigmax):
     print(typ,gribno,err,sigmax)
     outpath="/scratch/cx/cxjb/EDA_SFC_PERT/"+spec+"/"+typ+"/"
     if not os.path.exists(outpath): os.makedirs(outpath)

     for imember in range(members) : 
      print(imember)
      
      scale=2.7 #2.7 is the factor to account for smoothing effect during the convolution, it seems it converges to that number for any err and sigma values
      out = BuildMemberPert (imember, scale*err, sigmax, nlon, nlat)

      new = dummy.set_values(numpy.reshape(out, (nlat*nlon)))
      new = new.grib_set_long(['paramId', int(1000*gribtable+gribno)])

      filename = "pert_m%3.3i_%s_%s.grib" % (imember+1,spec,typ)
      filename = os.path.join(outpath, filename)
      mv.write(filename, new)
