#!/usr/bin/env python3

#"""
#Make some smooth noise...
#"""

import os
import scipy as sp
import numpy as np
import math
import xarray as xr
import yaml as ym
import argparse

from scipy import signal, ndimage 
from matplotlib import pyplot


def gamma_scaled(err,size):
    """
    Sample from a gamma distribution with mean=1 and stdev=err
    low errors (below 0.3) approx gaussian shape, large errors (above 1.0) approx lognormal shape
    """

    s=err**2
    a=1/s
    GAMMA=np.random.gamma(a,s,size)
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
    
    m = gamma_scaled(math.e*err, ndim).T 

    m = np.reshape(m, (rnlat,rnlon))    

    x=np.arange(0,rnlat,float(rnlat)/nlat)[0:nlat]
    rx=np.arange(rnlat)
    y=np.arange(0,rnlon,float(rnlon)/nlon)[0:nlon]
    ry=np.arange(rnlon) 

    xx,yy = np.meshgrid(rx,ry)  
    f = sp.interpolate.RectBivariateSpline(rx,ry,m,kx=1,ky=1)
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
  k = sp.signal.gaussian(n,sigma)
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
  wlength = int(np.ceil(sigma * 6))
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
    data = np.reshape(sample, (nlat,nlon))
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
def main():

    parser = argparse.ArgumentParser(
        description=(
            'NOISE SMOOTHER: Create initial emission perturbations using gamma pdfs')
    )

    required = parser.add_argument_group(title='required arguments')
    required.add_argument(
        '-i', '--yaml_file',
        help="yaml input file",
        type=str, required=True)

    args = parser.parse_args()
    ymlist = ym.load(open(args.yaml_file),Loader=ym.FullLoader)

    filename = ymlist["emission file"]
    varlist = ymlist["sector list"]
    errmagn = ymlist["sector pert"]
    hcorlen = ymlist["sector hcor"]
    geodims = ymlist["geo dims"]
    timedim = ymlist["time dim"]
    members = ymlist["members"]
    issfout = ymlist["scaling factors out"]
    outpath = ymlist["outpath"]
    inpath = ymlist["inpath"]
    domain = ymlist["domain"]
    if domain != "global":
       print("not ready yet for limited area files")
       exit()

    if not os.path.exists(outpath): os.makedirs(outpath)

    #open/read/write/close for each member to avoid adding perturbations 
    for imember in range(members) :

      dsemi = xr.open_dataset(inpath+filename)
       
      nlon = int(np.shape(dsemi[geodims[0]])[0])
      nlat = int(np.shape(dsemi[geodims[1]])[0])
      ntime = int(np.shape(dsemi[timedim])[0])

      #get the cell size at equator
      res = 40075 / nlon
      print(nlon,nlat)

      #
      # Loop over the sectors
      #
      for var,err,hcor in zip(varlist,errmagn,hcorlen):
        hcor_gp = hcor / res
        print(var, err, hcor, hcor_gp)
      
        scal_fac = BuildMemberPert(imember, err, hcor_gp, nlon, nlat)
        scal_fac = np.float32(np.repeat(scal_fac[np.newaxis, :, :], ntime, axis=0)) # need to generalize this
        dsemi[var].values = dsemi[var].values * scal_fac
        if issfout:
          dsemi[var+"_scalfac"] = ([timedim, geodims[1], geodims[0]], scal_fac)
      
      outfile = outpath + filename.split('.')[0]+"_"+f'{imember+1:03d}'+"."+filename.split('.')[1]
      dsemi.to_netcdf(path=outfile,mode='w') 

if __name__ == "__main__":
    main()
