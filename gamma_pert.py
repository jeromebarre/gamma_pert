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
#from matplotlib import pyplot


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
    print("BuildRandomSample start")
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
    print("BuildRandomSample end")
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
    print("Kernel start")
    #
    k = sp.signal.windows.gaussian(n,sigma)
    for i in range(nd-1): k = k*k[:,None]
    k=k**0.5
    #normalize by the integral to keep the average close to 1
    k=k/k.sum()

    #
    print("Kernel end")
    #
    return k

def Smooth(x, sigmax):
    '''
    Applies smoothing on n-dimensional data x
    '''
    #
    print("Smooth start")
    #
    sigma = sigmax
    wlength = int(np.ceil(sigma * 6))
    knl = Kernel(n=wlength,sigma=sigma,nd=x.ndim)
    #
    print("Apply kernel")
    #
    #y = ndimage.convolve(x,knl, mode="wrap")
    y = signal.fftconvolve(x, knl, mode='same')
    #
    print("Smooth end")
    #
    return y


def BuildMemberPert (imember,err, sigmax, nlon, nlat) : 
    """
    """
    #
    print("BuildMemberPert start")
    #
    sample = BuildRandomSample(err,nlon, nlat, 2 * sigmax) 
    #
    data = np.reshape(sample, (nlat,nlon))
    out = Smooth(data, sigmax)
    #
    print("BuildMemberPert end")
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

    template = ymlist["template file"]
    filenames = ymlist["emission files"]
    spelist = ymlist["species list"]
    seclist = ymlist["sector list"]
    errmagn = ymlist["sector pert"]
    hcorlen = ymlist["sector hcor"]
    geodims = ymlist["geo dims"]
    timedim = ymlist["time dim"]
    members = ymlist["members"]
    issfout = ymlist["scaling factors out"]
    pertoly = ymlist["only scaling"]
    outpath = ymlist["outpath"]
    inpath = ymlist["inpath"]
    domain = ymlist["domain"]
    if domain != "global":
       print("not ready yet for limited area files")
       exit()

    if not os.path.exists(outpath): os.makedirs(outpath)

    for imember in range(members) :

      dstpl = xr.open_dataset(inpath+template)

      nlon = int(np.shape(dstpl[geodims[0]])[0])
      nlat = int(np.shape(dstpl[geodims[1]])[0])
      ntime = int(np.shape(dstpl[timedim])[0])
      res = 40075 / nlon
      print(nlon,nlat)

      # CHANGED: Instead of writing files each time, keep them in memory until all changes are done.
      ds_dict = {}  # CHANGED: Dictionary to hold datasets in memory

      # CHANGED: Load all datasets once before modifications
      for spe, efile in zip(spelist, filenames):  # CHANGED: Moved this block before sector loop
          file_parts = efile.split('.')
          base_name = '.'.join(file_parts[:-1])
          extension = file_parts[-1]
          new_filename = f"{base_name}_pert{imember+1:03d}.{extension}"
          outfile = outpath + new_filename

          if os.path.isfile(outfile):
              dsemi = xr.open_dataset(outfile)
          else:
              dsemi = xr.open_dataset(inpath+efile)

          if dsemi.dims[timedim] != ntime:
              raise ValueError(f"Expected {ntime} time steps in {outfile}, found {dsemi.dims[timedim]}")

          ds_dict[outfile] = dsemi

      for sec,err,hcor in zip(seclist,errmagn,hcorlen):
        hcor_gp = hcor / res
        print(sec, err, hcor, hcor_gp)

        scal_fac_list = []  
        for itime in range(ntime):
            hour_field = BuildMemberPert(imember, err, hcor_gp, nlon, nlat)
            scal_fac_list.append(hour_field)
        scal_fac = np.float32(np.array(scal_fac_list))

        # CHANGED: Apply changes to ds_dict datasets in memory
        for spe, efile in zip(spelist,filenames):
            file_parts = efile.split('.')
            base_name = '.'.join(file_parts[:-1])
            extension = file_parts[-1]
            new_filename = f"{base_name}_pert{imember+1:03d}.{extension}"
            outfile = outpath + new_filename
            dsemi = ds_dict[outfile]

            var = f'{spe}{sec}'
            print(f' writing {var}')
            if pertoly:
                if var in dsemi:
                    dsemi = dsemi.drop_vars(var)
            else:
                dsemi[var].values = dsemi[var].values * scal_fac
            if issfout:
                # CHANGED: Just set variable here in memory
                dsemi[var+"_pert"] = ([timedim, geodims[1], geodims[0]], scal_fac)

            ds_dict[outfile] = dsemi

      # CHANGED: After all modifications, write out once with compression
      for outfile, dsemi in ds_dict.items():
          encoding = {}
          for var_name in dsemi.data_vars:
              if issfout or var_name not in [f'{spe}{sec}' for spe in spelist for sec in seclist if pertoly]:
                  encoding[var_name] = {'zlib': True, 'complevel': 5, 'shuffle': True, 'chunksizes': (24, 360, 180)}  # CHANGED: Apply compression here

          temp_outfile = outfile + '.tmp'
          dsemi.to_netcdf(path=temp_outfile, encoding=encoding, mode='w')  # CHANGED: Single write with compression
          dsemi.close()
          os.rename(temp_outfile, outfile)

if __name__ == "__main__":
    main()
