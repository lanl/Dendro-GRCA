'''
@author: Milinda Fernando
@brief: Computes the strain of GWs. 
        Based on getStarin.py shared by David. 

For more details on GW extraction see[1].

[1]. Bishop, N. T., & Reisswig, C. (2014). The gravitational wave strain in the characteristic
formalism of numerical relativity. General Relativity and Gravitation, 46(1), 1643.
'''

import numpy as np
import matplotlib.pyplot as pl
from pylab import *
import scipy as scipy
from scipy import constants
from scipy.fftpack import fft
from scipy.interpolate import interp1d
from scipy.integrate   import cumtrapz
import pandas as pd
import sys as sys
from scipy.interpolate import make_interp_spline, BSpline

'''
space time metric used for propergation time t* computation. 
'''
def metric_g(r,theta,G=scipy.constants.G,M=1,c=1):

   '''
   schwarzchild metric 
   '''

   g=np.zeros((4,4))
   inv_g=np.zeros((4,4))
   g[0,0]=(1-(2 * G * M) / (c** 2 * r)) * c**2
   g[1,1]=-1/(1-( (2*G*M) / (c**2*r)))
   g[2,2]=(-r**2)
   g[3,3]=-r**2 * sin(theta)**2 

   inv_g[0,0]=1.0/g[0,0]
   inv_g[1,1]=1.0/g[1,1]
   inv_g[2,2]=1.0/g[2,2]
   inv_g[3,3]=1.0/g[3,3]

   return [g,inv_g]

'''
compute the propergation time. 
'''

def compute_tstar(t,r,M=1,c=1):
   '''
   compute t^*(t,s) =\int_{0}^{t^*} \frac{\sqrt{-g^{ss}/ g^{tt}}}{1-2M/r} dt^\prime =r -2M ln (r/2M -1)
   '''
   [_,inv_g]=metric_g(r,np.pi/2,c=c)
   rp=r*np.sqrt(c)
   ts=( np.sqrt(-inv_g[1,1]/inv_g[0,0])/(1-2*M/rp) ) * t -rp -2*M*np.log(rp/(2*M)-1)
   return ts   
# code taken from script shared by David. 
##################################################
##     strain --via FFT->-1/omega^2->iFFT       ##
##################################################
def get_strain(t, psi4_real, psi4_img,pts,cutoff):
   
   #m = (int)(len(t)/2)
   #t=t[m:]
   #psi4_real=psi4_real[m:]
   #psi4_img=psi4_img[m:]

   jR    = interp1d(t,psi4_real,kind='cubic')
   jI    = interp1d(t,psi4_img,kind='cubic')

   tU    = np.linspace(t.min(),t.max(),num=pts,endpoint=True)
   psi4_r = jR(tU)
   psi4_i = jI(tU)

   psi4     = psi4_r+1j*psi4_i
   # later we need to change this for propergation time. 
   ct=tU

   fft_psi4   = fft(psi4)
   tstep  = tU[1]-tU[0]
   N      = psi4.size
   #print(N)
   
   # Compute the frequencies of the FFT:
   fft_psi4_f = fftfreq(N,d=tstep)
   # Kill the DC offset:
   fft_psi4[0]= 0.0
   # Change the low-end frequencies
   freqs=fft_psi4_f
   # Make sure zero frequency does not cause error
   freqs[0]=1.0
   f_thres=max(abs(fft_psi4_f))
   #print(f_thres)
   for i in range(len(fft_psi4_f)):
      if (abs(fft_psi4_f[i]) < cutoff):
         freqs[i] = cutoff*sign(fft_psi4_f[i])
      else:
         freqs[i] = fft_psi4_f[i]
   # Do the integration and inverse FFT:
   strain = -ifft(fft_psi4/((2.0*pi*freqs)**2))
   return tU, strain, ct


# code taken from script shared by David. 
##################################################
##     strain2 --via two time integrations      ##
##               interpolates onto uniform grid ##
##################################################
def get_strain2(t, psi4_real, psi4_img,pts):

   # Interpolate to uniform grid:
   tU     = np.linspace(t[0],max(t),num=pts,endpoint=True)
   
   jR    = interp1d(t,psi4_real,kind='cubic')
   jI    = interp1d(t,psi4_img,kind='cubic')

   # coordinate time
   #ct_    = tU#interp1d(t,z)
   
   psi4_r = jR(tU)
   psi4_i = jI(tU)
   psi4     = psi4_r+1j*psi4_i
   ct     = tU#ct_(tU)
   #
   j2p2   = cumtrapz(psi4,tU, initial=0.0)
   #
   # Subtract linear trend of first integral:
   #  (This appears to mess things up and so don't do it)
   #avg    = average(j2p2)
   #tavg   = average(tU)
   #j2p2   = j2p2 - (avg/tavg)*tU
   strain2= cumtrapz(j2p2, tU, initial=0.0)
   #
   #
   # Subtract linear trend:
   #     compute slope only over middle 3/5 of trend
   #     assume intercept is 0
   #start  = len(strain2)/5
   #end    = 3*start
   #avg    = average(strain2[start:end])
   #tavg   = tU[end]-tU[start]
   #strain2= strain2 - (avg/tavg)*tU
   #
   # Use an actual fit:
   coeffs  =polyfit(tU,strain2,1)
   #print(coeffs)
   intercept = coeffs[1]
   slope     = coeffs[0]
   #print("R: intecept=",intercept.real,slope.real,"\n")
   strain2= strain2 - ( slope*tU+intercept )
   #
   return tU, strain2, ct
#####################################################


def plt_strain(fname,lmodes,fprefix):
   df=pd.read_csv(fname,sep='\t')
   df.columns = df.columns.str.replace(' ', '')
   t=np.array(df["t"])

   for l in lmodes:
      fig = plt.figure()
      for m in range(-l,l+1):
         if m >=0: 
            colabs="abs(l_%dm_p_%d)"%(l,abs(m))
            colarg="arg(l_%dm_p_%d)"%(l,abs(m))
         else:
            colabs="abs(l_%dm_m_%d)"%(l,abs(m))
            colarg="arg(l_%dm_m_%d)"%(l,abs(m))

         r=np.array(df[colabs])
         theta=np.array(df[colarg])

         psi4_real=r*cos(theta)
         psi4_img=r*sin(theta)

         [tU,strain,ct]=get_strain(t, psi4_real, psi4_img,10000,4.0)
         #tU=tU[2000:6000]
         #strain=strain[2000:6000]
         #[tU,strain,ct]=get_strain2(t, psi4_real, psi4_img,10000)
         subplot(2*l+1,1,(l+m)+1)
         plt.plot(tU,strain)
         plt.title("strain l=%d,m=%d"%(l,m))
         plt.gcf().set_size_inches(40, 15)
         
      fig.savefig("%s_l%d.png"%(fprefix,l))
      plt.close(fig)

   

def plt_psi4(fname,lmodes,fprefix):
   df=pd.read_csv(fname,sep='\t')
   df.columns = df.columns.str.replace(' ', '')
   t=np.array(df["t"])
   
   for l in lmodes:

      fig = plt.figure()
      for m in range(-l,l+1):
         if m >=0: 
            colabs="abs(l_%dm_p_%d)"%(l,abs(m))
            colarg="arg(l_%dm_p_%d)"%(l,abs(m))
         else:
            colabs="abs(l_%dm_m_%d)"%(l,abs(m))
            colarg="arg(l_%dm_m_%d)"%(l,abs(m))

         r=np.array(df[colabs])
         theta=np.array(df[colarg])

         psi4_real=r*cos(theta)
         psi4_img=r*sin(theta)

         fR = interp1d(t, psi4_real, kind='cubic')
         fI = interp1d(t, psi4_img, kind='cubic')
         tnew=np.linspace(t.min(),t.max(),num=10000,endpoint=True)

         subplot(2*l+1,2,2*(l+m)+1)
         plt.plot(tnew,fR(tnew))
         plt.title("real l=%d,m=%d"%(l,m))
         plt.gcf().set_size_inches(50, 18)
         
         subplot(2*l+1,2,2*(l+m)+2)
         plt.plot(tnew,fI(tnew))
         plt.title("imag l=%d,m=%d"%(l,m))
         plt.gcf().set_size_inches(50, 18)

      fig.savefig("%s_l%d.png"%(fprefix,l))
      plt.close(fig)

   

def main():
   # print command line arguments
   if(len(sys.argv) == 0):
      print("Error: psi4 data file is not provided ")
      sys.exit(0)

   fname=sys.argv[1]
   plt_psi4(fname,[2,3,4],"r1")
   plt_strain(fname,[2,3,4],"r1_strain")

if __name__ == "__main__":
    main()





