import math
import numpy as np

#This script calculates the maximum STR3 that can be used in AFiD to satisfy the specified number of grid points in the thermal boundary layer. It is valid for ISTR3 = 6: A clipped chebyshev grid distribution.
#EP van der Poel 2015.

#### Settings
# Number of grid points in boundary layer
nbl = 11
# Grid size
nxm = 256
# Estimated Nusselt number
nu = 33
####

blsize = 1.0/(2.0*nu)
nx = nxm + 1
etazm = np.zeros([nx + 300,1])
etaz = np.zeros([nx + 300,1])
xc = np.zeros([nx,1])
xm = np.zeros([nx,1])
maxi = 0
tryuntil = 100
for i in xrange(tryuntil + 1):
   nxmo = nx+i*2
   for kc in xrange(nxmo):
      etazm[kc + 1]=math.cos(math.pi*(float(kc + 1)-0.5)/nxmo)
   for kc in xrange(nx):
      etaz[kc + 1]=etazm[kc + 1 + i]
   delet = etaz[1]-etaz[nx]
   for kc in xrange(nx):
      etaz[kc + 1]=etaz[kc + 1]/(0.5*delet)
   xc[0] = 0
   for kc in xrange(nxm):
      xc[kc + 1] = (1.-etaz[kc + 2])*0.5
   xc[nxm] = 1
   for kc in xrange(nxm):
      xm[kc] = (xc[kc] + xc[kc + 1])*0.5
   for j in xrange(nxm):
      if xm[j] > blsize:
          if j >= nbl:
              maxi = i
          break
if maxi == tryuntil:
   print "No maximum STR3 found lower than " + str(tryuntil) + "."
else:
   print "For N_bl = " + str(nbl) + ", NXM = " + str(nxm) + " and Nu = " + str(nu) + ";"
   print "the maximum STR3 that can be used is " + str(maxi) + "."
