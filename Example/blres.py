import math

#### Input parameters
nu = 33
pr = 1.0
a = 0.922
ra = 1e8
####

# Boundary layer resolution based on Shishkina 2010
pi = math.pi
A = 0.332
E = 0.982
if pr < 3e-4:
  nt = math.sqrt(2*nu)*a*pr**(-3.0/4.0)*A**(3.0/2.0)*pi**(3.0/4.0)
  nv = math.sqrt(2*nu)*a*pr**(-1.0/4.0)*A**(1.0/2.0)*pi**(1.0/4.0)

elif pr >= 3e-4 and pr <= 1:
  nt = math.sqrt(2*nu)*a*pr**(-0.5355+0.033*math.log(pr))
  nv = math.sqrt(2*nu)*a*pr**(-0.1785+0.011*math.log(pr))

elif pr > 1 and pr <= 3:
  nt = math.sqrt(2*nu)*a*pr**(-0.0355+0.033*math.log(pr))
  nv = math.sqrt(2*nu)*a*pr**(0.3215+0.011*math.log(pr))

elif pr > 3:
  nt = math.sqrt(2*nu)*a*E**(3.0/2.0)
  nv = math.sqrt(2*nu)*a*E**(1.0/2.0)*pr**(1.0/3.0)


nt = math.ceil(nt)
print "Number of grid points in thermal BL\t %i" % (nt)

lt = (2*nu)**-1
print "Normalized thermal BL size \t \t %e" % (lt)

nv = math.ceil(nv)
print "Number of grid points in viscous BL\t %i" % (nv)

pi = math.pi
if pr <= 1:
#  h = pi*(pr**2/(ra*nu))**0.25
   h = pr**(1.0/2.0)/(ra*(nu-1))**0.25
else:
#  h = pi*(1/(ra*nu*pr))**0.25
   h = 1.0/(ra*(nu-1))**0.25

print "Maximum grid spacing \t\t \t %.4f" % (h)
n = math.ceil(1/h)
print "Uniform grid resolution for L=1 \t %i" %(n)
