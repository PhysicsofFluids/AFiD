#### Settings
# Min and max number of gridpoints
minn = 100
maxn = 500
####

print "------------------------------------------"
print "Optimal FFTW resolution" 
print "Between %i and %i grid points" % (minn,maxn)
print "------------------------------------------"

for i in range(0,20):
    for j in range(0,15):
        for k in range(0,10):
            n = 2**i*3**j*5**k 
            if n >= minn and n <= maxn:
                 print "n = " + str(n)
print "------------------END---------------------"

