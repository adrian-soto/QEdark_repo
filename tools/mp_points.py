# Adrian Soto
# 21-03-2014
# Stony Brook University
#
#
# Generate the Monkhorst-Pack coefficients
# for a k-mesh
#

import numpy as np

# Constants
pi=3.14159265359
bohrtoang = 0.529177249

# Lattice info
a=10.200000 * bohrtoang

# Reciprocal lattice vectors
b1=np.array([ -1.000000 ,-1.000000 , 1.000000 ])
b1=(2.0*pi/a)*b1

b2=np.array([  1.000000 , 1.000000 , 1.000000 ])
b2=(2.0*pi/a)*b2

b3=np.array([ -1.000000 , 1.000000 , 1.000000 ])
b3=(2.0*pi/a)*b3


# Monkhorst-Pack parameters
q1=6
q2=6
q3=6

s1=0
s2=0
s3=0



kvectors=np.zeros(3*q1*q2*q3).reshape(q1*q2*q3,3)

ik=0
for r1 in range(q1):
    for r2 in range(q2):
        for r3 in range(q3):
            kvectors[ik,0] = (2.0*(r1+1) - q1 - 1)/(2.0*q1)
            kvectors[ik,1] = (2.0*(r2+1) - q2 - 1)/(2.0*q2)
            kvectors[ik,2] = (2.0*(r3+1) - q3 - 1)/(2.0*q3)
            
#            print ik+1, kvectors[ik,0]*b1 + kvectors[ik,1]*b2 + kvectors[ik,2]*b3
            print ik+1, kvectors[ik,:]

            ik = ik + 1
            

