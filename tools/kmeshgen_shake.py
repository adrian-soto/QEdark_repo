# Adrian Soto 
# April 13, 2015
# Stony Brook University
#
#
################################################
#
# Create uniform k-point mesh coordinates 
# and give it a random shake with maximum
# amplitude dk
#
################################################

import random

# number of k-points in each BZ direction
n1=4
n2=4
n3=4

k1min=-1.0
k1max=1.0
k2min=-1.0
k2max=1.0
k3min=-1.0
k3max=1.0
dk1=k1max-k1min
dk2=k2max-k2min
dk3=k3max-k3min


# Max shake amplitude
Dk=0.00


# k-point shifts
s1=0.5
s2=0.5
s3=0.5


# Number of k-points by hand
nhand=13

# Total number of kpoints in mesh
n=n1*n2*n3+nhand


# k-point weight
w=2.0/float(n)



# Print mesh

print "          ", n


# Points added by hand
k1, k2, k3 = 0.0, 0.0, 0.0
sh1, sh2, sh3 = Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1)
print " %.6f \t %.6f \t %.6f    \t    %.8f" % (k1+sh1, k2+sh2, k3+sh3, w)
k1, k2, k3 = 0.0, 0.5, 0.5
sh1, sh2, sh3 = Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1)
print " %.6f \t %.6f \t %.6f    \t    %.8f" % (k1+sh1, k2+sh2, k3+sh3, w)
k1, k2, k3 = 0.5, 0.0, 0.5
sh1, sh2, sh3 = Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1)
print " %.6f \t %.6f \t %.6f    \t    %.8f" % (k1+sh1, k2+sh2, k3+sh3, w)
k1, k2, k3 = 0.5, 0.5, 0.0
sh1, sh2, sh3 = Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1)
print " %.6f \t %.6f \t %.6f    \t    %.8f" % (k1+sh1, k2+sh2, k3+sh3, w)
k1, k2, k3 = 0.0, -0.5, -0.5
sh1, sh2, sh3 = Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1)
print " %.6f \t %.6f \t %.6f    \t    %.8f" % (k1+sh1, k2+sh2, k3+sh3, w)
k1, k2, k3 = -0.5, 0.0, -0.5
sh1, sh2, sh3 = Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1)
print " %.6f \t %.6f \t %.6f    \t    %.8f" % (k1+sh1, k2+sh2, k3+sh3, w)
k1, k2, k3 = -0.5, -0.5, 0.0
sh1, sh2, sh3 = Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1)
print " %.6f \t %.6f \t %.6f    \t    %.8f" % (k1+sh1, k2+sh2, k3+sh3, w)
k1, k2, k3 = 0.0, 0.5, -0.5
sh1, sh2, sh3 = Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1)
print " %.6f \t %.6f \t %.6f    \t    %.8f" % (k1+sh1, k2+sh2, k3+sh3, w)
k1, k2, k3 = 0.0, -0.5, 0.5
sh1, sh2, sh3 = Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1)
print " %.6f \t %.6f \t %.6f    \t    %.8f" % (k1+sh1, k2+sh2, k3+sh3, w)
k1, k2, k3 = 0.5, 0.0, -0.5
sh1, sh2, sh3 = Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1)
print " %.6f \t %.6f \t %.6f    \t    %.8f" % (k1+sh1, k2+sh2, k3+sh3, w)
k1, k2, k3 = -0.5, 0.0, 0.5
sh1, sh2, sh3 = Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1)
print " %.6f \t %.6f \t %.6f    \t    %.8f" % (k1+sh1, k2+sh2, k3+sh3, w)
k1, k2, k3 = 0.5, -0.5, 0.0
sh1, sh2, sh3 = Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1)
print " %.6f \t %.6f \t %.6f    \t    %.8f" % (k1+sh1, k2+sh2, k3+sh3, w)
k1, k2, k3 = -0.5, 0.5, 0.0
sh1, sh2, sh3 = Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1), Dk* (2.0*random.random() -1)
print " %.6f \t %.6f \t %.6f    \t    %.8f" % (k1+sh1, k2+sh2, k3+sh3, w)



# Generate regular mesh with a shake

for i1 in range(0,n1):
    for i2 in range(0,n2):
        for i3 in range(0,n3):

            # Calculate shakes
            sh1= Dk* (2.0*random.random() -1)
            sh2= Dk* (2.0*random.random() -1)
            sh3= Dk* (2.0*random.random() -1)
            
            # Regular k-point mesh with shakes
            k1=k1min + (i1 + s1)*dk1/n1 
            k2=k2min + (i2 + s2)*dk2/n2
            k3=k3min + (i3 + s3)*dk3/n3
            
            print " %.6f \t %.6f \t %.6f    \t    %.8f" % (k1+sh1, k2+sh2, k3+sh3, w)
            #print s1, s2, s3





