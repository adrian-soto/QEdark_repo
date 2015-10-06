# Adrian Soto 
# April 13, 2015
# Stony Brook University
#
#
################################################
#
# Create uniform k-point mesh coordinates.
#
################################################


# number of k-points in each BZ direction
n1=4
n2=4
n3=4

k1min=-0.5*2
k1max= 0.5*2
k2min=-0.5*2
k2max= 0.5*2
k3min=-0.5*2
k3max= 0.5*2

dk1=k1max-k1min
dk2=k2max-k2min
dk3=k3max-k3min

# k-point shifts
s1=0.5
s2=0.5
s3=0.5

# k-point weight
w=2.0/float(n1*n2*n3)


print "          ", n1*n2*n3

for i1 in range(0,n1):
    for i2 in range(0,n2):
        for i3 in range(0,n3):
            
            k1=k1min + (i1 + s1)*dk1/n1 
            k2=k2min + (i2 + s2)*dk2/n2
            k3=k3min + (i3 + s3)*dk3/n3
            
            print " %.6f \t %.6f \t %.6f    \t    %.8f" % (k1, k2, k3, w)
            
