# Adrian Soto 
# 17-09-2015
# Stony Brook University
#
#
# Generate k-point path for band structure
# for simple cubic lattice system in 
# units of 2*pi/a.
#
#

import math


def eucldist(x1, x2):
    aux = (x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2
    return math.sqrt(aux)




class kpath():

    def __init__(self, numkpointssegment, symmetrypoints):
        
        # total number of k-points in path
        self.nksseg = numkpointssegment 
        self.nsp = len(symmetrypoints)
        
        # Symmetry point labels
        self.label = []
        self.xksym = []
        
        for i in range (0, self.nsp):
            self.label.append(symmetrypoints[i][0])
            self.xksym.append(symmetrypoints[i][1:4])
        
        self.xk=[]
            


    def printsympoints(self):
        # Print symmetry points and their position on the
        # k-point path to screen.
        for i in range(0, self.nsp):
            print self.label[i], self.xksym[i], ' at point ', 1 + i * self.nksseg


        
    def totallength(self):
        # Length of entire k-point path in units of 2*pi/a
        # This can be useful if the band structure horizontal
        # axis needs to be scaled to the lenght of the path
        #
        #
        totlen=0.0
        for i in range(0, self.nsp-1):
           totlen=totlen+eucldist(self.xksym[i], self.xksym[i+1])
      
        return totlen


    def makekpath(self, precision):
        # Create k-point path in BZ for the given symmetry points
        # and the given number of k-points between symmetry points
        #
        # First calculate all points in k-path except the last point
        # (which is a symmetry point) and append. Then append last 
        # symmetry point.
        #
        # precision is set to round off k-point coordinates
        
        numpoints=0
        for isym in range(0, self.nsp-1):
            
            for jinterp in range (0, self.nksseg):

                currentxk=[]
                delta = float(jinterp)/self.nksseg
                for ix in range(0,3):
                    aux = self.xksym[isym][ix] + delta *(self.xksym[isym+1][ix]-self.xksym[isym][ix])
                    currentxk.append(round(aux, precision))
                    
                self.xk.append(currentxk)
                
                numpoints = numpoints+1

        # Append last k-point (symmetry point)
        currentxk = self.xksym[self.nsp-1]
        self.xk.append(currentxk)
        


#############################################
#############################################
#############################################
#############################################




kG = ['G', 0.0, 0.0, 0.0]
kR = ['R', 0.5, 0.5, 0.5]
kX = ['X', 0.0, 0.5, 0.0]
kM = ['M', 0.5, 0.5, 0.0]



symkpts = [kM, kG, kX, kR, kG]

nks_btw_sympoints = 50
kpathfilename="kpath_SC.dat"
precision=5

# Create instace class and initialize
mypath = kpath(nks_btw_sympoints, symkpts)
print " "
print "Creating k-point path by linear interpolation between symmetry points"
mypath.makekpath(precision)
print " "


# Print info to screen
nks = len(mypath.xk) # == len(symkpts)*nks_btw_sympoints + 1
print "The number of k-points in the path is ", nks
print " "

mypath.printsympoints()
print " "

# Formats
floatfmt='%' + str(precision+5) + '.' + str(precision) + 'f '
intfmt='%5d '
fmt = floatfmt + floatfmt + floatfmt + intfmt + '\n'

print "Writing k-point path to ", kpathfilename, '\n'
f1=open(kpathfilename, "w+" )

f1.write('\t' + str(nks) + '\n')
for ik in range(0, nks):
    f1.write(fmt % (mypath.xk[ik][0], mypath.xk[ik][1], mypath.xk[ik][2], ik+1) )
f1.close()


print "EXIT"
