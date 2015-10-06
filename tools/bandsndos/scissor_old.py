#
# Adrian Soto
# 22-12-2014
# Stony Brook University
# 
# Scissor correction to bands structure in
# Quantum Espresso format.
#

class band:

    def __init__(self, numkpoints, bandenergies):
        self.nks = numkpoints
        
        if (len(bandenergies) != numkpoints):
            print "ERROR: list of band energies has wrong length. Setting band to 0."
            self.nrg = [0] * numkpoints
        else:
            self.nrg = bandenergies
        
        

    def shiftband(self, energies, shift):
        energies = map(lambda x : x+shift, energies) # watch for scope here.
        return 




class kpoints:
    
    def __init__(self):
        self.klist = []




def ReadBandStructure(bandsfile, nks, xkpt, nbnd, bsflt): 
    #
    # This function reads the band structure as written 
    # to output of the bands.x program. It returns the bs
    # as a flat list with all energies and another list with
    # the k-point coordinates.
    #
    
    global nks 
    global nbnd 
    global xkpt
    global bsflt

    f = open(bandsfile, 'r')
    
    # First line contains nbnd and nks. Read.
    currentline = f.readline()
    nks = int(currentline[22:26])
    nbnd = int(currentline[12:16])

    print nbnd, nks
    

    # Following lines contain the k-point coordinates
    # and the band energies. 
    
    # Calculate number of lines containing band structure:
    # nks k-point lines
    # At each k-point there are (1+nbnd/10) energy values.
    
    nlpkp = 1+nbnd/10 # Number of Lines Per K-Point 
    nlines = nks + nks * nlpkp

    bsaux = []
    xkpt = []
    
    for ik in range (0, nks):

        currentline = f.readline()
        #kpoint = currentline[12:40]
        kpoint = [float(x) for x in currentline.split()]

        xkpt.append(kpoint)

        auxenerg = []
        for ibnd in range(0, nlpkp):

            currentline = f.readline()
        # append current line to auxiliary list
            auxenerg.append( float(x) for x in currentline.split() )

    # flatten list of lists containing energies for a given kpoint 
    # (each sublist corresponds to one line in the bands.dat file)
        energ = [item for sublist in auxenerg for item in sublist] 
        
        
        # append to band structure
        bsaux.append(energ)


    f.close() 


    
    # Flatten whole band structure and rearrange as a list of band objects
    bsflt = [item for sublist in bsaux for item in sublist]

    return



def ArrangeByBands(numkpts, numbands, bsflat, bs):
    
    global bs
    
    # Ensure bs is empty list
    if not bs:
        for ibnd in range (0, numbands):
            currentband=[]

            for ik in range (0, numkpts):
                currentband.append(bsflat[ik*nbnd+ibnd])    
                bs.append( band(nks, currentband) )

        for ibnd in range (0, nbnd):
            print bs[ibnd].nrg
                
    else: 
        print "ERROR: bs list is not empty"


    return


def GetBands(bandsfile, xk, bs): 
    nks=0
    nbnd=0
    bsflat=[]
    bs=[]
    ReadBandStructure(bandsfile, nks, xk, nbnd, bsflat)
    
    print nks, nbnd
    print bsflat
    
    ArrangeByBands(nks, nbnd, bsflat, bs)
    
    return







#def WriteBandStructure():
#    print ("          %10.6f%10.6f%10.6f" % (kpoint[0], kpoint[1], kpoint[2]) )


#######################

filename="bands.dat"


xk=[]
bs=[]

GetBands(filename, xk, bs) 



#filecontent = f.read()
#print filecontent[1]

#testband = band()
#testband.readband(filecontent)
