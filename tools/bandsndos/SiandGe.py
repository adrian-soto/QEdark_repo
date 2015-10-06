#
# Adrian Soto
# 22-12-2014
# Stony Brook University
# 
################################################
# Plot band structure and DOS from the
# output of the bands.x program in the 
# Quantum Espresso package.
#
# Features:
#  1) Allows for scissor correction (band shift)
#  2) 
#
################################################
import math
import matplotlib.pyplot as plt
from matplotlib import rcParams
#rcParams['font.family'] = 'serif'
#rcParams['font.serif'] = ['Times']
#rcParams['text.usetex'] = True
#rcParams['font.size'] = 14

class band:
    
    def __init__(self, numkpoints, bandenergies):
        self.nks = numkpoints
        
        if (len(bandenergies) != numkpoints):
            print "ERROR: list of band energies has wrong length. Setting band to 0."
            self.nrg = [0] * numkpoints
        else:
            self.nrg = bandenergies
        

    def printband(self):
        print self.nrg


    def shift(self, delta):
        self.nrg = map(lambda x : x+delta, self.nrg) # watch for scope here.
        return 

    


class kpoints:
    
    def __init__(self):
        self.klist = []



def w0gauss(x):
    # As in flib/w0gauss.f90 in the QE package

    pi = 3.141592653589793
    w0 = 1.0/math.sqrt(pi)*math.exp(-(x-1.0/math.sqrt(2.0))**2)*(2.0-math.sqrt(2.0)*x) 
    return w0




def ReadBandStructure(bandsfile): 
    #
    # This function reads the band structure as written 
    # to output of the bands.x program. It returns the bs
    # as a flat list with all energies and another list with
    # the k-point coordinates.
    #
    
    f = open(bandsfile, 'r')
    
    # First line contains nbnd and nks. Read.
    currentline = f.readline()
    nks = int(currentline[22:26])
    nbnd = int(currentline[12:16])


    # Following lines contain the k-point coordinates
    # and the band energies. 
    
    # Calculate number of lines containing band structure:
    # nks k-point lines
    # At each k-point there are (1+nbnd/10) energy values.
    
    nlpkp = 1+nbnd/10 # Number of Lines Per K-Point 
    nlines = nks + nks * nlpkp

    bsaux = []
    xk = []
    
    for ik in range (0, nks):

        currentline = f.readline()
        #kpoint = currentline[12:40]
        kpoint = [float(x) for x in currentline.split()]

        xk.append(kpoint)

        auxenerg = []
        for ibnd in range(0, nlpkp):

            currentline = f.readline()
        # append current line to auxiliary list
            auxenerg.append( float(x) for x in currentline.split() )

        # flatten list of lists containing energies for a given kpoint 
        # (each sublist corresponds to one line in the bands.dat file)
        energ = [item for sublist in auxenerg for item in sublist] 
        

        # Sort ascendingly band energies for current k-point (to
        # prevent artificial level crossings if QE bands.x output
        # does not sort them correctly) and append to band structure
        bsaux.append(sorted(energ)) 
        
        
    f.close() 
    
    # Flatten bs list
    bsflat = [item for sublist in bsaux for item in sublist]

    return nks, nbnd, xk, bsflat




def SortByBands(nks, nbnd, bsflat):
    # Rearrarange bs from k-points to bands

    bs = []
    for ibnd in range (0, nbnd):
        currentband=[]
        for ik in range (0, nks):
            #currentband.append(bsflat[ik*nbnd+ibnd])    
            bs.append(bsflat[ik*nbnd+ibnd])    
        #bs.append( currentband )

    return bs



def FindHLGap(nks, hvb, lcb):
    #
    # Find HOMO and LUMO energies and energy gap
    # 
    # hvb = highest valence band
    # lcb = lowest conduction band
    #
    # Ehvb = highest valence energy or HOMO energy
    # Elcb = lowest conduction energy or LUMO energy
    #

    gap = lcb[0] - hvb[0]
    for ik1 in range (0, nks):
        auxcond = lcb[ik1]
        for ik2 in range (0, nks):
            auxval = hvb[ik2] 
            currentgap = auxcond-auxval
            if (currentgap < 0.0):
                print "ERROR: negative gap"
            elif (currentgap < gap):
                gap = currentgap



    Ehvb = max(hvb)
    Elcb = min(lcb)

    return Ehvb, Elcb, gap



def Scissor(newgap, bands, nks, nbnd, shifttype):
    #
    # shifttype == 0 : shift valence bands by -0.5*delta and
    #                  conduction bands by 0.5*delta
    # shifttype == 1 : as in 0 but placing the highest valence
    #                  energy at 0.0
    # shifttype == 2 : as in 0 but placing the gap center at 0.0
    #

    
    EHOMO, ELUMO, oldgap = FindHLGap(nks, bands[nval-1].nrg , bands[nval].nrg)
    delta=(newgap-oldgap)/2.0

    # Apply scissor to band structure
    for ibnd in range (0, nbnd):
        if (ibnd < nval):
            bands[ibnd].shift(-1.0*delta)
        else:
            bands[ibnd].shift(delta)



    if (shifttype==0):
        print "Scissor correction to band energies has been applied."

        return

    elif (shifttype==1):

        EHOMO, ELUMO, gap = FindHLGap(nks, bands[nval-1].nrg , bands[nval].nrg)
        delta = -1.0*EHOMO
        #print "delta=", delta
        
        
        for ibnd in range (0, nbnd):
            bands[ibnd].shift(delta)


        print "Scissor correction to band energies has been applied."
        print "Highest valence energy has been set to 0.0 eV"

        

        return
    
    elif (shifttype==2):
        EHOMO, ELUMO, gap = FindHLGap(nks, bands[nval-1].nrg , bands[nval].nrg)
        delta = -0.5*(EHOMO+ELUMO)

        for ibnd in range (0, nbnd):
            bands[ibnd].shift(delta)

        print "Scissor correction to band energies has been applied."
        print "Gap center has been set to 0.0 eV"

        return

    else:
        print "ERROR: shifttype has an non-valid value. Default value shifttype==0."
        print "Scissor correction to band energies has been applied."
                
        return




def CreateDOS(nks, nbnd, deltaE, bnd):

    # ATTENTION: bnd must be an object of the class band

    print "Creating DOS"

    Emin = min(bnd[0].nrg)
    Emax = max(bnd[nbnd-1].nrg)
    ndos = int((Emax - Emin)/deltaE) + 1 # int always rounds to lower integer

    dosE = []
    dosG = []
    intg=0.0


# Create DOS
    for idos in range (0, ndos):
        E = Emin + idos * deltaE
        dosg = 0.0
        for ik in range(0, nks):
            for ibnd in range (0, nbnd):
                dosg = dosg + w0gauss ( (E - bnd[ibnd].nrg[ik] ) / deltaE ) # * wk(ik)
            
        dosg = dosg/deltaE
        intg = intg + dosg*deltaE # integrated DOS

        dosE.append(E)
        dosG.append(dosg)

# Normalize DOS
    dosGnorm=dosG
    for idos in range (0, ndos):
        dosGnorm[idos]=dosGnorm[idos]/intg

    return dosE, dosGnorm





def PlotBandStructure(nbnd, nval, bnd, plotfile, Ef, sympoints, nks_btw_sympoints ):
    #
    # ATTENTION: bnd must be an object of the class band
    #
    # nval: number of valence bands
    # Ef: Fermi Energy. If false then it won't print horizontal line
    # sympoints: list containing labels of symmetry points
    # nks_btw_sympoints: number of k-points between symmetry points
    #
    # NOTE: this function assumes that the number of points
    #       between symmetry points is constant
    #


    plt.clf() #clear figure

    print "Plotting band structure to", plotfile

    col = 'k'
    
    for ibnd in range (0, nbnd):
        #if (ibnd < nval):
        #    col='b'
        #else:
        #    col='r'
        plt.plot(bnd[ibnd].nrg,  markersize=2, linestyle='-', color=col) #marker = 'o')



    y_min = min(bnd[0].nrg)
    y_max = min(bnd[nbnd-1].nrg)


    plt.xlabel("Brillouin zone path")
    plt.ylabel("band energies (eV)")


    

    numsympoints = len(sympoints)

    kpath=[]

    xticks = range(0, numsympoints*nks_btw_sympoints + 1, nks_btw_sympoints)

    for i in range(0, numsympoints):
        kpath.append(sympoints[i])
        if (i < numsympoints-1): 
            for j in range (0, nks_btw_sympoints-1):
                kpath.append('')


                
#    plt.axvline(x=xticks, ymin=0, ymax=1, hold=None, **kwargs)

    # Ticks and vertical lines across BS plot
    plt.xticks(xticks, sympoints)
    for i in range(0,numsympoints):
        plt.axvline(x=xticks[i], ymin=y_min, ymax=y_max, hold=None, color='k', linewidth=0.25)
    

    if (not Ef): 
        plt.axhline(Ef, color="black", linestyle="--")

    plt.xlim( 0, len(bnd[0].nrg)-1 )
    plt.savefig(plotfile)

    return




def PlotDOS(dosE, dosG, plotname):

    # ATTENTION: dosG and dosE must be lists of reals
    
    plt.clf() #clear figure
    plt.plot(dosG, dosE)

    plt.xlabel("Density Of States")
    plt.ylabel("band energies (eV)")

    plt.gca().set_xlim(left=0)
    plt.savefig(plotname)

    return




def PlottwoDOS(dosE1, dosG1, dosE2, dosG2, plotname):

    # ATTENTION: dosG and dosE must be lists of reals
    
    plt.clf() #clear figure
    plt.plot(dosG1, dosE1, color='r', label='Si')
    plt.plot(dosG2, dosE2, color='b', label='Ge')

    plt.xlabel("Density Of States")
    plt.ylabel("band energies (eV)")

    plt.legend(loc=1)
    plt.gca().set_xlim(left=0)
    plt.savefig(plotname)

    return



def PlotMultipleDOS(dosE, dosG, plotname):

    # ATTENTION: dosG and dosE must be lists of lists of reals
    
    Ndos=len(dosE[:])
    
    for i in range(0, Ndos):
        plt.plot(dosG[i], dosE[i])

    plt.xlabel("Density Of States")
    plt.ylabel("band energies (eV)")

    plt.savefig(plotname)

    return



#def WriteBandStructure():
#    print ("          %10.6f%10.6f%10.6f" % (kpoint[0], kpoint[1], kpoint[2]) )





############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################     PROGRAM STARTS HERE     ###################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################


# Datafiles containing band structure
filename=["bands.si.301.dat", "bands.ge.301.dat"]


sympoints=['$\Gamma$', '$X$', '$W$', '$L$', '$\Gamma$', '$K$', '$X$']  

nks_btw_sympoints=50 # To be set by user

    


nks0 = 0
nbnd0=0
xk0=[]
bsflt0=[]
bs0=[]
nks1 = 0
nbnd1=0
xk1=[]
bsflt1=[]
bs1=[]


nks0, nbnd0, xk0, bsflt0 = ReadBandStructure(filename[0]) # Si
nks1, nbnd1, xk1, bsflt1 = ReadBandStructure(filename[1]) # Ge



########################
# For scissor correction
nval = 4
ncond0 = nbnd0 - nval    
ncond1 = nbnd1 - nval    
exptgap0 = 1.1 # eV # Si
exptgap1 = 0.74  # Ge
#########################

# For DOS plot
deltaE = 0.1 #eV
    

# For Si -- 0
if(nbnd0 == 0):
    print "%% ERROR READING BANDS. EXIT %%"
else:
    bs0 = SortByBands(nks0, nbnd0, bsflt0)

    print "nks=", nks0
    print "nbnd=", nbnd0

# Create band objects
    bands=[]
    for ibnd in range (0, nbnd0):
        ledge = ibnd*nks0
        redge = ledge+nks0
        currentband = bs0[ledge:redge]
        bands.append( band(nks0, currentband) )
        print "band ", ibnd+1, " created"
    

    
    Scissor(exptgap0, bands, nks0, nbnd0, 1) # 3rd argument is 1. Then set 3rd argument of PlotBandStructure to 0.0
    print "Scissor correction with gap set to", exptgap0


# Generate DOS to be plotted later
    print "ASC-- checkpoint: creating DOS"
    dosE0, dosG0 = CreateDOS(nks0, nbnd0, deltaE, bands)

# Plot BS
    PlotBandStructure(nbnd0, nval, bands, "Si_BS.pdf", 0.0, sympoints, nks_btw_sympoints)





# For Ge -- 0
if(nbnd0 == 0):
    print "%% ERROR READING BANDS. EXIT %%"
else:
    bs1 = SortByBands(nks1, nbnd1, bsflt1)

    print "nks=", nks1
    print "nbnd=", nbnd1

# Create band objects
    bands=[]
    for ibnd in range (0, nbnd0):
        ledge = ibnd*nks0
        redge = ledge+nks0
        currentband = bs1[ledge:redge]
        bands.append( band(nks1, currentband) )
        print "band ", ibnd+1, " created"
    

    
    Scissor(exptgap1, bands, nks1, nbnd1, 1) # 3rd argument is 1. Then set 3rd argument of PlotBandStructure to 0.0
    print "Scissor correction with gap set to", exptgap1


# Generate DOS to be plotted later
    dosE1, dosG1 = CreateDOS(nks1, nbnd1, deltaE, bands)

# Plot BS
    PlotBandStructure(nbnd1, nval, bands, "Ge_BS.pdf", 0.0, sympoints, nks_btw_sympoints)





PlottwoDOS(dosE0, dosG0, dosE1, dosG1, "DOS.pdf")
