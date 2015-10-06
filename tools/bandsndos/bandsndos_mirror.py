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
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import AutoMinorLocator
import matplotlib.gridspec as gridspec
import csv

plt.rcParams['font.family'] = 'Serif'
plt.rcParams['font.serif'] = 'Times New Roman'

#rcParams['text.usetex'] = True
rcParams['font.size'] = 24

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


################################################
# End of class band
################################################




    
class kpoints:
    
    def __init__(self):
        self.klist = []




class dos:
    
    def __init__(self): #, numE, dosE, dosG, dosI):
        self.numE = 0
        self.dosE = []
        self.dosG = []
        self.dosI = []
        

    def Load(self, dosfile):
        #
        # Load DOS from dos.x output
        #


        print " "
        print "Loading DOS from ", dosfile
        print " "

        
        # Count lines in file
        self.numE=sum(1 for line in open(dosfile))


        # Read file line by line and process
        f=open(dosfile, 'r')
       

        # First line is header. Discard
        data=f.readline()
        
        # Iterate over file lines
        for ilin in range(1,self.numE):
            data=f.readline()
            
            E=float(data[0:7])
            self.dosE.append(E)

            G=float(data[9:19])
            self.dosG.append(G)

            I=float(data[21:31])
            self.dosI.append(I)

        f.close()
        
        return

################################################
# End of class dos
################################################




#
# Global functions
#


def w0gauss(x):
    # As in flib/w0gauss.f90 in the QE package

    pi = 3.141592653589793
    sqrt2=math.sqrt(2)

    arg = min([200.0, (x - 1.0 / sqrt2 ) **2])

    w0 = (1.0/math.sqrt(pi)) * math.exp(-1.0 * arg )*(2.0 - sqrt2*x) 

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



def Scissor(nks, newgap, bands, shifttype):
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




def CreateDOS(nks, nbnd, bzv, Emin, Emax, deltaE, bnd, normalize):

    # ATTENTION: bnd must be an object of the class band



    Emin = min(bnd[10].nrg)
    Emax = max(bnd[nbnd-1].nrg)
    ndos = int((Emax - Emin)/deltaE + 0.50000001) # int always rounds to lower integer 

    dosE = []
    dosG = []
    intg=0.0

    
    deltaEgauss=5.0*deltaE
    d3k=(1.0/nks)*bzv
    wk=2.0/nks


    print "Creating DOS with uniform k-point weights"

# Create DOS
    for idos in range (0, ndos):
        E = Emin + idos * deltaE
        dosg = 0.0
        for ik in range(0, nks):
            for ibnd in range (0, nbnd):
                dosg = dosg + w0gauss ( (E - bnd[ibnd].nrg[ik] ) / deltaEgauss ) * wk
                ###dosg = dosg + w0gauss ( (E - bnd[ibnd].nrg[ik] ) / deltaE ) * wk
            
        dosg = dosg/deltaEgauss
        intg = intg + dosg*deltaE # integrated DOS

        dosE.append(E)
        dosG.append(dosg)



    print "\n Integrated DOS=", intg, "\n"



    
    # Normalize DOS
    if (normalize == 1):
        print "Normalizing DOS to 1.0 \n"
        dosGnorm=dosG
        for idos in range (0, ndos):
            dosGnorm[idos]=dosGnorm[idos]/intg

        return dosE, dosGnorm

    if(normalize==0):
        return dosE, dosG

    else:
        print " ERROR!! in CreateDOS function: wrong DOS normalization choice."
        return 




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
    
    plt.plot(dosG, dosE)

    plt.xlabel("Density Of States")
    plt.ylabel("band energies (eV)")

    plt.gca().set_xlim(left=0)
    plt.savefig(plotname)

    return




def PlotBnD(nbnd, nval, bnd, Ef, sympoints, nks_btw_sympoints, dosE, dosG, plotname):



    col = 'k'


    # Two subplots, unpack the axes array immediately
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

    for ibnd in range (0, nbnd):
        ax1.plot(bnd[ibnd].nrg,  markersize=2, linestyle='-', color=col) #marker = 'o') 

    ax1.set_title('Sharing Y axis')

    ax2.plot(dosG, dosE)

    ax2.set_xlim([0.0, 0.1])
    plt.ylim([-15.0, 20.0])
    

    #plt.subplots_adjust(left=0.0, right=0.8)
    plt.subplots_adjust(wspace = 0.0)



    plt.show()

    return







def PlotBnDD(nbnd, nval, bnd, Ef, sympoints, nks_btw_sympoints, sym_pt_dists, dosE1, dosG1, dosE2, dosG2, plotname):

    ######################################
    # Plot generation and formatting
    ######################################


    # Two subplots, unpack the axes array immediately
    gs = gridspec.GridSpec(1, 2,width_ratios=[1,4])
    f = plt.figure()
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[0])


    # Formatting
    col = 'k'

    ax1.set_xlabel("Brillouin zone path")
    ax1.xaxis.set_label_position("bottom")
    ax1.set_ylabel("E [eV]", rotation=270)
    ax1.yaxis.set_label_position("right")



    ax1.text(3.50-0.12, -12.50, 'Ge', fontsize=28)
    ###ax2.text(0.07, 18.00, 'Si', fontsize=18)

    ax2.set_xlabel("DOS \n [eV$^{-1}$]")
    ax2.xaxis.set_label_position("top")
    #ax2.set_ylabel("E [eV]", rotation=270)



    #y_min = -32.0     
    y_min = -13.0 
    y_max = 20.0

    
    x2_min = 0.00
    x2_max = 5.00
    
    # Mirror
    x2_min = 0.12
    x2_max = 0.00
    

    ax1.set_ylim([y_min, y_max])
    ax2.set_xlim([x2_min, x2_max])
    #ax2.set_xlim([0.0, 10.0])
    ax2.set_ylim([y_min, y_max])
        


    # Ticks
    #minor_locator = AutoMinorLocator(2)
    #ax2.xaxis.set_minor_locator(minor_locator)



    # Number of symmetry points
    numsympoints = len(sympoints)
    


    # Generate horizontal axis containing k-path accumulated length (for BS plot)
    x=0.0
    klen=[x]
    dx=1.0/((numsympoints-1)*nks_btw_sympoints)
    for isym in range(0, numsympoints-1):
        dx=sym_pt_dists[isym]/nks_btw_sympoints

        for ipt in range(1, nks_btw_sympoints+1):
            x=x+dx
            klen.append(x)
    

            
    #xticks = range(0, numsympoints*nks_btw_sympoints + 1, nks_btw_sympoints)
    xticks=[]
    for isym in range(0, numsympoints):
        j = isym * nks_btw_sympoints    
        xticks.append(klen[j])



    x1_min=min(xticks)
    x1_max=max(xticks)
    ax1.set_xlim(x1_min, x1_max)

    

    # Plot bands    
    col = '0.4'   
    for ibnd in range (0, nbnd):
        ax1.plot(klen , bnd[ibnd].nrg,  markersize=2, linestyle='-', color=col) #marker = 'o') 



                
#    plt.axvline(x=xticks, ymin=0, ymax=1, hold=None, **kwargs)



    # Ticks and vertical lines across BS plot
    ax1.set_xticks(xticks)
    ax1.set_xticklabels(sympoints)



    # Plot DOSs
    ax2.plot(dosG1, dosE1, linestyle='-', linewidth=1.0, color='b')
    ax2.plot(dosG2, dosE2, linestyle='-', color='r')

    

    #dosticks=[0.0, 0.05, 0.1, 0.15]
    dosticks=[5, 0] # Mirror
    ax2.set_xticks(dosticks)
    ax2.set_xticklabels(dosticks)
    
    
    #minor_locator = AutoMinorLocator(5)    
    #ax2.xaxis.set_minor_locator(minor_locator)
    minorx2ticks=[4, 3, 2, 1]
    ax2.set_xticks(minorx2ticks, minor = True)



    # BS ticks
    yticks=[-10, -5, 0, 5, 10, 15, 20]
    minor_locator = AutoMinorLocator(5)
    ax1.yaxis.set_minor_locator(minor_locator)
    ax2.yaxis.set_minor_locator(minor_locator)
    ax1.xaxis.tick_top()
    
    
    #ax1.set_yticks(yticks)
    #ax1.set_yticklabels(yticks)


    # Mirror
    ax1.yaxis.tick_right()
    ax1.set_yticks(yticks)
    ax1.set_yticklabels(yticks)
    ax2.set_yticklabels([])

    
    
    #plt.subplots_adjust(left=0.0, right=0.8)
    plt.subplots_adjust(wspace = 0.0)



   
    # Attempt to fill the area to the left of the DOS
    # split values into positive and negative

    alpha_fill=0.5

    dosE1neg=[]
    dosG1neg=[]
    dosE1pos=[]
    dosG1pos=[]
    for i in range(0, len(dosE1)):
        if(dosE1[i]<0.0):
            dosE1neg.append(dosE1[i])
            dosG1neg.append(dosG1[i])
        else:
            dosE1pos.append(dosE1[i])
            dosG1pos.append(dosG1[i])
    

    dosE1new =[y_min]+dosE1+[y_max]
    dosG1new =[0.0]+dosG1+[0.0]
    ax2.fill_between(dosG1new, 0, dosE1new, alpha=alpha_fill, linewidth=0.0, edgecolor='w')



    # Vertical lines across BS plot
    for i in range(0,numsympoints):
        ax1.axvline(x=xticks[i], ymin=y_min, ymax=y_max, color='k', linewidth=0.25)
    

    # Horizontal line at top of valence band
    if (not Ef): 
        ax1.axhline(Ef, color="black", linestyle="--")
        ax2.axhline(Ef, color="black", linestyle="--")


    #plt.show()
    plt.savefig(plotname, bbox_inches='tight')
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



bohr2ang=0.52918



############
# Band structure
############

filename="ge.bands.dat"
nks = 0
nbnd=0
xk=[]
bsflt=[]
bs=[]


sympoints=['$L$','$\Gamma$', '$X$', '$W$', '$K$', '$\Gamma$']  
sym_pt_dists=[0.5*math.sqrt(3), 1.0, 0.5, 0.25*math.sqrt(2), 0.75*math.sqrt(2)] ## distances between symmetry points (by hand)



nks_btw_sympoints=50

# Read from file and sort bs by bands
nks, nbnd, xk, bsflt = ReadBandStructure(filename)



if(nbnd==0):
    print "%% ERROR READING BANDS. EXIT %%"
else:
    bs = SortByBands(nks, nbnd, bsflt)

    print "nks=", nks
    print "nbnd=", nbnd

# Create band objects
    bands=[]
    for ibnd in range (0, nbnd):
        ledge = ibnd*nks
        redge = ledge+nks
        currentband = bs[ledge:redge]
        bands.append( band(nks, currentband) )
    


# Scissor correction
   
    # Si
    ###alat = 10.330495  # Bohr
    ###nval = 4 # for Si
    ###exptgap = 1.11 # eV # Si

    # Ge
    alat = 10.8171069 # Bohr
    nval = 14 # for Ge with semicore
    exptgap = 0.67  # Ge
    

    
    # Convert to ANG and calculate BZV
    alat=alat*bohr2ang
    V=(alat**3)/4.0 # Good for FCC 
    bzv = (2.0*math.pi)**3/V
    ncond = nbnd - nval  
    



    Scissor(nks, exptgap, bands, 1) # 3rd argument is 1. Then set 3rd argument of PlotBandStructure to 0.0
    print "Scissor correction with gap set to", exptgap




#############
# DOS
#############

filename='ge.bands_full.dat'
nks1, nbnd1, xk1, bsflt1 = ReadBandStructure(filename)


if(nbnd==0):
    print "%% ERROR READING BANDS. EXIT %%"
else:
    bs1 = SortByBands(nks1, nbnd1, bsflt1)

    print "nks=", nks1
    print "nbnd=", nbnd1


# Create band objects
    bands1=[]
    for ibnd in range (0, nbnd1):
        ledge1 = ibnd*nks1
        redge1 = ledge1+nks1
        currentband1 = bs1[ledge1:redge1]
        bands1.append( band(nks1, currentband1) )
    


# Scissor correction   
    Scissor(nks1, exptgap, bands1, 1) # 3rd argument is 1. Then set 3rd argument of PlotBandStructure to 0.0
    print "Scissor correction with gap set to", exptgap





filename='ge.bands_243.dat'
nks2, nbnd2, xk2, bsflt2 = ReadBandStructure(filename)


if(nbnd==0):
    print "%% ERROR READING BANDS. EXIT %%"
else:
    bs2 = SortByBands(nks2, nbnd2, bsflt2)

    print "nks=", nks2
    print "nbnd=", nbnd2


# Create band objects
    bands2=[]
    for ibnd in range (0, nbnd2):
        ledge2 = ibnd*nks2
        redge2 = ledge2+nks2
        currentband2 = bs2[ledge2:redge2]
        bands2.append( band(nks2, currentband2) )
    


# Scissor correction   
    Scissor(nks2, exptgap, bands2, 1) # 3rd argument is 1. Then set 3rd argument of PlotBandStructure to 0.0
    print "Scissor correction with gap set to", exptgap




# Generate DOSs
    deltaE = 0.03 #eV
    dosE1, dosG1 = CreateDOS(nks1, nbnd1, bzv, -13.0, 25.0, deltaE, bands1, 0)
    dosE2, dosG2 = CreateDOS(nks2, nbnd2, bzv, -13.0, 25.0, deltaE, bands2, 0)
    

# Plot
    #PlotDOS(dosE, dosG, "DOS.pdf")
    #PlotBandStructure(nbnd, nval, bands, "BS.pdf", 0.0, sympoints, nks_btw_sympoints)

    PlotBnDD(nbnd, nval, bands, 0.0, sympoints, nks_btw_sympoints, sym_pt_dists, dosE1, dosG1, dosE2, dosG2, "BSnDOS.pdf")





# DOS
#mydos=dos()
#mydos.Load('dos_full.dat')
#mydos.Printout()
