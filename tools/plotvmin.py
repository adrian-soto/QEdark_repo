#!/usr/local/bin/python
import math
import numpy as np
import matplotlib.pyplot as plt


def vmin(dE, mx, q):
    return dE/q + q/(2.0*mx)


def qmax(ecutwfc):
    me = 0.510998910E+06
    return 2.0*math.sqrt(2.0 * me * ecutwfc)



def format_e(n):
    #
    # Convert real numbers to scientific notation. For printing.
    # from http://stackoverflow.com/questions/6913532/python-how-to-convert-decimal-to-scientific-notation
    #
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]



#######################################################

# Parameters
c = 299792458.0
Ry2eV = 13.605698066
me = 0.510998910E+06
ecutwfc_Ry = [100.0, 250.0, 500.0]  # in Ry!!
ecutwfc=[]
for i in range (0, len(ecutwfc_Ry)):
    ecutwfc.append(ecutwfc_Ry[i] * Ry2eV)





q_max= qmax( max(ecutwfc) )



dE=[1.0, 10.0, 100.0]  
mxmin=0.10E06 
mxmax=1.0E09
q=np.arange(0.001*q_max , q_max, 0.001*q_max)

vearth = 2.40E+05/c
vesc = 6.0E+05/c
vmax = vearth+vesc


linecolors=['blue', 'red', 'green', 'cyan', 'magenta', 'yellow']
mylabel=['', '', ''] 

fig = plt.figure()

i=0
mx=mxmin
while (mx <= mxmax):
    mx = mxmin * 10**i
    
    mylabel[0] = 'mx=' + str(format_e(mx/1.0E+06)) + 'MeV' + ' , dE=', str(dE[0]) + 'eV' 
    mylabel[1] = 'mx=' + str(format_e(mx/1.0E+06)) + 'MeV' + ' , dE=', str(dE[1]) + 'eV' 
    mylabel[2] = 'mx=' + str(format_e(mx/1.0E+06)) + 'MeV' + ' , dE=', str(dE[2]) + 'eV' 
    
    
    plt.plot(q, vmin(dE[0], mx, q), color=linecolors[i], linestyle='-', label=mylabel[0])
    plt.plot(q, vmin(dE[1], mx, q), color=linecolors[i], linestyle='-.', label=mylabel[1])
    plt.plot(q, vmin(dE[2], mx, q), color=linecolors[i], linestyle=':', label=mylabel[2])

    i=i+1

# Plot boundaries
x_min=0.0
x_max=1.1*q_max
y_min=0.0
y_max=17.0*vmax
plt.xlim([x_min, x_max]) 
plt.ylim([y_min, y_max]) 


# Plot labels, legend and title
plt.xlabel('q (1/c)')
plt.ylabel('vmin (c)')
plt.title('vmin(q)')
plt.legend(loc=1,prop={'size':6})


# Subplot for parameter space boundaries
ax = fig.add_subplot(111)

# Horizontal line showing the max allowed value of v
ax.axhline(y=vmax, xmin=0.0, xmax=x_max, linewidth=2, color = 'k')
note='$v_{esc} + v_{Earth} =' + str(format_e(vmax)) + '$'
xtext = 0.7*q_max
ytext = 1.2*vmax
ax.annotate(note, xy=(xtext, ytext), xytext=(xtext, ytext), size=10)#, arrowprops=dict(facecolor='black', shrink=0.02))

# Fill kinematically excluded regions
plt.axvspan(x_min, x_max, ymin=vmax/y_max, ymax=1.0, alpha=0.70, color='grey')




# Vertical line showing the max q for a given ecutwfc. Fill excluded zone
opacity=0.10
for i in range (0, len(ecutwfc)):
    currentq_max = qmax(ecutwfc[i])
    ax.axvline(x=currentq_max, ymin=0.0, linewidth=2, color = 'k', linestyle='--')

    note= str(ecutwfc_Ry[i]) + 'Ry'
    xtext = 1.02*currentq_max
    ytext = 0.4*y_max
    ax.annotate(note, xy=(xtext, ytext), xytext=(xtext, ytext), size=10, rotation=90.0 ) #, arrowprops=dict(facecolor='black', shrink=0.02))
    
    xshade=np.arange(currentq_max, x_max, 0.001*x_max)
    plt.fill_between(xshade, 0.0, x_max-q_max, facecolor='red', alpha=opacity)


    


###yshade=np.arange(vmax, y_max, 0.001*y_max)
###plt.fill_between(yshade, y_max, vmax, where=None, facecolor='yellow', alpha=opacity)

plt.show()
