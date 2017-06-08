# Python image-generating code...
# The user must correctly specify the parameters (below). The code then
# searches the current working directory for files ending in ".ascii".
# Each file is read in, and a PNG plot file is produced.
# JMP 23/05/17 - Working!
# To run: python image.py

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib
import numpy as np
import os
import glob

problem = "CWB"    # problem name
geom = "ZRP"       # geometry
nd   = 2           # dimensionality of grid
imax = 200          # number of real zones in 1st direction
jmax = 200          # number of real zones in 2nd direction (1 if 1D)
kmax = 1           # number of real zones in 3rd direction (1 if 2D)
smax = 1           # total number of advected scalars (min=1)
cmax = 1           # number of advected colours (min=1)
iqmax = 2+nd+smax  # hydro: density, pressure, nd velocity, smax
nghost = 2         # number of ghost zones at each boundary, surounding the grid

# Total number of zones in each direction (real plus ghost). These must be set by hand!
# Default values are 1
imaxd = imax + 2*nghost
jmaxd = jmax + 2*nghost
#jmaxd = 1
kmaxd = 1

# Set range for real cells on grid (irs->ire-1 inclusive, so can use range(irs,ire)
# These must also be set by hand!
irs = nghost
ire = imax+nghost
jrs = nghost
jre = jmax+nghost
krs = 0
kre = 1

# Global grid extent
xming = -1.0e12
xmaxg = 1.0e12
yming = 0.0
ymaxg = 2.0e12
zming = 0.0
zmaxg = 1.0

# Local grid extent
xmin = xming
xmax = xmaxg
ymin = yming
ymax = ymaxg
zmin = zming
zmax = zmaxg
 
# Grid arrays
P0 = np.zeros((iqmax,kmaxd,jmaxd,imaxd),dtype=float)
zxa = np.zeros((imax+1+2*nghost),dtype=float)
zxc = np.zeros((imax+1+2*nghost),dtype=float)
zxg = np.zeros((imax+1+2*nghost),dtype=float)
zdx = np.zeros((imax+1+2*nghost),dtype=float)
zya = np.zeros((jmax+1+2*nghost),dtype=float)
zyc = np.zeros((jmax+1+2*nghost),dtype=float)
zyg = np.zeros((jmax+1+2*nghost),dtype=float)
zdy = np.zeros((jmax+1+2*nghost),dtype=float)
zza = np.zeros((kmax+1+2*nghost),dtype=float)
zzc = np.zeros((kmax+1+2*nghost),dtype=float)
zzg = np.zeros((kmax+1+2*nghost),dtype=float)
zdz = np.zeros((kmax+1+2*nghost),dtype=float)

# Array variable indices
iqd = 0  # density
iqu0 = 1 # velocity0/mtm0 (in 2D, the velocity components are iqu0 and iqu0+1)
iqe = iqu0+nd # pressure/energy
iqal0 = iqe+1 # advected scalar 0 (with smax=2, the scalar components are iqal0 and iqal0+1)

# Global variables
nfile2d = 0
t = 0.0

def readIn2D(filename):
    """Read a 2D data file."""

    global t
    print("Reading in file: ",filename)
    
    #Read data line-by-line from an ascii .txt file
    with open(filename, 'rt') as input_file:
        header = input_file.readline()
        #print(header)
        words = header.split()
        t = float(words[-1])
        #print(t)
        header = input_file.readline()
        #print(header)
        for j in range(jrs,jre):
            for i in range(irs,ire):
                line = input_file.readline()
                #print(line,j,jrs,jre,i,irs,ire)
                floats = [float(x) for x in line.split()]
                zxc[i] = floats[0]
                zyc[j] = floats[1]
                P0[iqd][0][j][i] = floats[2]
                P0[iqe][0][j][i] = floats[3]
                P0[iqu0][0][j][i] = floats[4]
                P0[iqal0][0][j][i] = floats[5]
            pad = input_file.readline()
    
    
    return


def plot2DData():
    """Plot the data as a 2D image."""

    # 2D SOD
    #qmin = 0.0
    #qmax = 15.0

    # 2D SNR
    #qmin = 0.1*gVar.pamb
    #qmax = 1.0e5*gVar.pamb

    # 2D WBB
    #qmin = 0.01*gVar.pamb
    #qmax = 1.0e5*gVar.pamb

    # 2D CWB
    qmin = 1.0e-15
    qmax = 1.0e-10

    filename = problem.lower() + "_" + str(nd) + "d_rho_" + str(nfile2d).zfill(4)
    #filename = problem.lower() + "_" + str(nd) + "d_pre_" + str(nfile2d).zfill(4)

    # It is not possible to use a Python view to extract the data without the ghost zones,
    # because the plot function needs contiguous data. Therefore we have to copy the data
    # in the real cells to a temporary array
    data2d = np.zeros((jmax,imax),dtype=float)
    for j in range(0,jmax):
        jg = j+nghost
        for i in range(0,imax):
            ig = i+nghost
            data2d[j,i] = P0[iqd][0][jg][ig]
            
    #print("plot2DData: P0.size = ",P0.shape)
    #print("plot2DData: data2d.size = ",data2d.shape)
    #sys.exit(1)
    
    plt.figure(figsize=(2,2))
    matplotlib.rcParams.update({'font.size': 4})
    #matplotlib.rc('font', size=4) # this works too...
    #imgplot = plt.imshow(data2d, clim=(qmin,qmax), origin='lower', extent=[0,1,0,1], interpolation='none')
    imgplot = plt.imshow(np.log10(data2d), clim=(np.log10(qmin),np.log10(qmax)), origin='lower', interpolation='none')
    imgplot.set_cmap('nipy_spectral')
    plt.colorbar()
    #plt.title('SOD 2D XY - Pressure - t = {:.2e}'.format(t))    
    #plt.title('WBB 2D ZRP - Pressure - t = {:.2e}'.format(t))
    plt.title('CWB 2D ZRP - Density - t = {:.2e}'.format(t))    
    if geom == "ZRP":
        plt.xlabel('Z')
        plt.ylabel('R')
    elif geom == "XYZ":
        plt.xlabel('X')
        plt.ylabel('Y')        
    #plt.show()
    plt.savefig(filename+".png",dpi=300,bbox_inches="tight") # 4x4 inches
    plt.close('all')
    return


def main():
    """Main program..."""
    cwd = os.getcwd()    
    for filename in os.listdir(os.getcwd()):
        if filename.endswith(".ascii"):
            readIn2D(filename)
            plot2DData()
            global nfile2d
            nfile2d = nfile2d + 1

main()
