# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 14:58:26 2016

@author: bmanubay
"""

# todo -- simplify the total energy read in, the kinetic energy read-in, temperature read-in
#=========================================================
#IMPORTS
#=========================================================
import sys
import numpy
import pandas as pd
import pylab as pl
import pymbar # for MBAR analysis
from pymbar import timeseries # for timeseries analysis
import os
import os.path
import optparse
from optparse import OptionParser


#===================================================================================================
# INPUT PARAMETERS
#===================================================================================================
parser = OptionParser()

parser.add_option("-d", "--directory", dest="simulation",
                  help="the directory of the energies we care about")
parser.add_option("-b", "--nbootstraps", dest="nBoots",type = int, default=0,
                  help="Number of samples at the beginning that are ignored") 
parser.add_option("-s", "--spacing", dest="NumIntermediates",type = "int", default = 100,
                  help="Number of intermediate simulations used to calculate finite differences (default 100)")
parser.add_option("-f", "--finitedifftype", dest="dtype", default="temperature",
                  help="the type of finite difference energy, choice is \"temperature\" or \"beta\" [default = %default]")
(options, args) = parser.parse_args()

simulation = options.simulation
nBoots = options.nBoots 
NumIntermediates = options.NumIntermediates
dtype = options.dtype

#========================================================
# CONSTANTS
#=========================================================

kB = 0.008314462  #Boltzmann constant (Gas constant) in kJ/(mol*K)
TE_COL_NUM = 11  #The column number of the total energy in ener_box#.output

NumTemps = 16  	# Last TEMP # + 1 (start counting at 1)
NumIterations = 1000  # The number of energies to be taken and analyzed, starting from the last
			# Extra data will be ignored
NumIntermediates = 200  # Number of temperatures over which to calculate properties

if (dtype == 'temperature'):
	types = ['var','dT','ddT']
elif (dtype == 'beta'):
	types = ['var','dbeta','ddbeta']
else:
	print 'type of finite difference not recognized must be \'beta\' or \'temperature\''
	quit()
ntypes = len(types)

###########################################################
# For Cv vs T    _____
#               /     \
#  ____________/       \____________
#
############################################################

#=========================================================
# SUBROUTINES
#=========================================================

def read_col(filename,colname,numsteps):
	"""Reads in the TEMP#/ener_box#.output file and parses it, returning an array of energies 
	ARGUMENTS
		filename (string) - the path to the folder of the simulation
	"""

	print "--Reading total energies from %s/..." % filename


	# Read in file output as pandas df
	df = pd.read_csv(filename, sep= ',')
	
	# Read values direct from column into numpy array
        E_kn = df.as_matrix(columns = colname)
        E_kn = E_kn[-numsteps:]


	return E_kn	

def read_simulation_temps(pathname,NumTemps):
	"""Reads in the various temperatures from each TEMP#/simul.output file by knowing
		beforehand the total number of temperatures (parameter at top)
	"""
	
	print "--Reading temperatures from %s/..." % filename	
	
	# Read in file output as pandas df
	df = pd.read_csv(filename, sep = ',')
	
	# Read values direct from column into numpy array
        temps_from_file = df.colname.values 


	return temps_from_file

def PrintResults(string,E,dE,Cv,dCv,types):

	print string
	print "Temperature    dA        <E> +/- d<E>  ",
	for t in types:
		print "    Cv +/- dCv (%s)" % (t),    
	print ""	
	print "------------------------------------------------------------------------------------------------------"
	for k in range(originalK,K):
		print "%8.3f %8.3f %9.3f +/- %5.3f" % (Temp_k[k],mbar.f_k[k]/beta_k[k],E[k],dE[k]),
		for i in range(len(types)):
			if Cv[k,i,0] < -100000.0:
				print "         N/A          ", 
			else:
				print "    %7.4f +/- %6.4f" % (Cv[k,i,0],dCv[k,i]), 
		print ""	
			
#========================================================================
# MAIN
#========================================================================

#------------------------------------------------------------------------
# Read Data From File
#------------------------------------------------------------------------
filename = "wm10000.csv"
print("")
print("Preparing data:")
E_from_file = -0.001 * read_col(filename,["Total Energy (kJ/mole)"],200)
T_from_file = read_col(filename,["Temperature (K)"],200)
Vol_from_file = read_col(filename,["Box Volume (nm^3)"],200)
Dens_from_file = read_col(filename,["Density (g/mL)"],200)

NumIntermediates = len(T_from_file) + 1  # Number of temperatures over which to calculate properties
K = len(T_from_file)
N_k = numpy.zeros(K,numpy.int32)
g = numpy.zeros(K,numpy.float64)

# Consult MRS about issues with subsampling energies. Covariance == 0 cannot compute statistical inefficiency 

for k in range(K):  # subsample the energies
   g[k] = 1
   indices = numpy.array(timeseries.subsampleCorrelatedData(E_from_file[k],g=g[k])) # indices of uncorrelated samples
   N_k[k] = len(indices) # number of uncorrelated samples

print(N_k)
print(len(N_k))

# Plot energy vs time to look for reasons why covariance might be zero
t_chk = 2. * read_col(filename,["Step"],450000) # in fs
E_chk = -0.001 * read_col(filename,["Total Energy (kJ/mole)"],450000) 

pl.figure()
pl.plot(t_chk, E_chk, color = 'red', linestyle = '-', linewidth = 0.4)
pl.xlabel('Time (fs)')
pl.ylabel('Total Energy (kJ/mole)')
pl.ylim(24.5,29)
pl.show()

#------------------------------------------------------------------------
# Insert Intermediate T's and corresponding blank U's and E's
#------------------------------------------------------------------------
Temp_k = T_from_file
minT = T_from_file[0]
maxT = T_from_file[len(T_from_file) - 1]
#beta = 1/(kB*T)
#T = 1/(kB*beta)
if dtype == 'temperature':
	minv = minT
	maxv = maxT
elif dtype == 'beta':   # actually going in the opposite direction as beta for logistical reasons
	minv = 1/(kB*minT)
	maxv = 1/(kB*maxT)
delta = (maxv-minv)/(NumIntermediates-1)
originalK = len(Temp_k)

print("--Adding intermediate temperatures...")

val_k = []
currentv = minv
if dtype == 'temperature':	
	# Loop, inserting equally spaced T's at which we are interested in the properties
	while (currentv <= maxv) :
		val_k = numpy.append(val_k, currentv)
		currentv = currentv + delta
	val_k = numpy.reshape(val_k, (len(val_k),1))
	Temp_k = numpy.concatenate((Temp_k,numpy.array(val_k)))
elif dtype == 'beta':
	# Loop, inserting equally spaced T's at which we are interested in the properties
	while (currentv >= maxv) :
		val_k = numpy.append(val_k, currentv)
		currentv = currentv + delta
	val_k = numpy.reshape(val_k, (len(val_k),1))
	Temp_k = numpy.concatenate((Temp_k,(1/(kB*numpy.array(val_k)))))

# Update number of states
K = len(Temp_k)
#print(Temp_k)
#print(K)

# Loop, inserting E's into blank matrix (leaving blanks only where new Ts are inserted)
Nall_k = numpy.zeros([K], numpy.int32) # Number of samples (n) for each state (k) = number of iterations/energies
E_kn = numpy.zeros([K, NumIterations], numpy.float64)

for k in range(originalK):
	E_kn[k,0:N_k[k]] = E_from_file[k,0:N_k[k]]
	Nall_k[k] = N_k[k]

#------------------------------------------------------------------------
# Compute inverse temperatures
#------------------------------------------------------------------------
beta_k = 1 / (kB * Temp_k)

#------------------------------------------------------------------------
# Compute reduced potential energies
#------------------------------------------------------------------------

print "--Computing reduced energies..."

u_kln = numpy.zeros([K,K,NumIterations], numpy.float64) # u_kln is reduced pot. ener. of segment n of temp k evaluated at temp l
E_kn_samp = numpy.zeros([K,NumIterations], numpy.float64) # u_kln is reduced pot. ener. of segment n of temp k evaluated at temp l

nBoots_work = nBoots + 1 # we add +1 to the bootstrap number, as the zeroth bootstrap sample is the original

allCv_expect = numpy.zeros([K,ntypes,nBoots_work], numpy.float64)
dCv_expect = numpy.zeros([K,ntypes],numpy.float64)
allE_expect = numpy.zeros([K,nBoots_work], numpy.float64)
allE2_expect = numpy.zeros([K,nBoots_work], numpy.float64)
dE_expect = numpy.zeros([K],numpy.float64)

for n in range(nBoots_work):
	if (n > 0):
		print "Bootstrap: %d/%d" % (n,nBoots)
	for k in range(K):
	# resample the results:
		if Nall_k[k] > 0:
			if (n == 0):  # don't randomize the first one
				booti = numpy.array(range(N_k[k]))
			else:
				booti=numpy.random.randint(Nall_k[k],size=Nall_k[k])
			E_kn_samp[k,0:Nall_k[k]] = E_kn[k,booti] 

	for k in range(K):
		for l in range(K):
			u_kln[k,l,0:Nall_k[k]] = beta_k[l] * E_kn_samp[k,0:Nall_k[k]]
   
#------------------------------------------------------------------------
# Initialize MBAR
#------------------------------------------------------------------------

# Initialize MBAR with Newton-Raphson
	if (n==0):  # only print this information the first time		
		print ""
		print "Initializing MBAR:"
		print "--K = number of Temperatures with data = %d" % (originalK)
		print "--L = number of total Temperatures = %d" % (K) 
		print "--N = number of Energies per Temperature = %d" % (numpy.max(Nall_k))

        # Use Adaptive Method (Both Newton-Raphson and Self-Consistent, testing which is better)
	if (n==0):
		initial_f_k = None # start from zero 
	else:
		initial_f_k = mbar.f_k # start from the previous final free energies to speed convergence
		
	mbar = pymbar.MBAR(u_kln, Nall_k, verbose=False, relative_tolerance=1e-12, initial_f_k=initial_f_k)   