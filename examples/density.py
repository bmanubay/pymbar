
#--------------------------------------
# OpenMM Minimization, NVT, NPT, PROD
#
# Author: Dr Gaetano Calabro'
# University of California, Irvine
# ver 0.0 06/23/2016
#--------------------------------------


import numpy as np
import os, sys
import math
import simtk.openmm as mm
from simtk.openmm import app
from simtk.openmm import Platform
from simtk.unit import *

from pymbar import timeseries as ts
import pandas as pd



#----------------USER INFO-----------------
#------------------------------------------

identifier = "methanol_cyclohexane_1_800"

DATA_PATH = "./data/"
RESULT_PATH = "./results/"

prmtop_filename = DATA_PATH + "tleap/" + identifier + ".prmtop"
inpcrd_filename = DATA_PATH + "tleap/" + identifier + ".inpcrd"


#--------------MINIMIZATION----------------
MIN_TIME_STEP = 0.5 * femtoseconds
MIN_STEPS =  0 # 0=Until convergence is reached
MIN_TOLERANCE  = 10.0 * kilojoule_per_mole
MIN_PLATFORM = Platform.getPlatformByName('Reference')
MIN_FRICTION = 1.0 / picoseconds

#-------------------NVT--------------------
NVT_TIME_STEP = 0.5 * femtoseconds
NVT_STEPS = 100000
NVT_FRICTION = 1.0 / picoseconds
NVT_PLATFORM = Platform.getPlatformByName('CUDA')
NVT_PROPERTIES = {'CudaPrecision': 'mixed'}
NVT_OUTPUT_FREQ = 10000
NVT_DATA_FREQ = 10000

#------------------NPT--------------------
NPT_TIME_STEP = 1.0 * femtoseconds
NPT_STEPS = 10000000
NPT_FRICTION = 1.0 / picoseconds
BAROSTAT_FREQUENCY = 25
NPT_PLATFORM = Platform.getPlatformByName('CUDA')
NPT_PROPERTIES = {'CudaPrecision': 'mixed'}
NPT_OUTPUT_FREQ = 500000
NPT_DATA_FREQ = 500000

#--------------PRODUCTION------------------
PROD_TIME_STEP = 1.0 * femtoseconds
PROD_STEPS = 1000
PROD_FRICTION = 1.0 / picoseconds
PROD_PLATFORM = Platform.getPlatformByName('CUDA')
PROD_PROPERTIES = {'CudaPrecision': 'mixed'}
PROD_OUTPUT_FREQ = 500
PROD_DATA_FREQ = 500

#------------GEN PARAMETERS--------------
CUTOFF = 0.95 * nanometers
TEMPERATURE = 300 * kelvin
PRESSURE = 1.0 * atmospheres
#Density STD tolerance
STD_ERROR_TOLERANCE = 0.0002 # g/mL


#---------------END USER INFO---------------
#-------------------------------------------




nvt_dcd_filename = RESULT_PATH + "nvt/" + identifier + "_nvt.dcd"
nvt_data_filename = RESULT_PATH + "nvt/" + identifier + "_nvt.csv"

npt_dcd_filename = RESULT_PATH + "npt/" + identifier + "_npt.dcd"
npt_data_filename = RESULT_PATH + "npt/" + identifier + "_npt.csv"

prod_data_filename = RESULT_PATH + "prod/" + identifier + "_prod.csv"
prod_dcd_filename = RESULT_PATH + "prod/" + identifier + "_prod.dcd"



def make_path(filename):
    try:
        path = os.path.split(filename)[0]
        os.makedirs(path)
    except OSError:
        pass


def minimization():
    # Load Amber files 
    prmtop = app.AmberPrmtopFile(prmtop_filename)
    inpcrd = app.AmberInpcrdFile(inpcrd_filename)
        
    # Create OpenMM System
    system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=CUTOFF, constraints=app.HBonds)
        
    # Select Integrator
    integrator = mm.LangevinIntegrator(TEMPERATURE, MIN_FRICTION, MIN_TIME_STEP)
        
    # Set Simulation
    simulation = app.Simulation(prmtop.topology, system, integrator, MIN_PLATFORM)

    # Set Position
    simulation.context.setPositions(inpcrd.positions)

    state = simulation.context.getState(getEnergy=True)

    if np.isnan(state.getPotentialEnergy() / kilojoule_per_mole):
        raise ValueError("The Potential Energy before minimization is NaN")

    # Minimization
    print('Minimizing...\n')
    simulation.minimizeEnergy(tolerance=MIN_TOLERANCE, maxIterations=MIN_STEPS)

    state = simulation.context.getState(getPositions=True, getEnergy=True)

    if np.isnan(state.getPotentialEnergy() / kilojoule_per_mole):
        raise ValueError("The Potential Energy after minimization is NaN")

    coords = state.getPositions()
        
    return coords
        

def nvt(coords):
    # Load Amber files 
    prmtop = app.AmberPrmtopFile(prmtop_filename)
    
    # Create OpenMM System
    system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=CUTOFF, constraints=app.HBonds)
        
    # Select Integrator
    integrator = mm.LangevinIntegrator(TEMPERATURE, NVT_FRICTION, NVT_TIME_STEP)
        
    # Set Simulation
    simulation = app.Simulation(prmtop.topology, system, integrator, NVT_PLATFORM, NVT_PROPERTIES)

    # Set Position and velocities
    simulation.context.setPositions(coords) 
    simulation.context.setVelocitiesToTemperature(TEMPERATURE)
        
    # Set Reporter 
    simulation.reporters.append(app.DCDReporter(nvt_dcd_filename, NVT_OUTPUT_FREQ))
    simulation.reporters.append(app.StateDataReporter(nvt_data_filename, NVT_DATA_FREQ, step=True, potentialEnergy=True, temperature=True, density=True))
         
    state = simulation.context.getState(getEnergy=True)
    
    if np.isnan(state.getPotentialEnergy() / kilojoule_per_mole):
        raise ValueError("The Potential Energy before NVT is NaN")

    print('NVT...\n')
    simulation.step(NVT_STEPS)

    state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True)

    if np.isnan(state.getPotentialEnergy() / kilojoule_per_mole):
        raise ValueError("The Potential Energy after NVT is NaN")


    coords = state.getPositions()
    velocities = state.getVelocities()

    return coords, velocities


        
def npt(coords, velocities):
    # Load Amber files 
    prmtop = app.AmberPrmtopFile(prmtop_filename)
        
    # Create OpenMM System
    system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=CUTOFF, constraints=app.HBonds)
        
    # Select Integrator
    integrator = mm.LangevinIntegrator(TEMPERATURE, NPT_FRICTION, NPT_TIME_STEP)

    # Set Barostat
    system.addForce(mm.MonteCarloBarostat(PRESSURE, TEMPERATURE, BAROSTAT_FREQUENCY))
    
    # Set Simulation
    simulation = app.Simulation(prmtop.topology, system, integrator, NPT_PLATFORM, NPT_PROPERTIES)

    # Set Position and velocities
    simulation.context.setPositions(coords) 
    simulation.context.setVelocities(velocities)
        
    # Set Reporter 
    simulation.reporters.append(app.DCDReporter(npt_dcd_filename, NPT_OUTPUT_FREQ))
    simulation.reporters.append(app.StateDataReporter(npt_data_filename, NPT_DATA_FREQ, step=True, potentialEnergy=True, temperature=True, density=True))
            
    state = simulation.context.getState(getEnergy=True)
    
    if np.isnan(state.getPotentialEnergy() / kilojoule_per_mole):
        raise ValueError("The Potential Energy before NPT is NaN")

    print('NPT...\n')
    simulation.step(NPT_STEPS)

    state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
    
    if np.isnan(state.getPotentialEnergy() / kilojoule_per_mole):
        raise ValueError("The Potential Energy after NPT is NaN")


    coords = state.getPositions()
    velocities = state.getVelocities()
    box = state.getPeriodicBoxVectors()

    return coords, velocities, box

        

def production(coords, velocities, box):
    # Load Amber files 
    prmtop = app.AmberPrmtopFile(prmtop_filename)
        
    # Create OpenMM System
    system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=CUTOFF, constraints=app.HBonds)
        
    # Select Integrator
    integrator = mm.LangevinIntegrator(TEMPERATURE, PROD_FRICTION, PROD_TIME_STEP)

    # Set Barostat
    system.addForce(mm.MonteCarloBarostat(PRESSURE, TEMPERATURE, BAROSTAT_FREQUENCY))
    
    # Set Simulation
    simulation = app.Simulation(prmtop.topology, system, integrator, PROD_PLATFORM, PROD_PROPERTIES)

    # Set Position and velocities
    simulation.context.setPositions(coords) 
    simulation.context.setVelocities(velocities)

    # Set Box
    #box = box.in_units_of(nanometer)
    simulation.context.setPeriodicBoxVectors(box[0], box[1], box[2])

        
    # Set Reporter 
    simulation.reporters.append(app.DCDReporter(prod_dcd_filename, PROD_OUTPUT_FREQ))
    simulation.reporters.append(app.StateDataReporter(prod_data_filename, PROD_DATA_FREQ, step=True, potentialEnergy=True, temperature=True, density=True))
      
      
    state = simulation.context.getState(getEnergy=True)
    
    if np.isnan(state.getPotentialEnergy() / kilojoule_per_mole):
        raise ValueError("The Potential Energy before Production is NaN")

    print('PRODUCTION...\n')
    
    converged = False
    
    while not converged:
        
        simulation.step(PROD_STEPS)
        
        d = pd.read_csv(prod_data_filename, names=["step", "U", "Temperature", "Density"], skiprows=1)
        density_ts = np.array(d.Density)
        [t0, g, Neff] = ts.detectEquilibration(density_ts, nskip=1000)
        density_ts = density_ts[t0:]
        density_mean_stderr = density_ts.std() / np.sqrt(Neff)
                
        print("Current density mean std error = %f g/mL" % density_mean_stderr)
        
        if density_mean_stderr < STD_ERROR_TOLERANCE:
            converged = True
            print("...Convergence is OK\n")

    state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True)

    
    if np.isnan(state.getPotentialEnergy() / kilojoule_per_mole):
        raise ValueError("The Potential Energy after Production is NaN")


    coords = state.getPositions()
    velocities = state.getVelocities()
    box = state.getPeriodicBoxVectors()

    return coords, velocities, box


if __name__=='__main__':

    make_path(RESULT_PATH + 'nvt/')
    make_path(RESULT_PATH + 'npt/')
    make_path(RESULT_PATH + 'prod/')

    
    coords = minimization()
    coords, velocities = nvt(coords)
    coords, velocities, box = npt(coords, velocities)
    production(coords, velocities, box)
