import pylab
import numpy as numpy
import math
import copy

# import classes to make MD
from geometries import CubicLattice
from ensembles import Ensemble
from solvers import VeloVerlet
from forces import lennardjones
from sampler import Sampler
from analyse import Analyse
from diffusion import Diffusion

class MDSim:
    def __init__(self, numpyart, rho, temp, force, stepsize=0.001):
        """Initializes the various parts of an MD simulation.

        Parameters:
        numpyart: Number of particles
        rho  : Density
        temp : Temperature
        force: A force model
        stepsize: Time step
        """ 
        self.geometry = CubicLattice( numpyart, rho )               # Initialize a geometry,
                                                                 # in this case a cubic lattice.
        self.system = Ensemble( self.geometry, temp )            # Create an ensemble to
                                                                 # store information about
                                                                 # the simulation.
        self.solver = VeloVerlet( self.system, force, stepsize ) # Create the sampler class.
        self.sampler = Sampler()
    
    def run( self, nsteps, tSet):
        """Simulates N steps of the simulation.""" 
        # Equilibration Run.
        print "equilibrating %i steps" % (nsteps/2)
        for n in range(0,nsteps/2,2):
            self.solver.step()
            if( n % 1000 == 0 ):
                #self.system.setT(tSet)
                print "step = %05i, temp = %12.7f" % (n, self.system.getT()) 
        print "Equilibration done." 
        # Production Run.  
        for n in range(nsteps/2, nsteps+1): 
            # sample every 100 steps but only after nsteps / 2 equilibration steps
            if( n % fil == 0 and n > nsteps / 2 ):
                U = copy.deepcopy(self.system.getPos()) 
                #Sampling
                self.sampler.sampleData( U[0], U[1], U[2],
                                         self.system.getEpot(),
                                         self.system.getEkin(),
                                         self.system.getT(),
                                         self.system.getP(),
                                         copy.deepcopy(self.system.getVel()) ) 
            # print every 500 steps
            if( n % 1000 == 0 ):
                print "step = %05i, temp = %12.7f" % (n, self.system.getT()) 
            self.solver.step()  # take a simulation step

def velDistFigure(lists):
    # Overlay several MB-distributions within a single figure.
    pylab.clf()
    for list in lists:
        pylab.hist(list, bins=40)
    pylab.savefig('velDistLists.png') 

# run the simulation
if( __name__ == '__main__' ):

  numpyart  = 108 # The total number of particles
  R0     = [0.7]#, 0.4, 0.9]     # Initial density
  fil = 50    # Sampling every 'filter' steps.
  T     = [0.5]#, 0.8, 1.0, 1.2]
  dt    = 0.001    # Stepsize
  stepnumber = 20000

  velDistributions = [] 
  for t in T:
      for r in R0:
          print "*"*50, 'r=', r, 'T=', t
          print 'r=', r, 'T=', t
          simulation = MDSim( numpyart, r, t, lennardjones, dt ) 
          simulation.run( stepnumber, t) 
          diff = Diffusion(simulation, fil, dt) 
          diff.calcMSD()
          diff.calcVACF() 

          #diff.plotMSD()
          #diff.plotVACF() 

          diff.output() 
          print 'diff.msd', diff.msd 
          print "diff.vacf", diff.vacf 

          linreg = diff.linReg()
          int = diff.integrate() 

          print "diff.msd", linreg[0]/6.0
          print "diff.vacf", int/3.0 

          #diff.writeData()
          #analyse = Analyse(simulation)
          #analyse.figures()
          #print "Number of samples taken:", len(simulation.sampler.getData())
          #analyse.showData()
          #velDistributions.append(analyse.getVelocityDistribution())

