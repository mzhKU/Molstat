from numpy import *
from numpy.random import rand,randn,normal
import copy
from math import sqrt

# Solvers solve the equation of motion for ensembles.
# They contain information about time and they update
# the ensembles physical properties. This file contains
# 1 solver.

# VeloVerlet
#   * Solves the equations of motion using the VeloVerlet algorithm
#   * This is an NVE solver

class VeloVerlet:
  # Creates the class instance. Takes and ensemble, a force method
  # and a time interval as parameters.
  def __init__(self,ensemble,force,dt=0.001):
    self.ensemble = ensemble
    self.force  = force
    self.dt = dt

    U = self.ensemble.getPos()
    box = self.ensemble.getBox()

    (epot,F,Vir) = self.force(U,box)

    self.F = F
    self.ensemble.setEpot(epot)
    self.ensemble.setVir(Vir)

  # Makes the actual VeloVerlet step and updates
  # the ensemble
  def step(self):
    U = self.ensemble.getPos()
    V = self.ensemble.getVel()
    box = self.ensemble.getBox()
    (ndim,n) = shape(U)

    F = self.F
    Dt = self.dt

    U += V*Dt + 0.5*F*Dt*Dt
    F0 = F
    (epot,F,Vir) = self.force(U,box)
    V += 0.5*(F + F0)*Dt
    self.F = F

    self.ensemble.setPos(U)
    self.ensemble.setVel(V)
    self.ensemble.setEpot(epot)
    self.ensemble.setVir(Vir)
