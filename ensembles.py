from numpy import *
from numpy.random import rand
from math import pi

# Ensembles hold information about
#  * Potential energy
#  * Kinetic energy
#  * The Virial
#  * Temperature
#
# It is the job of the solvers to update the relevant
# properties via the getX and setX methods.

# An ensemble.
class Ensemble:
  # Creates the class instance. Takes a geometry and
  # a temperature as parameters.
  def __init__(self,geo,temp=0.0):
    self.geometry = geo
    self.T = temp

    self.Ekin = 0.0
    self.Epot = 0.0
    self.Vir  = 0.0
    
    # Calculate random velocities
    #Make sure that total velocity is zero.
    (ndim,n) = shape(self.geometry.getPos())
    V = rand( ndim,n ) - 0.5
    Vm = sum(V,1)/n
    V -= reshape(Vm, (ndim,1))

    #Scale the velocity with temperature
    fs = dot( ravel(V), ravel(V) ) / n
    fs = sqrt(ndim * self.T/fs)
    V *= fs

    # set the velocity
    self.setVel(V)

  # Sets the velocities of all the particles to
  # newvel and calculates the temperature and
  # kinetic energy while doing so.
  def setVel(self,newvel):
    (ndim,n) = shape(newvel)
    self.V = newvel
    self.Ekin = 0.0
    self.T = 0.0
    
    #self.getVel() returns [[vx1, vx2, ...], [vy1, vy2, ...], [vz1, vz2, ...]]
    #The kinetic energy is calculated by summing over all products
    #of velocity-vector components:
    #E_kin(t) = vx1*vx1 + vx2*vx2 + ... + vxn*vxn + vy1*vy1 + vy2*vy2 + ...
    #           ... + vyn*vyn + vz1*vz1 + vz2*vz2 + ... + vzn*vzn
    self.Ekin = 0.5*dot(ravel(self.getVel()), ravel(self.getVel()))

    #T = 2*E_kin/(N_dim*N_part)
    #
    #>>> A = np.array([[1,2,3], [4,5,6]])
    #>>> np.shape(A)
    #(2, 3)
    #
    #newvel = self.getVel()
    #shape(self.getVel())[0] = N_dim
    #shape(self.getVel())[1] = N_part
    self.T    = 2*self.Ekin/(shape(newvel)[0]*shape(newvel)[1]) 


  # Returns the velocity of all the particles
  def getVel(self):
    return self.V

  def getRho(self):
    return self.geometry.getRho()

  # Sets the positions of all the particles.
  # NOTE: it actually updates the geometry class
  # since it is where we store the positions
  def setPos(self,newpos):
    self.geometry.setPos(newpos)

  # Returns the particle positions from the geometry class
  def getPos(self):
    return self.geometry.getPos()

  # Returns the box sizes
  def getBox(self):
    return self.geometry.getBox()

  # Returns the number of particles
  def getN(self):
    (ndim,n) = shape(self.geometry.getPos())
    return n

  # Returns the kinetic energy
  def getEpot(self):
    return self.Epot

  # Sets the kinetic energy
  def setEpot(self,epot):
    self.Epot = epot

  # Sets the kinetic energy
  def getEkin(self):
    return self.Ekin

  # Gets the temperature
  def getT(self):
    return self.T

  # Sets the temperature. Remember to scale it!
  def setT(self, temp):
    V = array(self.V)
    (ndim,n) = shape(V)

    #Reset velocity to zero.
    Vm = sum(V,1)/n
    V -= reshape(Vm, (ndim,1))

    #Rescaling temperature
    fs = dot( ravel(V), ravel(V) ) / n
    fs = sqrt(ndim * temp/fs)
    V *= fs

    #E_kin
    V2 = dot( ravel(V), ravel(V) )

    self.Ekin = 0.5*V2
    self.T = temp
    self.V = V

  # Sets the virial
  def setVir(self, vir):
    self.Vir = vir

  # Gets the virial
  def getVir(self):
    return self.Vir

  # Get the pressure
  def getP(self):
    rc = 2.5  # The cutoff distance
    pressure = 0.0
    correction = 0.0

    #P(t) = rho*T(t) + 1/(3*V)*Sum_(i,j)(F(r(i,j))*r(i,j)
    #     = rho*T(t) + 1/(3*V)*K_vir
    #
    #DeltaP^tail = 16/3*pi*rho^2*(2/3*(1/r_c)^9-(1/r_c)^3) 
    pressure = self.getRho()*self.getT() + 1.0/(3*self.getBox()[0]**3)*self.getVir()
    correction = 16.0/3.0*pi**2*(2.0/3.0*(1.0/rc)**9 - (1.0/rc)**3) 
    return pressure + correction








