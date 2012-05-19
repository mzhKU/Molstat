from numpy import *

# Geometries define starting geometries and
# keeps positions stored for time > 0

# A simple cubic lattice
class CubicLattice:
  # Creates the class instance. Takes the total number
  # of particles and arranges them in a cubic lattice.
  def __init__(self,natms=10,rho=1.0):

    n = natms
    p = int(floor(n**(1.0/3.0))) + 1

    box = (n/rho)**(1.0/3.0)
    du = box/p
    
    X = multiply.outer(ones(p),arange(p))
    X = multiply.outer(ones(p),X)
    
    Y = multiply.outer(arange(p),ones(p))
    Y = multiply.outer(ones(p),Y)

    Z = multiply.outer(arange(p),ones(p))
    Z = multiply.outer(Z,ones(p))

    U = reshape(ravel((X,Y,Z)), (3,p**3))
    U = U[:,:n]*du + du/2

    self.pos = U
    self.box = array([box]*3)
    self.rho = rho

  def getRho(self):
    return self.rho

  # Returns the box size in all dimensions
  def getBox(self):
    return self.box

  # Returns the positions of the particles
  # inside the box.
  def getBoxPos(self):
    X = self.U[0]; bx = self.Box[0]
    Y = self.U[1]; by = self.Box[1]
    Z = self.U[2]; bz = self.Box[2]

    X = where(X>=0.0, fmod(X,bx), fmod(X,bx)+bx)
    Y = where(Y>=0.0, fmod(Y,by), fmod(Y,by)+by)
    Z = where(Z>=0.0, fmod(Z,bz), fmod(Z,bz)+bz)
    U = array((X,Y,Z))

    return U


  # Returns the positions of all the particles
  def getPos(self):
    return self.pos

  # Sets the positions of all the particles
  def setPos(self,newpos):
    self.pos = newpos
