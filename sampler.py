class Sampler:
  """Sampler class to sample data during MD simulations
  """
  def __init__(self):
    """Initialize variables needed to store data
    """
    self.count = 0
    self.data = []

  def sampleData(self, X, Y, Z, ePot, eKin, temp, pre, vel):
    """Actual sampling takes place here. Just store as
       consecutive elements in a list
    """ 
    self.data.append( [X, Y, Z, ePot, eKin, temp, pre, vel] )

  def getData(self):
    """Return the entire data set
    """
    return self.data
  
