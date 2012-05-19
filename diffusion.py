import numpy
import pylab

# The module for calculating the diffusion coefficient.

# Initializing a variable to store the MSD values

class Diffusion:
    
    def __init__(self, simulation, fil, dt):
        self.data  = simulation.sampler.getData() 
        self.nsamples = len(self.data)
        print 'self.nsamples', self.nsamples
        self.positions = self.getPositions()
        self.velocities = numpy.array(self.getVelocities())

        # Number of correlation steps
        self.ncs = self.nsamples/10
        print 'self.ncs', self.ncs
        #print 'ncs', self.ncs
        self.ito = self.ncs/10
        print 'self.ito', self.ito
        #print 'ito', self.ito
        self.nto = (self.nsamples - self.ncs)/self.ito
        print 'self.nto', self.nto
        #print 'nto', self.nto
        
        # Required for the determination of the time.
        self.fil = fil
        self.dt = dt
        self.time = self.getTime()

        self.npart= len(self.positions[0][0])
        self.msd  = numpy.zeros(self.ncs, float)
        self.vacf = numpy.zeros(self.ncs, float)

    def getPositions(self):
        positions = []
        for sample in self.data:
            positions.append([sample[0], sample[1], sample[2]])
        return positions 

    def getVelocities(self):
        velocities = []
        for sample in self.data:
            velocities.append(sample[7])
        return velocities

    def calcMSD(self): 
        for n in range(self.nto):
            # Auxiliary index.
            index = 0
            # Setting the time origin.
            t0 = n*self.ito
            # Iterating over correlation steps.
            for t in range(t0, t0 + self.ncs):

                # Iterating over particles indexes.
                distT = 0.0
                for i in range(self.npart):
                    # Distance of the particle i at time t
                    # from the position of the particle i at
                    # time 0.
                    #disti = numpy.sqrt(
                    distT += numpy.sqrt(
                              (self.positions[t][0][i] - self.positions[t0][0][i])**2
                            + (self.positions[t][1][i] - self.positions[t0][1][i])**2
                            + (self.positions[t][2][i] - self.positions[t0][2][i])**2 )**2
                    #self.msd[index] += disti**2
                self.msd[index] += distT/self.npart
                index += 1
        #self.msd /= (self.npart*self.nto)
        self.msd /= (self.nto)

    def calcVACF(self):
        for n in range(self.nto):
            # Auxiliary index.
            index = 0
            # Setting the time origin.
            t0 = n*self.ito
            # Iterating over correlation steps.
            for t in range(t0, t0 + self.ncs): 
                # Iterating over particles indexes.
                vT = 0.0
                for i in range(self.npart):
                    vT += self.velocities[t][0][i] * self.velocities[t0][0][i]+\
                         self.velocities[t][1][i] * self.velocities[t0][1][i]+\
                         self.velocities[t][2][i] * self.velocities[t0][2][i] 
                    #v0 = self.velocities[t0][0][i] * self.velocities[t0][0][i]+\
                    #     self.velocities[t0][1][i] * self.velocities[t0][1][i]+\
                    #     self.velocities[t0][2][i] * self.velocities[t0][2][i] 
                    #self.vacf[index] += vi#/v0
                self.vacf[index] += vT/self.npart
                index += 1
        #self.vacf /= (self.npart*self.nto)
        self.vacf /= (self.nto)


    def output(self):
        data_msd = open('data_msd.dat', 'w')
        data_vacf = open('data_vacf.dat', 'w')
        for i in range(len(self.msd)):
            data_msd.write(str(self.time[i]) + ' ' + str(self.msd[i]) + '\n') 
        for i in range(len(self.vacf)):
            data_vacf.write(str(self.time[i]) + ' ' + str(self.vacf[i]) + '\n')
        data_msd.close()
        data_vacf.close()

    def linReg(self):
        x = numpy.array(self.time)
        y = numpy.array(self.msd)
        S = len(x)
        Sx = numpy.sum(x)
        Sy = numpy.sum(y)
        Sxx = numpy.sum(x*x)
        Sxy = numpy.sum(x*y)
        Delta = Sxx*S - Sx**2
        a = (S*Sxy - Sx*Sy)/Delta
        b = (Sy*Sxx - Sx*Sxy)/Delta
        return a, b 

    def integrate(self):
        I = 0.0
        for i in range(len(self.time)-1):
           s = 0.5*(self.vacf[i+1] + self.vacf[i])*(self.time[i+1] - self.time[i])
           I += s
        return I 

    def getTime(self):
        return [i*self.dt*self.fil for i in range(self.ncs)]

    def plotMSD(self):
        pylab.clf()
        pylab.plot(self.time, self.msd, linewidth=2)
        pylab.savefig('msd.png') 

    def plotVACF(self):
        pylab.clf()
        pylab.plot(self.time, self.vacf, linewidth=2)
        pylab.savefig('vacf.png') 

    def writeData(self):
        dataFile = open('data.dat', 'w')
        for sample in self.positions:
            for comp in sample: 
                dataFile.write(" ".join([str(i) for i in comp]) + '\n') 
            dataFile.write('\n')
        dataFile.close()


