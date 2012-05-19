import numpy as np
import pylab as plt


class Analyse:
    def __init__(self, simulation):
        self.sampler = simulation.sampler
        self.system  = simulation.system

    def figures( self ):
        print "calculating results"
        # create list of snapshots (10 steps pr. snapshot)
        slist = [50*i for i in range(len(self.sampler.getData()))]
        # create a list of snapshots to represent the averaged results
        slist_avg = [500*i for i in range(len(self.sampler.getData())/10)]
  
        # create energy lists
        epot = [sample[3] for sample in self.sampler.getData()]
        ekin = [sample[4] for sample in self.sampler.getData()]
        etot = [epot[i]+ekin[i] for i in range(len(self.sampler.getData()))]

        #print 'epot', epot, 'len(epot)', len(epot)
        #print 'ekin', ekin, 'len(ekin)', len(ekin)

        # make averages over each 10 elements
        #epot_avg = [np.average(epot[i:i+10]) for i in range(len(epot)/10.0)]
        #ekin_avg = [np.average(ekin[i:i+10]) for i in range(len(ekin)/10.0)]
        epot_avg = [np.average(epot[i:i+10]) for i in range(len(epot)/10)]
        ekin_avg = [np.average(ekin[i:i+10]) for i in range(len(ekin)/10)]
        # since the total energy is constant, we dont worry about that
  
        # plot the energy
        plt.clf()
        plt.suptitle('Potential and Kinetic Energy')
        plt.xlabel('N simulation steps')
        plt.ylabel('Energy [a.u]')
        plt.plot( slist, epot, color='r' )
        plt.plot( slist, ekin, color='b' )
        plt.plot( slist, etot, color='g' )
        plt.plot( slist_avg, epot_avg, color='black' )
        plt.plot( slist_avg, ekin_avg, color='black' )
        plt.savefig('energy.png')
  
        # create temperature and pressure lists
        t = [sample[5] for sample in self.sampler.getData()]
        p = [sample[6] for sample in self.sampler.getData()]
        # create averages
        t_avg = [np.average(t[i:i+10]) for i in range(len(t) / 10)]
        p_avg = [np.average(p[i:i+10]) for i in range(len(p) / 10)]
  
        # plot temperature
        plt.clf()
        plt.suptitle('Temperature')
        plt.xlabel('N simulation steps')
        plt.ylabel('Temperature [a.u]')
        plt.plot( slist, t, color='b')
        plt.plot( slist_avg, t_avg, color='black')
        plt.savefig('temperature.png')
  
        # plot pressure
        plt.clf()
        plt.suptitle('Pressure')
        plt.xlabel('N simulation steps')
        plt.ylabel('Pressure [a.u]')
        plt.plot( slist, p, color='b')
        plt.plot( slist_avg, p_avg, color='black')
        plt.savefig('pressure.png')
        print "Figures done."
  
    def rdf(self):
        print "creating rdf"
        nr = 75
        g = np.zeros( nr, float )
        dr = self.system.getBox()[0] / (2.0 *nr)
        cnt = 0
        for snapshot in self.sampler.getData():
            cnt += 1
            print "working [%i/%i] %6.2f %%" % (cnt, len(self.sampler.getData()), float(cnt) / len(self.sampler.getData()) * 100)
            X = snapshot[0]
            Y = snapshot[1]
            Z = snapshot[2]
            for i in range( self.system.getN() ):
                for j in range( self.system.getN() ):
                    if( i > j ):
                        x = X[i]-X[j]
                        y = Y[i]-Y[j]
                        z = Z[i]-Z[j]
                        x = x - self.system.getBox()[0] * round( x / self.system.getBox()[0] )
                        y = y - self.system.getBox()[1] * round( y / self.system.getBox()[1] )
                        z = z - self.system.getBox()[2] * round( z / self.system.getBox()[2] )
                        rij = np.sqrt( x**2 + y**2 + z**2 )
                        if( rij < min(self.system.getBox()) / 2.0 ):
                            k = int( rij / dr )
                            g[k] += 2
    
        print "normalizing rdf"
        r = []
        for i in range( 0, nr ):
            r.append( dr*(i+0.5) )
            vb = (dr**3) * ((i+1)**3 - (i)**3)
            nid = 4.0 * np.pi*vb*self.system.getRho()
            g[i] = g[i] / (nid * self.system.getN() * len(self.sampler.getData()))
  
        plt.clf()
        plt.plot( r, g, label = "g(r), rho = %6.2f, T = %6.2f" % (self.system.getRho(), self.system.getT()) )
        plt.legend()
        plt.savefig( "rdf.png" )

        print "RDF done."

    def showData(self):
        print "Displaying the sampled data."
        for i in self.sampler.getData():
            print "*"*50
            print "X", i[0]
            """
            print "Y", i[1]
            print "Z", i[2]
            print "ePot:       ", i[3]
            print "eKin:       ", i[4]
            print "Temperature:", i[5]
            print "Pressure:   ", i[6]
            print "Velocities: ", i[7]
            """
            print

    def getVelocityDistribution(self):
        # Produce a single MB-distribution figure.
        velocities       = []
        for i in self.sampler.getData():
            velocities.append(i[7])
        npart = self.system.getN()
        velMagnitudes = []

        for sample in range(len(velocities)):
            # Skipping the component iteration because there
            # are always three components.  
            #vlist = []
            for part in range(npart): 
                # Absolute velocity of particle i for sample sample.
                absvi = np.sqrt( velocities[sample][0][part]*velocities[sample][0][part]
                        + velocities[sample][1][part]*velocities[sample][1][part]
                        + velocities[sample][2][part]*velocities[sample][2][part] )
                #vlist.append(absvi) 
                velMagnitudes.append(absvi)
            #velMagnitudes.append(vlist)

        velMagnitudes = np.ravel( np.array(velMagnitudes) )
        #plt.hist(velMagnitudes, bins=40) 
        #plt.savefig('mb-%02.1f.png' % temp)
        return velMagnitudes


