#!/usr/bin/env python
""" First steps in a particle-in-cell routine

Author: Bob Wimmer, Inst. f. Exp. & Appl. Physics,
                    Univ. Kiel, Germany

Date: June 14, 2009

Version: 0.0

Description: Follow ES1 code as described in C.K.Birdsall and A.B.Langdon,
"Plasma Physics via Computer Simulation" , Institute of Physics, Series in PlasmPhysics, Taylor & Francis Group, New York, NY, USA, 2005, (ISBN 0-7503-1035-1)

Particle positions and velocities are continuous in phase space, fields are derived at grid points. Time advances are performed as a leapfrog scheme.

1.) Initialize
2.) Temporal evlution:
2a.) integrate equations of motion for particles and derive v and x
2b.) weight grid points to derive charge and current densities
2c.) derive fields (E,B) at grid points from charge and current densities
2d.) derive fields at particle positions so we can derive forces again in step 2a.
3.) end

Conventions:
number particles with i
number grid points with j

classes:
one class for particles
one class for grid

"""

from numpy import *
from math import *
import Gnuplot, Gnuplot.funcutils
import time
from numpy.linalg import norm

#Define global values
ng = 32      #number of grid points (power of two for FFT efficiency)
L = 3.e-3       #physical length to be gridded
nsp = 2      #number of species
np  = 256 #number of particles per species (np > ng)
duration = 1024 #duration of simulation
ampl = 1.e-2*L   #amplitude of sinusoidal perturbation
k = 2.*pi/L     #Wave number of sinusoidal perturbation
q_e_2_m_e = -1.758820150e11 #electron charge-to-mass ratio
q_e_2_m_p = 9.57883392e7    #proton charge-to-mass ratio
m_e = 9.10938215e-31
m_e_2_m_p = 5.4461702177e-4
m_p_2_m_e = 1836.15267247
one_amu = 1.660538782e-27
m_p = 1.672621637e-27
m_p_amu = 1.00727646677
m_4He = 6.64465620e-27
q_e = 1.602176487e-19
c_light = 299792458.
mu_0 = 4.*pi*1.e-7
eps_0 = 8.85418781762e-12
dt = 2.e-0/sqrt(np/L**3*q_e**2/eps_0/m_e)
#dt = 1.

def lorentz_gamma(vel):
    """return the lorentz correction factor"""
    speed =norm(vel)
    if speed > c_light: return "speed has to be less than speed of light!"
    else:
       return 1./sqrt(1.-(speed/c_light)**2)

def summary_plot(E_kin,eE_thermal,pE_thermal,E_field,E_tot,gs):
    """plot speeds and thermal energy history"""
    Ekin = array(E_kin)
    eEt = array(eE_thermal)
    pEt = array(pE_thermal)
    Ef = array(E_field)
    Ett = array(E_tot)
    t= Ekin[:,0]
    #print t
    Ek = Ekin[:,1]
    pt = pEt[:,1]
    et = eEt[:,1]
    ef = Ef[:,1]
    ett = Ett[:,1]
    datEk = Gnuplot.Data(t,Ek)
    datpt = Gnuplot.Data(t,pt)
    datet = Gnuplot.Data(t,et)
    datef = Gnuplot.Data(t,ef)
    datett = Gnuplot.Data(t,ett)
    gs('set xlabel time')
    gs('set ylabel energy')
    maxe = max(pt.max(),et.max())
    mine = min(pt.min(),et.min())
    #print 'mine and maxe: ',  mine, maxe
    #gs('set yrange ['+str(mine)+':'+str(maxe)+']')
    gs.plot(datEk,datef,datpt,datet,datett)

def plot_state(distrib,grid,graph,time):
    """plot the current distribution"""
    xx = []
    vv = []
    xx2 = []
    dd2 = []
    for p in distrib.particles:
        xx.append(p.pos)
        vv.append(p.vel)
        #print 'particle kind, position, and speed is: \t', p.kind,',\t', p.pos, ',\t', p.vel
    graph.reset()
    string = 'time = ' + str(time)
    graph('set title "'+string+'"')
    graph('set ylabel "y-axis [m/s]"')
    graph('set xlabel "x-axis [m]"')
    #g('set term postscript eps color')
    #g('set output "dipole_drift.eps"')
    graph('set style data points')
    #maxf = max(grid.field)
    maxf = fabs(max(vv))
    minf = fabs(min(vv))
    maxf = max(maxf,minf)
    #print 'max plot range = ', maxf
    graph('set xrange [0:'+str(grid.L)+']')
    graph('set yrange [-'+str(maxf)+':'+str(maxf)+']')
    dat = Gnuplot.Data(xx,vv)
    dat2 = Gnuplot.Data(grid.X,grid.field*1.e3)
    graph.plot(dat,dat2)
    #z = time.sleep(0.2)

def save_state(filename,distrib,fields):
    """write current state to file filename"""
    outfile = open(filename,'w')
    #string = 'output data'
    #print string
    #outfile.write(string)
    #print '-------------------------------------------------------------------------------'
    for p in distrib.particles:
        string = str(p.kind) + '\t, pos = ' + str(p.pos) + '\t, vel =' + str(p.vel)
        #print string 
        outfile.write(string + '\n')
    outfile.close()

class particle(object):
     """ particles are defined by mass, charge, velocity, and position (vectors)

     """
     KINDS = ["electron", "proton", "positron"]     
     def __init__(self, kind, charge, vel,pos,counter):
        """mass is rest mass, energy is relativistic"""
	self.kind = kind
        self.number = counter
	if kind == "electron":
	   self.mass = m_e_2_m_p*one_amu
	elif kind == "positron":
	   self.mass = m_e_2_m_p*one_amu
        elif kind == "proton":
	   self.mass = m_p_amu*one_amu
        else:
	   print "kind not yet implemented - defaulting to proton"
	   self.mass = m_p
        
	self.charge = charge*q_e       #particle charge
	self.vel = vel   #particle velicity
	self.pos = pos   #particle position 

     def accel(self,grid,dt):
         """Advance particle velocity one time step dt using quantities on grid"""
         #Determine force on particle
         j = int(self.pos/grid.L*ng)    #index of particles position on grid
         #print 'in outer loop of accel, j = ', j, self.pos, grid.X[j], grid.X[j+1]-self.pos
         if j < ng:
             #print '1 in accel, j = ', j, ng, self.pos, grid.X[j], grid.X[j+1]-self.pos, self.pos - grid.X[j]
             E = (grid.X[j+1]-self.pos)/grid.dX*grid.field[j] + (self.pos - grid.X[j])/grid.dX*grid.field[j+1]
         else:
             print 'in accel, j = ', j, self.pos, grid.X[j], grid.X[j+1]-self.pos, self.pos - grid.X[j]
             E = (grid.X[0]-self.pos)/grid.dX*grid.field[ng] + (self.pos - grid.X[ng])/grid.dX*grid.field[0] 
         F = self.charge * E
         #print '1:', grid.X[j], grid.X[j+1], grid.field[j], grid.dX, E
         acc = F/self.mass/lorentz_gamma(self.vel)
         #print self.kind, self.pos, self.vel, ', accel = ', accel, ', dt = ', dt, accel*dt
         self.vel = self.vel + acc*dt
         #print self.vel, self.pos, dt, self.vel*dt
         #print '2:', self.kind, j, self.pos, F, self.vel

     def move(self,grid,dt):
         """Advance particle one time step dt using quantities on grid"""
         oldpos = self.pos
         self.pos = self.pos + self.vel*dt
         #if self.number == 8: print self.kind,' 8 position = \t', oldpos, self.pos, oldpos - self.pos
         #print self.vel
         while self.pos > grid.L:
             self.pos = self.pos - grid.L
             print 'just moved ' +self.kind +' back into range to pos = ' + str(self.pos) + ' from pos =' + str(self.pos+grid.L) + ' > L'
         while self.pos < 0.:
             self.pos = self.pos + grid.L
             print 'just moved ' + self.kind +' back into range to pos = ' + str(self.pos) + ' from pos =' + str(self.pos-grid.L) + ' < 0'

     def advance(self,grid,dt):
         """Advance particle one time step dt using quantities on grid"""
         self.accel(grid,dt)
         self.move(grid,dt)

class distribution(particle):
    """Generate and proagate an ensemble of particles which are distributed in phase space
       according to a specified distribution.
    """

    def __init__(self,grid,N,kind,charge,pos_dist,vel,vel_wid=0.):
        """distribute particles according to pos_dist and vel"""
        self.particles = []
        self.fill_phase_space(grid,N,kind,charge,pos_dist,vel,vel_wid=0.)

    def fill_phase_space(self,grid,N,kind,charge,pos_dist,vel,vel_wid=0.):
        """distribute particles according to pos_dist and vel"""
        for i in xrange(0,N,1):
            #initialize a particle
            p = particle(kind,charge,0.,0.,i) 
            if pos_dist == 'uniform':
                p.pos = grid.L/(N)*(i + 0.5)
                #print i, p.pos
            elif pos_dist == 'sin':
                x = grid.L/N*(i+0.5)
                p.pos = x + ampl*grid.dX*sin(k*x)
                while p.pos > grid.L:
                    p.pos -= grid.L
                    #print p.pos
                while p.pos < 0:
                    p.pos += grid.L
                    #print p.pos
                #print 'kind and pos: ', p.kind, p.pos,x, ampl*grid.dX*sin(k*x)
            elif pos_dist == 'cos':
                x = grid.L/N*(i+0.5)
                p.pos = x + ampl*grid.dX*cos(k*x)
                while p.pos > grid.L:
                    p.pos -= grid.L
                    #print p.pos
                while p.pos < 0:
                    p.pos += grid.L
                    #print p.pos
                #print 'kind and pos: ', p.kind, p.pos,x, ampl*grid.dX*sin(k*x)
            elif pos_dist == 'random':
                p.pos = random.uniform(0.,grid.L)
            elif pos_dist == 'linear':
                x = grid.L/N*(i+0.5)
                p.pos = x**2/grid.L
            else:
                print 'this position distribution not yet implemented, using uniform'
                p.p = grid.L(N)*i

            if vel_wid > 0.:
                p.vel = vel + random.standard_normal()*vel_wid
            else:
                p.vel = vel
            self.particles.append(p)

    def add(self,grid,N,kind,charge,pos_dist,vel,vel_wid=0.):
        """distribute particles according to pos_dist and vel"""
        self.fill_phase_space(grid,N,kind,charge,pos_dist,vel,vel_wid)

    def __clear__(self):
	self.particles = []

    def __str__(self):
	for particle in self.particles:
	    print particle
        rep= "-----------------"
	return rep

    def add_particle(self,particle):
	"""add a particle to the distribution"""
	self.particles.append(particle)

    def delete_particle(self,particle):
	"""delete a particle from the distribution"""
	self.particles.remove(particle)

    def advance(self,grid,dt):
        """Advance all particles by one time step dt"""
        for p in self.particles:
            #print 'advance: ' + repr(p.pos) + '\t' + repr(p.vel)
            p.advance(grid,dt)

    def advance_vel(self,grid,dt):
        """Advance all particle velocities by one time step dt"""
        for p in self.particles:
            p.accel(grid,dt)

    def print_stat(self):
        """print statiscs"""
        psumpos = 0.
        psumvel = 0.
        pcnt = 0
        esumpos = 0.
        esumvel = 0.
        ecnt = 0
        for p in self.particles:
            if p.kind == 'proton':
                psumpos += p.pos
                psumvel += p.vel
                pcnt += 1
            if p.kind == 'positron':
                psumpos += p.pos
                psumvel += p.vel
                pcnt += 1
            if p.kind == 'electron':
                esumpos += p.pos
                esumvel += p.vel
                ecnt += 1           
        pavepos = psumpos/pcnt
        pavevel = psumvel/pcnt
        eavepos = esumpos/ecnt
        eavevel = esumvel/ecnt
        print 'proton/electron average position/velocity: \t', pavepos, pavevel, eavepos, eavevel

    def print_particle(self,number):
        """print position and velocity of particle number number"""
        print 'pos and vel of particle ', self.particles[number].pos, self.particles[number].vel


    def E_kin(self):
        """return bulk kinetic energy of distribution"""
        E_kin = 0.
        for p in self.particles:
            E_kin += p.mass*p.vel**2
        E_kin *= 0.5
        return E_kin

    def E_therm(self):
        pE_t = 0.
        eE_t = 0.
        pv = 0.
        pcnt = 0
        ev = 0.
        ecnt = 0
        for p in self.particles:
            if p.kind == 'proton':
                pv += p.vel
                pcnt +=1
            elif p.kind == 'positron':
                pv += p.vel
                pcnt +=1
            elif p.kind == 'electron':
                ev += p.vel
                ecnt += 1
        if pcnt > 0: pv /= pcnt
        if ecnt > 0: ev /= ecnt
        for p in self.particles:
            if p.kind == 'proton':
                pE_t += p.mass*(p.vel -pv)**2
            elif p.kind == 'positron':
                pE_t += p.mass*(p.vel -pv)**2
            elif p.kind == 'electron':
                eE_t += p.mass*(p.vel -ev)**2
        pE_t *= 0.5
        eE_t *= 0.5
        return pE_t, eE_t

    def moment(self,n,subp=0., sube =0.):
        """return various moments"""
        pmom = 0.
        pcnt = 0
        emom = 0.
        ecnt = 0
        for p in self.particles:
            if p.kind == 'proton':
                pmom += p.vel**n - subp
                pcnt +=1
            elif p.kind == 'positron':
                pmom += p.vel**n - subp
                pcnt +=1
            elif p.kind == 'electron':
                emom += p.vel**n - sube
                ecnt += 1
        return pmom, pcnt, emom, ecnt

class grid(object):
    """define grid properties and methods"""

    def __init__(self):
        """initialize grid"""
        self.X = linspace(0.,L,ng+1)   #grid (array)
        self.L = L
        self.dX = L/float(ng)
        self.field = zeros(ng+1)       #the field, e.g. electric or magnetic
        self.dens = zeros(ng+1)        #the corresponding density, e.g., charge or current
        self.init = True

    def update(self,distribution):
        """update densities and fields"""
        #update density
        self.dens[0:ng+1] = 0.
        #print 'self.dens = \n', self.dens
        for p in distribution.particles:
            j = int(p.pos/self.L*ng)    #index of particle position on grid, int always gives lower integer
#            self.dens[j] = self.dens[j] + p.charge
            #print p.pos, p.kind, j, 0, ng, self.dens[j]
            if (j < ng-1 and j >= 0):
                #print 'update interior: ', j, p.pos, p.vel, p.kind, self.dens[j], self.dens[j+1]
                self.dens[j] = self.dens[j] + p.charge*(self.X[j+1] - p.pos)/self.dX
                self.dens[j+1] = self.dens[j+1] + p.charge*(p.pos - self.X[j])/self.dX
            elif j == ng-1:
                #print 'update upper edge: ', j, p.pos, p.vel, p.kind, self.dens[j], self.dens[0]
                self.dens[j] = self.dens[j] + p.charge*(self.X[j+1] - p.pos)/self.dX
                self.dens[0] = self.dens[0] + p.charge*(p.pos - self.X[j])/self.dX
            #elif j == 0:
            #    #print 'update lower edge: ', j, p.pos, p.vel, p.kind, self.dens[j], self.dens[0]
            #    self.dens[j] = self.dens[j] + p.charge*(p.pos - self.X[0])/self.dX
            #    self.dens[ng] = self.dens[ng] + p.charge*(p.pos - self.X[ng])/self.dX
            else:
                print '?????? j = ', j, p.pos, p.kind
        self.dens[ng] = self.dens[0]

        #Solve for field using trapezoidal rule and setting self.field[0] = 0.
        self.field[0:ng] = 0.
        j = 0
        for j in xrange(0,ng):
            #print j, self.field[j], self.dens[j], self.dens[j+1], 0.5*(self.dens[j+1] + self.dens[j])*self.dX/eps_0
            self.field[j+1] = self.field[j] + 0.5*(self.dens[j+1] + self.dens[j])*self.dX/eps_0
        #this can be done quicker the following way
        #self.field[1:ng] = self.field[0:ng-1] + 0.5*(self.dens[1:ng] + self.dens[0:ng-1])*self.dX/eps_0
        sum = 0.
        for j in xrange(0,ng):
            sum += self.field[j]
            #print '-----', j, ' ---:', self.dens[j]
        self.field[0:ng+1] -= sum/ng
        if self.init:
            gd = Gnuplot.Gnuplot()
            gf = Gnuplot.Gnuplot()
            densdat = Gnuplot.Data(self.X,self.dens)
            fielddat = Gnuplot.Data(self.X,self.field)
            gd('set ylabel "dens"')
            gd('set title "density"')
            gf('set ylabel field')
            gf('set title "field"')
            gd.plot(densdat)
            gf.plot(fielddat)
            raw_input('press return to continue')
            self.init=False
        #print 'density and field at point 4 and 12:\t', self.dens[4], self.field[4],self.dens[12], self.field[12] 

    def print_stat(self):
        """print statistics"""
        avedens = 0.
        avefield = 0.
        avecnt = 0
        for j in xrange(0,ng):
            avedens += self.dens[j]
            avefield += self.field[j]
            avecnt +=1
        print 'average density and field:\t', avedens/avecnt, avefield/avecnt

    def energy(self):
        """return energy density stored in field"""
        en = 0.
        for j in xrange(0,ng):
            en += self.field[j]**2
        return 0.5*eps_0*en

#----------------------------------------------------------------------------------------------------------------------#

#this allows us to use this file as a module, but also to call it as a standalone script
if __name__ == "__main__":

    speeds = []
    eE_therm = []
    pE_therm = []
    E_field = []
    E_kin = []
    E_tot = []
    #initialize and save initial conditions to file
    #initialize the field grid
    field = grid()
    #now we have the field in the initial state
    f = distribution(field,np,'proton',1000.,'uniform',0.,0.)    #populate with protons
    #f = distribution(field,np,'electron',-1000.,'uniform',1.,0.)    #populate with electrons
    f.add(field,1.*np,'electron',-1000.,'sin',0.0,0.e-0)  #add electrons
    #f.add(field,0.5*np,'electron',-1000.,'sin',-20.0,0.5e-1)  #add electrons
    #f = distribution(field,np,'positron',1000.,'uniform',1.,0.)    #populate with electrons
    #f.add(field,np,'electron',-1000.,'uniform',1.0,0.0)  #add some streaming electrons
    #f.add(field,np,'electron',-1000.,'uniform',-1.0,0.0)  #add some streaming electrons
    g = Gnuplot.Gnuplot()
    #plot_state(f,field,g,0.)
    #raw_input('Please press the return key to continue...\n')
    print '1111'
    field.update(f)
    print '222222'
    #plot_state(f,field,g,0.)
    #raw_input('Please press the return key to continue...\n')    

    #next need to move velocities (and only velocities) backwards half a time step
    f.advance_vel(field,-0.5*dt)

    filename = 'initial_state.dat'
    save_state(filename,f,field)
    plot_state(f,field,g,0.)
    raw_input('Please press the return key to continue...\n')    

    #some statistics stuff
    Ef = field.energy()
    Ek = f.E_kin()
    pEt, eEt = f.E_therm()
    Et = Ef + Ek + pEt + eEt
    E_field.append([0,Ef])
    E_kin.append([0,Ek])
    pE_therm.append([0,pEt])
    eE_therm.append([0,eEt ])
    E_tot.append([0,Ef+Ek])

    #loop over time steps
    step = 0
    output_steps = [0,4,8]
    t = dt
    t_end = duration*t
    while (t < t_end):
        #advance velocities and then positions
        #print dt
        f.advance(field,dt)

        #update fields
        field.update(f)
        plot_state(f,field,g,t)
        
        #some more statistics stuff
        Ef = field.energy()
        Ek = f.E_kin()
        pEt, eEt = f.E_therm()
        Et = Ef + Ek + pEt + eEt
        E_field.append([t,Ef])
        E_kin.append([t,Ek])
        pE_therm.append([t,pEt])
        eE_therm.append([t,eEt ])
        E_tot.append([t,Ef+Ek])
        
        #once in a while save current state to file
        step += 1
        t = step*dt
        if step in output_steps:
            filename = 'state_at_t_=_' + repr(t)+'.dat'
            save_state(filename,f,field)
        #f.print_stat()
        #field.print_stat()
        #f.print_particle(6)

    #save final state to file
    filename = 'final_state.dat'
    save_state(filename,f,field)
    plot_state(f,field,g,t)
    print '-------------'
    gs = Gnuplot.Gnuplot()
    summary_plot(E_kin,pE_therm,eE_therm,E_field,E_tot,gs)
    raw_input('Please press the return key to continue...\n')
