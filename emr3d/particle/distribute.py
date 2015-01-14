#!/usr/bin/env python
# -*- coding: utf-8 -*-

from base.pic_2 import *

#from initialization.constant import *
#from grid.general import *
from particle.general import *

from numpy import *
from math import *
import Gnuplot, Gnuplot.funcutils
import time
from numpy.linalg import norm

class Distribution(Particle):
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
            p = Particle(kind,charge,0.,0.,i) 
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
            E_kin += p.mass*p.vel*p.vel_old
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
