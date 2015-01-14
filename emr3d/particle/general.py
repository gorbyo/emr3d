#!/usr/bin/env python
# -*- coding: utf-8 -*-

from initialization.constant import *

class Particle(object):
    """ particles are defined by mass, charge, velocity, and position (vectors)"""
    KINDS = ["electron", "proton", "positron"]
    def __init__(self, kind, charge, vel,pos,counter):
        """mass is rest mass, energy is relativistic"""
        self.kind = kind
        self.number = counter
        self.charge = charge*q_e       #particle charge
        self.vel = vel   #particle velocity
        self.vel_old = 0. #old particle velocity
        self.pos = pos   #particle position
        
        if kind == "electron":
            self.mass = m_e_2_m_p*one_amu
        elif kind == "positron":
            self.mass = m_e_2_m_p*one_amu
        elif kind == "proton":
            self.mass = m_p_amu*one_amu
        else:
            print "kind not yet implemented - defaulting to proton"
            self.mass = m_p
        
     

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
        acc = F/self.mass/gamma_lorentz_factor(self.vel)
        #print self.kind, self.pos, self.vel, ', accel = ', accel, ', dt = ', dt, accel*dt
        self.vel_old = self.vel
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

