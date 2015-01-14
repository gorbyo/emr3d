#!/usr/bin/env python
# -*- coding: utf-8 -*-

from base.pic_2 import *

#from initialization.constant import *
#from particle.general import *

class Grid(object):
    """define grid properties and methods"""


    def __init__(self):
        """initialize grid"""
        self.X = linspace(0.,L,ng+1)   #grid (array)
        self.L = L
        self.dX = L/float(ng)
        self.field = zeros(ng+1)       #the field, e.g. electric or magnetic
        self.dens = zeros(ng+1)        #the corresponding density, e.g., charge or current
        self.field_ft = zeros(ng+1)       #the field, e.g. electric or magnetic, but Fourier transformed
        self.dens_ft = zeros(ng+1)        #the corresponding density, e.g., charge or current, but Fourier transformed
        self.freq = zeros(ng+1)
        self.pot  = zeros(ng+1)        #the potential
        self.pot_ft  = zeros(ng+1)        #Fourier transform of the potential
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
        #Once all is known compute Fourier transforms of the quantities
        self.dens_ft = fft(self.dens)
        self.field_ft = fft(self.field)
        step = L/ng
        self.freq = L*fftfreq(ng, d=step)
        self.freq[0] = 0.01
        self.pot_ft[0:ng] = self.dens_ft[0:ng]/self.freq[0:ng]**2



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
            en += self.pot_ft[j]*self.dens_ft[j]
            #en += self.field[j]**2
        #return 0.5*eps_0*en
        return 1./L*eps_0*en