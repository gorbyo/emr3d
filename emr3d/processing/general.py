#!/usr/bin/env python
# -*- coding: utf-8 -*-

import Gnuplot, Gnuplot.funcutils
from numpy import *

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
    gs.plot(datEk,datef,datett)#,datpt,datet)#,datett)
    

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
