#!/usr/bin/env python
# -*- coding: utf-8 -*-

from base.pic_2 import *


# this allows us to use this file as a module, but also to call it as a standalone script
if __name__ == "__main__":

    speeds = []
    eE_therm = []
    pE_therm = []
    E_field = []
    E_kin = []
    E_tot = []
    # initialize and save initial conditions to file
    # initialize the field grid
    field = grid()
    # now we have the field in the initial state
    f = distribution(field,np,'proton',1000.,'uniform',0.,0.)    # populate with protons
    # f = distribution(field,np,'electron',-1000.,'uniform',1.,0.)    # populate with electrons
    f.add(field,1*np,'electron',-1000.,'sin',0.0,0.e-0)  #add electrons
    # f.add(field,0.5*np,'electron',-1000.,'sin',-20.0,0.5e-1)  #add electrons
    # f = distribution(field,np,'positron',1000.,'uniform',1.,0.)    #populate with electrons
    # f.add(field,np,'electron',-1000.,'uniform',1.0,0.0)  #add some streaming electrons
    # f.add(field,np,'electron',-1000.,'uniform',-1.0,0.0)  #add some streaming electrons
    g = Gnuplot.Gnuplot()
    # plot_state(f,field,g,0.)
    # raw_input('Please press the return key to continue...\n')
    # print '1111'
    field.update(f)
    # print '222222'
    # #plot_state(f,field,g,0.)
    # raw_input('Please press the return key to continue...\n')

    # next need to move velocities (and only velocities) backwards half a time step
    f.advance_vel(field,-0.5*dt)

    filename = 'initial_state.dat'
    save_state(filename,f,field)
    plot_state(f,field,g,0.)
    raw_input('Please press the return key to continue...\n')    

    # some statistics stuff
    Ef = field.energy()
    Ek = f.E_kin()
    pEt, eEt = f.E_therm()
    Et = Ef + Ek + pEt + eEt
    E_field.append([0,Ef])
    E_kin.append([0,Ek])
    pE_therm.append([0,pEt])
    eE_therm.append([0,eEt ])
    E_tot.append([0,Ef+Ek])

    # loop over time steps
    step = 0
    output_steps = [0,4,8]
    t = dt
    t_end = duration*t
    while (t < t_end):
        # advance velocities and then positions
        # print dt
        f.advance(field,dt)

        # update fields
        field.update(f)
        plot_state(f,field,g,t)
        
        # some more statistics stuff
        Ef = field.energy()
        Ek = f.E_kin()
        pEt, eEt = f.E_therm()
        Et = Ef + Ek + pEt + eEt
        E_field.append([t,Ef])
        E_kin.append([t,Ek])
        pE_therm.append([t,pEt])
        eE_therm.append([t,eEt ])
        E_tot.append([t,Ef+Ek])
        
        # once in a while save current state to file
        step += 1
        t = step*dt
        if step in output_steps:
            filename = 'state_at_t_=_' + repr(t)+'.dat'
            save_state(filename,f,field)
        # f.print_stat()
        # field.print_stat()
        # f.print_particle(6)


    # save final state to file
    filename = 'final_state.dat'
    save_state(filename,f,field)
    plot_state(f,field,g,t)
    print '-------------'
    gs = Gnuplot.Gnuplot()
    summary_plot(E_kin,pE_therm,eE_therm,E_field,E_tot,gs)
    raw_input('Please press the return key to continue...\n')