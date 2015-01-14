emr3d
=====

Preamble

    Project based on First steps in a particle-in-cell routine. Follow ES1
    code from "Plasma Physics via Computer Simulation" [C.K.Birdsall and A.B.Langdon]

    Author: Bob Wimmer, Inst. f. Exp. & Appl. Physics, Univ. Kiel, Germany.
    Date: July 6, 2009. Version: 0.1

Introduction

    Particle In Cell Method - ElectroMagnetic & Relativistic Code In Three Dimensions  (EMR3D)

Description

    1. Initialize
	    Module "initialization" sets global variants and calculate gamma factor
	 
    2. Evolution
	    Modules "particle" and "cell"

    3. Solution
	    Calculate forces and fields
		
    4. Processing
	    Module "processing" store results in file and plot results
