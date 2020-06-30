#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 14:21:55 2020

@author: maha
"""

# Initialize
import logging
import numpy as np
import capytaine as cpt

# Set up logging
logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s:\t%(message)s")

###########################
## PEARL geometry inputs ##
###########################
# These are design variables, all can be changed

D_f = 5 # effective diameter of top float
t_f = 2 # thickness of top float
alpha = 0.5 # % of top float that is submerged

N_s = 1 # number of column supports
D_s = 1 # diameter of column support(s)
t_s = 6 # length of column support
r_c = 2 # radius of circle that centerpoints of columns trace out if N_s > 1

D_d = 5 # effective diameter of damping plate
t_d = 0.5 # thickness of damping plate

tol = 1e-6 # tolerance to ensure top float submerged
mult = 2 # multiplier to distance columns away from top float to ensure algorithm convergence

rho_structure = 650 # kg/m^3 average density of platform (voids and material)

######################
## Mesh generation ##
#####################
# This automatically generates the mesh from the design variables chosen above

mass = rho_structure*np.pi*(t_d*(D_d/2)**2 + t_s*(D_s/2)**2 + t_f*(D_f/2)**2)
rho_sw = 1023 #kg/m^3
## top float ##
r_f = D_f/2 # effective radius of top float
s_f = r_f/np.sqrt(3) # side length of top float hexagon
h_f = alpha*t_f # submerged thickness of top float


# Initialize floating body by generating a geometric mesh
top_float = cpt.VerticalCylinder(
    length=h_f, radius=r_f,  # Dimensions
    center=(0, 0,-h_f/2-tol),        # Position
    nr=25, nx=20, ntheta=30,   # Fineness of the mesh
    name = "top_float")
#top_float.add_all_rigid_body_dofs()
#top_float.add_translation_dof(name="Heave")
# top_float.show()

## support column(s) ##
r_s = D_s/2; # effective radius of support column(s)

if N_s == 1:
    support_column = cpt.VerticalCylinder(
        length=t_s, radius=r_s,  # Dimensions
        center=(0, 0, -(h_f+t_s/2)-mult*tol),        # Position
        nr=10, nx=50, ntheta=30,   # Fineness of the mesh
        name = "support_column")
    #support_column.add_all_rigid_body_dofs()
    #support_column.add_translation_dof(name="Heave")

else:
    support_coldict = {}
    for ii in range(1,N_s+1):
        support_coldict[ii] = cpt.VerticalCylinder(
            length=t_s, radius=r_s,  # Dimensions
            center=(r_c*np.cos(2*np.pi*(ii-1)/N_s), r_c*np.sin(2*np.pi*(ii-1)/N_s), -(h_f+t_s/2)-mult*tol),        # Position
            nr=10, nx=50, ntheta=30,   # Fineness of the mesh
            name = "support_column" + str(ii))
        #support_coldict[ii].add_all_rigid_body_dofs()
        #support_coldict[ii].add_translation_dof(name="Heave")
    support_column = support_coldict[1]
    for ii in range(2,N_s+1):
        support_column = support_column + support_coldict[ii]
#support_column.show()

## damping plate ##
r_d = D_d/2; # effective radius of the damping plate

damping_plate = cpt.VerticalCylinder(
    length=t_d, radius=r_d,  # Dimensions
    center=(0, 0, -(h_f+t_s+t_d/2)-mult*tol),        # Position
    nr=25, nx=10, ntheta=30,   # Fineness of the mesh
    name = "damping_plate")
#damping_plate.add_all_rigid_body_dofs()
#damping_plate.add_translation_dof(name="Heave")
#damping_plate.show()

## Combine into one body ##
PEARL_body = top_float + support_column + damping_plate;
#PEARL_body = top_float + support_column;
PEARL_body.show()
PEARL_body.keep_immersed_part()
#PEARL_body.show()

#########################
## Problem formulation ##
#########################
# This generates the problem to be solved using BEMSolver

# Automatically add the six degrees of freedom of a rigid body
#PEARL_body.add_all_rigid_body_dofs()
PEARL_body.add_translation_dof(name="Heave")

# Define the range of frequencies as a Numpy array
omega_range = np.linspace(0.1, 6.0, 20)
#omega_range = np.linspace(0.1, 1.0, 10) # [rad/s] (?)

# Set up the problems: we will solve a radiation problem for each
# degree of freedom of the body and for each frequency in the
# frequency range.

problems = [
    cpt.RadiationProblem(body=PEARL_body, radiating_dof=dof, omega=omega, rho=rho_sw)
    for dof in PEARL_body.dofs
    for omega in omega_range
]

# Water density, gravity and water depth have not been specified.
# Default values are used.

###################
## Solve problem ##
###################
# This uses BEMSolver to solve the problem(s) generated above

# Solve all radiation problems
#solver = cpt.BEMSolver()
direct_linear_solver = cpt.BasicMatrixEngine(linear_solver='direct')
solver = cpt.BEMSolver(engine=direct_linear_solver)
results = [solver.solve(pb) for pb in sorted(problems)]
# The 'sorted' function ensures that the problems are sequentially
# treated in an optimal order.

# Gather the computed added mass into a labelled array.
data = cpt.assemble_dataset(results)

##################
## Plot results ##
##################
# This plots the results of the added mass computation

dof = "Heave"
A33 = data['added_mass'].sel(radiating_dof=dof, influenced_dof=dof)
C33 = np.pi*data.rho*data.g/4*D_f**2
omega_0 = np.sqrt(C33/(A33+mass))

# Plot the added mass of each dofs as a function of the frequency
import matplotlib.pyplot as plt
plt.figure()
#for dof in PEARL_body.dofs:
#    plt.plot(
#        omega_range,
#        data['added_mass'].sel(radiating_dof=dof, influenced_dof=dof),
#        label=dof,
#        marker='o',
#    )

plt.plot(
        omega_range,
        data['added_mass'].sel(radiating_dof=dof, influenced_dof=dof),
        label=dof,
        marker='o',
    )

plt.xlabel('omega [rad/s]')
plt.ylabel('added mass [kg]')
plt.legend()
plt.tight_layout()
plt.show()

plt.plot(
        omega_range,
        data['radiation_damping'].sel(radiating_dof=dof, influenced_dof=dof),
        label=dof,
        marker='o',
    )

plt.xlabel('omega [rad/s]')
plt.ylabel('radiation damping [Ns/m]')
plt.legend()
plt.tight_layout()
plt.show()

## Still need to:
# Add in mass inertia matrix computation
# Add in stiffness matrix computation
# Solve for RAO and plot

