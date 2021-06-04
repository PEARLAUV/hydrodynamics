# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 22:02:18 2021

@author: sharmis
"""

import logging
import capytaine as cpt

import numpy as np
np.set_printoptions(precision=3)
np.set_printoptions(linewidth=160)

import meshmagick.hydrostatics as hs
import meshmagick.mesh as mm

from scipy.linalg import block_diag
import matplotlib.pyplot as plt


bem_solver = cpt.BEMSolver()
body = cpt.FloatingBody.from_file("PEARLmesh_sharmi v1v4.stl", file_format="stl")

body.add_all_rigid_body_dofs()
body.keep_immersed_part()

hsd = hs.Hydrostatics(mm.Mesh(body.mesh.vertices, body.mesh.faces)).hs_data

# Inertial properties for neutrally buoyant constant density body
    #you can use meshmagick to find inertia matrix but they use constant density
    #alternatively use a CAD software to acqurie inertia matrix
m = 130.78
I = np.array([[34.298, 1.133E-10,0],
                      [1.133E-10,34.303,0],
                      [0,0,33.331]])
M = block_diag(m, m, m, I)
body.mass = body.add_dofs_labels_to_matrix(M)
# print(body.mass)

# Hydrostatics
kHS = block_diag(0,0,hsd['stiffness_matrix'],0)
body.hydrostatic_stiffness = body.add_dofs_labels_to_matrix(kHS)

# change omega values
omega=np.linspace(0.1, 2,10)
wave_direction=0.0
#body.keep_only_dofs(['Heave'])


problems = [cpt.RadiationProblem(omega=w, body=body, radiating_dof=dof) for dof in body.dofs for w in omega]
problems += [cpt.DiffractionProblem(omega=w, body=body, wave_direction=wave_direction) for w in omega]
results = [bem_solver.solve(problem) for problem in problems]
*radiation_results, diffraction_result = results
dataset = cpt.assemble_dataset(results)

    # COMPUTE RAO
dataset['RAO'] = cpt.post_pro.rao(dataset, wave_direction=wave_direction)
# print(dataset['RAO'])
# print(dataset['RAO'].data)
#plt.plot(omega,dataset['RAO'].data)

#change dof to heave, pitch, surge
dof="Heave"
plt.plot(omega,np.abs(dataset['RAO'].sel(radiating_dof=dof))) 
plt.xlabel('omega')
plt.ylabel('RAO ' + dof)
plt.show()
