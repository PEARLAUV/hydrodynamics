# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 20:37:22 2021

@author: Sharmi
"""
from capytaine import *
from numpy import infty
from capytaine.ui.vtk.animation import Animation
import logging

logging.basicConfig(level=logging.INFO)
cylinder= VerticalCylinder(length=10.0,radius= 1.0,center=(0,0,0),nx=20, ntheta=20,nr=10,clever=True,name='cyl')
# cylinder.show()

# cylinder.add_all_rigid_body_dofs()
cylinder.add_rotation_dof(name="Pitch")
# cylinder.add_translation_dof(name="Heave")
print(cylinder.dofs.keys())
cylinder.keep_immersed_part()

problem= RadiationProblem(body=cylinder, radiating_dof="Pitch",omega=1,rho=5)

solver=BEMSolver()
result=solver.solve(problem)

print(result.added_masses)
print(result.radiation_dampings)

# animation=Animation(loop_duration=result.period)
# animation.add_body(cylinder,faces_motion=cylinder.dofs["Heave"])
# animation.run()
# direct_linear_solver=cpt.BasicMatrixEngine(linear_solver)