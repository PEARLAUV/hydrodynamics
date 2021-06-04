# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 12:26:26 2021

@author: Sharmi
"""

#Plot heave added mass
import logging
import numpy as np
from capytaine import *
logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s:\t%(message)s")
body = FloatingBody.from_file('PEARLmesh_sharmi v1v5.stl', file_format='stl') #filename

body.keep_immersed_part() #clip body
body.show()


body.add_translation_dof(name="Heave") #add dof to body

# Define the range of frequencies as a Numpy array
omega_range = np.linspace(0.0, 1,10)
rho=3875
#define problems
problems = [
    RadiationProblem(body=body, radiating_dof=dof, omega=omega, rho=rho)
    for dof in body.dofs
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
direct_linear_solver = BasicMatrixEngine(linear_solver='direct')
solver = BEMSolver(engine=direct_linear_solver)
results = [solver.solve(pb) for pb in sorted(problems)]
# The 'sorted' function ensures that the problems are sequentially
# treated in an optimal order.

# Gather the computed added mass into a labeled array.
data = assemble_dataset(results)


###################
## Plot ##
###################
dof = "Heave"
#print(data['added_mass'])
import matplotlib.pyplot as plt
for dof in body.dofs:
    plt.plot(
        omega_range,
        data['added_mass'].sel(radiating_dof=dof, influenced_dof=dof),
        label=dof,
        marker='o',
    )
plt.xlabel('omega')
plt.ylabel('added mass')
plt.legend()
plt.tight_layout()
plt.show()