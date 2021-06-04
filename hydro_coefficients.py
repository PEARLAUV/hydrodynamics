import logging
import capytaine as cpt

import numpy as np
np.set_printoptions(precision=3)
np.set_printoptions(linewidth=160)

import meshmagick.hydrostatics as hs
import meshmagick.mesh as mm

from scipy.linalg import block_diag

body = cpt.FloatingBody.from_file("hexagonplate_small_diameterheight.stl", file_format="stl")

#body.translate([0,0,-.2032])
body.show()
body.add_all_rigid_body_dofs()
body.keep_immersed_part()



# --------------------------- added mass and hydrostatic coeff -------------------------
# Inertial properties for neutrally buoyant constant density body
    #you can use meshmagick to find inertia matrix but they use constant density
    #alternatively use a CAD software to acquire inertia matrix
    

hsd = hs.Hydrostatics(mm.Mesh(body.mesh.vertices, body.mesh.faces)).hs_data
m = hsd['disp_mass']
I = np.array([[hsd['Ixx'], -1*hsd['Ixy'], -1*hsd['Ixz']],
              [-1*hsd['Ixy'], hsd['Iyy'], -1*hsd['Iyz']],
              [-1*hsd['Ixz'], -1*hsd['Iyz'], hsd['Izz']]])
M = block_diag(m, m, m, I)
body.mass = body.add_dofs_labels_to_matrix(M)
print(body.mass) #mass matrix


# Hydrostatics
kHS = block_diag(0,0,hsd['stiffness_matrix'],0)
body.hydrostatic_stiffness = body.add_dofs_labels_to_matrix(kHS)

print(body.hydrostatic_stiffness) #hydrostatic coefficient




# ---------------------- find damping coeff + added mass @ frequency -----------------------------
omega=5.0
problem = cpt.RadiationProblem(body=body, radiating_dof="Heave", omega=omega)
solver = cpt.BEMSolver()
result = solver.solve(problem)
print(result.radiation_dampings) #damping coefficient
print(result.added_masses) #added mass coefficient
