import pygmsh
import gmsh
import matplotlib.pyplot as plt
import numpy as np
import capytaine as cpt

import meshmagick.hydrostatics as hs
import meshmagick.mesh as mm

from scipy.linalg import block_diag

import logging
logging.basicConfig(level=logging.INFO, format='%(levelname)-8s: %(message)s')

#hexagonplatemeterv3.stl --> right hexagon version
#Simple_Model_Scaledby0.001 v5.stl' --> right peral version
body = cpt.FloatingBody.from_file('Simple_Model_Scaledby0.001 v5.stl') #file name
#body.translate([0,0,-.2032])
body.keep_immersed_part() #clips top of body
#body.show() #uncomment to visualize

#hsd is hydro coeffs, I is inertia, M is mass matrix
hsd = hs.Hydrostatics(mm.Mesh(body.mesh.vertices, body.mesh.faces)).hs_data 

m = hsd['disp_mass']
I = np.array([[hsd['Ixx'], -1*hsd['Ixy'], -1*hsd['Ixz']],
              [-1*hsd['Ixy'], hsd['Iyy'], -1*hsd['Iyz']],
              [-1*hsd['Ixz'], -1*hsd['Iyz'], hsd['Izz']]])
M = block_diag(m, m, m, I)
# body.mass = body.add_dofs_labels_to_matrix(M)

kHS = block_diag(0,0,hsd['stiffness_matrix'],0)
# body.hydrostatic_stiffness = body.add_dofs_labels_to_matrix(kHS)

##[x,x] 0=Surge 1=Sway 2=Heave 3=Roll 4=Pitch 5=Yaw
omega0_estimate = np.sqrt(kHS[2,2]/M[2,2]) #omega_natural = sqrt(k/m)
print(omega0_estimate)