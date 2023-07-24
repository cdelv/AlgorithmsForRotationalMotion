#---------------------------------------------------------------------
# Chess pieces bouncing inside a box
# 2023 © Vasileios Angelidakis <vasileios.angelidakis@fau.de>
# 2023 © Carlos Andrés del Valle <cdelv@unal.edu.co>
# Institute for Multiscale Simulation
#---------------------------------------------------------------------
from yade import ymport, plot
import random
random.seed(12)
import numpy as np
import csv
import os

import numpy as np
from stl import mesh

#------------------------------------------------------------------------------
t_max = 1.5          # s
DT = 0.3*4.68812776e-06  # s
samples = 250
Vel = 2.0            # m/s
paraview = False

pieces = {}
scaling_factors = {}
frame = 0

# Integration Method
integrationType='Fincham1992'
integrationType='Omelyan1998'
integrationType='delValle2023'

#------------------------------------------------------------------------------
# New units for the simulation (to avoid numerical problems)
E_new = 1 #1e9    # new unit of pressure
r_new = 1 #1e-2   # new unit of lenght
s_new = 1 #1.9311703796025776e6 # new unit of time (Rayleigh timestep)^-1

#------------------------------------------------------------------------------
Energy_Err = []
def Energy_Error():
    Ef = utils.kineticEnergy() + O.engines[2].lawDispatcher.functors[0].normElastEnergy() 
    Energy_Err.append(100*abs(E0 - Ef)/E0)
    return Energy_Err[-1]/len(Energy_Err)

# Function to load the clumps
def read_x_y_z_r_intoList(fileName,shift,scale):
    if not shift: shift = Vector3.Zero
    infile = open(fileName,"r")
    lines = infile.readlines()
    infile.close()
    c1_list=[]
    for line in lines:
        data = line.split(",")
        if (data[0][0] == "#"):
            continue
        else:
            pos = Vector3(float(data[0]),float(data[1]),float(data[2]))
            c1_list.append( (shift+scale*pos,scale*float(data[3])) )
    return c1_list

def add_pieces(stl_file_path):
    for b in O.bodies:
        if isinstance(b.shape, Clump):
            if b.id not in pieces:
                pieces[b.id] = stl_file_path

def saveSTLs():
    n_mesh = 1
    global frame
    frame += 1
    for b_id in pieces:
        original_stl = mesh.Mesh.from_file(pieces[b_id])
        original_stl.rotate_using_matrix(scaling_factors[pieces[b_id]]*np.identity(3))
        center_of_mass = np.average(original_stl[:,:3], axis=0)
        original_stl.translate(-center_of_mass)
        original_stl.rotate_using_matrix(np.array(O.bodies[b_id].state.ori.conjugate().toRotationMatrix()))
        q = Quaternion()
        q.setFromTwoVectors((1,0,0),(0,0,-1))
        original_stl.rotate_using_matrix(np.array(q.toRotationMatrix()))
        original_stl.translate(np.array(O.bodies[b_id].state.pos))
        
        transformed_filename = f'paraview/transformed_{n_mesh}_{frame}.stl'
        original_stl.save(transformed_filename)
        n_mesh += 1

#------------------------------------------------------------------------------
# Material for the simulation

# Chess Championship pieces are made of Baswood -> https://www.youtube.com/watch?v=-Tg9xiJ6D6k
# young 10e9 -> https://www.engineeringtoolbox.com/timber-mechanical-properties-d_1789.html
# poisson 0.03 -> https://www.matweb.com/search/datasheet.aspx?matguid=1775e2148e704fd0b10d2c8af6c3efb7&n=1&ckck=1
# density 340 -> https://www.engineeringtoolbox.com/timber-mechanical-properties-d_1789.html
# In the video the wood seem pritty dry.

# However, this makes the time step quite small. Then, we have fancy chess pieces made of marble
# young, density -> https://www.matweb.com/search/datasheet.aspx?matguid=a1c4c37b55e24765bcb63b4665b44f05
# poisson -> https://www.engineeringtoolbox.com/poissons-ratio-d_1224.html
matID = O.materials.append(FrictMat(young=60e9/E_new, poisson=0.25, density=3200*(r_new**2/(E_new*s_new**2)), frictionAngle=0.0, label='spheresMat'))

#------------------------------------------------------------------------------
# Simulation Box

# Facet box
corner1 = Vector3(0.0,0.0,0.0)
corner2 = (0.25/r_new)*Vector3(1.0,1.0,1.0) # 25 x 25 x 25 cm box
id1 = O.bodies.append(geom.facetBox(0.5*(corner1 + corner2), corner2, material=O.materials['spheresMat']))

#------------------------------------------------------------------------------
# Create simulation bodies

corner1 = -0.5*corner2
corner2 = 1.5*corner2
unit = 1e-2/r_new # cm to m

def color_particles(b):
    clump_ids = []  # List to store unique clump IDs

    # Extract clump IDs from bodies
    for i in b:
        clump_id = O.bodies[i].clumpId
        if clump_id not in clump_ids:
            clump_ids.append(clump_id)

    # Assign colors based on clump IDs
    num_clumps = len(clump_ids)
    half_point = num_clumps // 2

    color_map = {}
    for index, clump_id in enumerate(clump_ids):
        if index < half_point:
            color = Vector3(0.1, 0.1, 0.1)  # Black color for the first half
        else:
            color = Vector3(1, 1, 1)  # White color for the second half
        color_map[clump_id] = color

    # Color bodies based on clump IDs
    for i in b:
        clump_id = O.bodies[i].clumpId
        color = color_map[clump_id]
        O.bodies[i].shape.color = color


# Chess piece dimentions come from FIDE guidelines: 
# https://www.fide.com/FIDE/handbook/Standards_of_Chess_Equipment_and_tournament_venue.pdf
# Add Kings
sp = pack.SpherePack()
file = "Clumps/clumps/Chess_King_100_spheres.txt"
King = pack.SpherePack(read_x_y_z_r_intoList(file, Vector3.Zero, 0.0861*unit)) # 9.5 cm
sp.makeClumpCloud(corner1, corner2, [King], num=2, seed=1, periodic=False)
b = sp.toSimulation(material='spheresMat')
scaling_factors["Clumps/stl_files/Chess_King.stl"] = 0.0861*unit
add_pieces("Clumps/stl_files/Chess_King.stl")
color_particles(b)

# Add Queens
sp = pack.SpherePack()
file = "Clumps/clumps/Chess_Queen_100_spheres.txt"
Queen = pack.SpherePack(read_x_y_z_r_intoList(file, Vector3.Zero, 0.0825*unit)) # 8.5 cm
sp.makeClumpCloud(corner1, corner2, [Queen], num=2, seed=2, periodic=False)
b = sp.toSimulation(material='spheresMat')
scaling_factors["Clumps/stl_files/Chess_Queen.stl"] = 0.0825
add_pieces("Clumps/stl_files/Chess_Queen.stl")
color_particles(b)


# Add Rooks
sp = pack.SpherePack()
file = "Clumps/clumps/Chess_Rook_100_spheres.txt"
Rook = pack.SpherePack(read_x_y_z_r_intoList(file, Vector3.Zero, 0.115*unit)) # 5.5 cm
sp.makeClumpCloud(corner1, corner2, [Rook], num=4, seed=3, periodic=False)
b = sp.toSimulation(material='spheresMat')
scaling_factors["Clumps/stl_files/Chess_Rook.stl"] = 0.115
add_pieces("Clumps/stl_files/Chess_Rook.stl")
color_particles(b)


# Add Bishops
sp = pack.SpherePack()
file = "Clumps/clumps/Chess_Bishop_100_spheres.txt"
Bishop = pack.SpherePack(read_x_y_z_r_intoList(file, Vector3.Zero, 0.0424*unit)) # 7 cm
sp.makeClumpCloud(corner1, corner2, [Bishop], num=4, seed=4, periodic=False)
b = sp.toSimulation(material='spheresMat')
scaling_factors["Clumps/stl_files/Chess_Bishop.stl"] = 0.0424
add_pieces("Clumps/stl_files/Chess_Bishop.stl")
color_particles(b)


# Add Knights
sp = pack.SpherePack()
file = "Clumps/clumps/Chess_Knight_100_spheres.txt"
Knight = pack.SpherePack(read_x_y_z_r_intoList(file, Vector3.Zero, 0.12*unit)) # 6 cm
sp.makeClumpCloud(corner1, corner2, [Knight], num=4, seed=5, periodic=False)
b = sp.toSimulation(material='spheresMat')
scaling_factors["Clumps/stl_files/Chess_Knight.stl"] = 0.12
add_pieces("Clumps/stl_files/Chess_Knight.stl")
color_particles(b)


# Add Pawns
sp = pack.SpherePack()
file = "Clumps/clumps/Chess_Pawn_100_spheres.txt"
Pawn = pack.SpherePack(read_x_y_z_r_intoList(file, Vector3.Zero, 0.062*unit)) # 5 cm
sp.makeClumpCloud(corner1, corner2, [Pawn], num=16, seed=8, periodic=False)
b = sp.toSimulation(material='spheresMat')
scaling_factors["Clumps/stl_files/Chess_Pawn.stl"] = 0.062
add_pieces("Clumps/stl_files/Chess_Pawn.stl")
color_particles(b)

O.bodies.updateClumpProperties(discretization = 45)

# Add Initial velocity
Vel = Vel*s_new/r_new
for b in O.bodies:
    if isinstance(b.shape, Clump):
        theta = random.uniform(0, np.pi)
        phi = random.uniform(0, np.pi)
        vv = Vector3(np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta))
        R = b.state.ori.toRotationMatrix()
        b.state.vel = Vel*vv          # m / s
        b.state.angVel = np.pi*vv   # rads / s
        b.state.angMom = (R*b.state.inertia.asDiagonal()*R.transpose())*b.state.angVel


#------------------------------------------------------------------------------
# Time parameters

# Set time step
O.stopAtTime = t_max/s_new 
O.dt = DT/s_new
steps = int(O.stopAtTime/O.dt)
vis_steps = int(steps/samples)

print("PWaveTimeStep = ", PWaveTimeStep()*s_new)
print("RayleighWaveTimeStep = ", RayleighWaveTimeStep()*s_new)


#------------------------------------------------------------------------------
directory = "paraview"
if paraview:
    if os.path.exists(directory):
        # Remove all files and subdirectories inside the directory
        for filename in os.listdir(directory):
            file_path = os.path.join(directory, filename)
            if os.path.isfile(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                os.rmdir(file_path)
    else:
        # Create the directory
        os.makedirs(directory)

#------------------------------------------------------------------------------
# Engines
O.wait()
O.trackEnergy=True
O.engines=[
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Facet_Aabb()]),
        InteractionLoop(
            [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
            [Ip2_FrictMat_FrictMat_MindlinPhys(betan=0.0, betas=0.0)],
            [Law2_ScGeom_MindlinPhys_Mindlin(calcEnergy=True, label='law2')]
        ),
        NewtonIntegrator(damping=0, gravity=[0,0,0], kinSplit=True,
            exactAsphericalRot=True, rotAlgorithm=integrationType),
        PyRunner(command='addPlotData()', iterPeriod=vis_steps),
        #PyRunner(command='saveSTLs()', iterPeriod=vis_steps),
        VTKRecorder(fileName='paraview/temp-', recorders=['all'], iterPeriod=vis_steps, dead=not paraview)]

#------------------------------------------------------------------------------
# Initial energy
O.step()
O.step()
E0 = utils.kineticEnergy() + O.engines[2].lawDispatcher.functors[0].normElastEnergy() 
O.wait()

#------------------------------------------------------------------------------
# Plot energies
plot.plots={'t':('total', 'kinTrans', 'kinRot', 'Elastic'), 'tt':('Energy_Error')}
plot.plot(subPlots=True)

def addPlotData(): 
    normElastEnergy = O.engines[2].lawDispatcher.functors[0].normElastEnergy() 
    plot.addData(t=O.time*s_new, tt=O.time*s_new, total=(utils.kineticEnergy() + normElastEnergy)*(E_new*r_new**3), kinTrans=O.energy['kinTrans']*(E_new*r_new**3), 
        kinRot=O.energy['kinRot']*(E_new*r_new**3), Elastic=normElastEnergy*(E_new*r_new**3), Energy_Error=Energy_Error())

#------------------------------------------------------------------------------
# Function: Start the user interface if YADE is ran with GUI (i.e. if "yadedaily -n program.py" is not used)
try:
    from yade import qt
    rndr = yade.qt.Renderer()
    v = qt.View()
except: pass

O.saveTmp()
O.run()