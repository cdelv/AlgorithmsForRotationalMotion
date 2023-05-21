#---------------------------------------------------------------------
# Chess pieces bouncing inside a box
# 2023 © Vasileios Angelidakis <vasileios.angelidakis@fau.de>
# 2023 © Carlos Andrés del Valle <cdelv@unal.edu.co>
# Institute for Multiscale Simulation
#---------------------------------------------------------------------
from yade import plot
import random
random.seed(12)
import numpy as np
import csv
import os

readParamsFromTable(dt=1e-5, integration_type="Carlos_2023")
from yade.params import table 

#---------------------------------------------------------------------
BATCH = True
t_max = 1.5          # s
Vel = 1.0            # m/s
samples = 250

#------------------------------------------------------------------------------
# New units for the simulation (to avoid numerical problems)
E_new = 1e9    # new unit of pressure
r_new = 1e-2   # new unit of lenght
s_new = 1.9311703796025776e6 # new unit of time (Rayleigh timestep)^-1

#------------------------------------------------------------------------------
Energy_Err = []
def Energy_Error():
    Ef = utils.kineticEnergy() + O.engines[2].lawDispatcher.functors[0].normElastEnergy() 
    Energy_Err.append(100*abs(E0 - Ef)/E0)
    return sum(Energy_Err)/len(Energy_Err)

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
matID=O.materials.append(FrictMat(young=60e9/E_new, poisson=0.25, density=3200*(r_new**2/(E_new*s_new**2)), frictionAngle=0.0, label='spheresMat'))

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

# Chess piece dimentions come from FIDE guidelines: 
# https://www.fide.com/FIDE/handbook/Standards_of_Chess_Equipment_and_tournament_venue.pdf
# Add Kings
sp = pack.SpherePack()
file = "Clumps/clumps/Chess_King_100_spheres.txt"
King = pack.SpherePack(read_x_y_z_r_intoList(file, Vector3.Zero, 0.0861*unit)) # 9.5 cm
sp.makeClumpCloud(corner1, corner2, [King], num=2, seed=1, periodic=False)
sp.toSimulation(material='spheresMat')

# Add Queens
sp = pack.SpherePack()
file = "Clumps/clumps/Chess_Queen_100_spheres.txt"
Queen = pack.SpherePack(read_x_y_z_r_intoList(file, Vector3.Zero, 0.0825*unit)) # 8.5 cm
sp.makeClumpCloud(corner1, corner2, [Queen], num=2, seed=2, periodic=False)
sp.toSimulation(material='spheresMat')

# Add Rooks
sp = pack.SpherePack()
file = "Clumps/clumps/Chess_Rook_100_spheres.txt"
Rook = pack.SpherePack(read_x_y_z_r_intoList(file, Vector3.Zero, 0.115*unit)) # 5.5 cm
sp.makeClumpCloud(corner1, corner2, [Rook], num=4, seed=3, periodic=False)
sp.toSimulation(material='spheresMat')

# Add Bishops
sp = pack.SpherePack()
file = "Clumps/clumps/Chess_Bishop_100_spheres.txt"
Bishop = pack.SpherePack(read_x_y_z_r_intoList(file, Vector3.Zero, 0.0424*unit)) # 7 cm
sp.makeClumpCloud(corner1, corner2, [Bishop], num=4, seed=4, periodic=False)
sp.toSimulation(material='spheresMat')

# Add Knights
sp = pack.SpherePack()
file = "Clumps/clumps/Chess_Knight_100_spheres.txt"
Knight = pack.SpherePack(read_x_y_z_r_intoList(file, Vector3.Zero, 0.12*unit)) # 6 cm
sp.makeClumpCloud(corner1, corner2, [Knight], num=4, seed=5, periodic=False)
sp.toSimulation(material='spheresMat')

# Add Pawns
sp = pack.SpherePack()
file = "Clumps/clumps/Chess_Pawn_100_spheres.txt"
Pawn = pack.SpherePack(read_x_y_z_r_intoList(file, Vector3.Zero, 0.062*unit)) # 5 cm
sp.makeClumpCloud(corner1, corner2, [Pawn], num=16, seed=8, periodic=False)
sp.toSimulation(material='spheresMat')

O.bodies.updateClumpProperties(discretization=100)

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
O.dt = table.dt/s_new
steps = int(O.stopAtTime/O.dt)
vis_steps = int(steps/samples)


#------------------------------------------------------------------------------
# Engines

O.engines=[
        ForceResetter(),
        InsertionSortCollider([Bo1_Sphere_Aabb(), Bo1_Facet_Aabb()]),
        InteractionLoop(
            [Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom(label='ig2')],
            [Ip2_FrictMat_FrictMat_MindlinPhys(betan=0.0, betas=0.0, label='ip2')],
            [Law2_ScGeom_MindlinPhys_Mindlin(label='law', calcEnergy=True)]
        ),
        NewtonIntegrator(damping=0.0, gravity=[0,0,0], label='newton', kinSplit=True,
            exactAsphericalRot=True, RotAlgorithm=table.integration_type),
        PyRunner(command='Energy_Error()', iterPeriod=vis_steps),
        ]

#------------------------------------------------------------------------------
# Initial energy
O.step()
O.step()
E0 = utils.kineticEnergy() + O.engines[2].lawDispatcher.functors[0].normElastEnergy() 

#------------------------------------------------------------------------------
# Run simulation
O.run()
if BATCH:
    O.wait()

#------------------------------------------------------------------------------

if BATCH:
    # Save results
    fileE = "Data/"+table.integration_type+"_energy.csv"
    if not os.path.exists(fileE):
        with open(fileE, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['dt','mean_error','std_error'])

    with open(fileE, 'a') as f:
        writer = csv.writer(f)
        row = [str(table.dt), str(np.mean(Energy_Err)), str(np.std(Energy_Err))]   
        writer.writerow(row)
else:
    print("dt = ", table.dt, ", err = ", np.mean(Energy_Err), " +- ", np.std(Energy_Err))