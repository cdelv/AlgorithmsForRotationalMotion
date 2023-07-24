#---------------------------------------------------------------------
# Rods bouncing inside a box
# 2023 © Vasileios Angelidakis <vasileios.angelidakis@fau.de>
# 2023 © Carlos Andrés del Valle <cdelv@unal.edu.co>
# Institute for Multiscale Simulation
#---------------------------------------------------------------------
from yade import qt
from yade import plot
import random
import numpy as np
import csv
O.trackEnergy=True

# PWaveTimeStep step fraction
DT = 2.0
samples = 500

# Integration Method
integrationType='Fincham1992'
integrationType='Omelyan1998'
integrationType='delValle2023'

# Material
matID=O.materials.append(FrictMat(young=1.0e9, poisson=0.1, density=2500.0, frictionAngle=0.0, label='spheresMat'))

# Facet box
id1 = O.bodies.append(geom.facetBox((Vector3(0.5,0.5,0.5)), (Vector3(1.0, 1.0, 1.0)), material=O.materials['spheresMat']))

# Create a sphere cloud
sp = pack.SpherePack()
sp.makeClumpCloud((0, 0, 0), (1, 1, 1), 
    [
    pack.SpherePack(
            [
            ((-0.04, -0.04, -0.04), 0.02*sqrt(3)), 
            ((  0.0,   0.0,   0.0), 0.02*sqrt(3)), 
            (( 0.04,  0.04,  0.04), 0.02*sqrt(3))
            ]
        )
    ], 
    periodic=False, num=30, seed=11)

# Add the sphere pack to the simulation
sp.toSimulation(material='spheresMat')

# Data files
"""
fileE = "Data/"+integrationType+"_Energy_"+str(DT)+".csv"
with open(fileE, 'w') as f:
    f.write('t,total,kinTrans,kinRot\n')

fileC = "Data/Clumps.csv"
with open(fileC, 'w') as f:
    f.write('x,y,z,vx,vy,vz,wx,wy,wz,qw,q1,q2,q3\n')
"""

# Ensure that all particles have the same kinetic energy
Vel = 6.0
random.seed(11)
for b in O.bodies:
    if isinstance(b.shape,Clump):
        theta = random.random()
        phi = random.random()
        vv = Vector3(np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta))
        b.state.vel=Vel*vv
        b.state.angVel=0.5*Vel*vv

for b in O.bodies:
    if isinstance(b.shape,Clump):
        x = b.state.pos
        v = b.state.vel
        w = b.state.angVel
        q = b.state.ori
        I = b.state.inertia
        m = b.state.mass
        
        # Save data
        """
        with open(fileC, 'a') as f:
            writer = csv.writer(f)
            data = [x[0], x[1], x[2], v[0], v[1], v[2], w[0], w[1], w[2], q[3], q[0], q[1], q[2]]   
            row = list(map(lambda t: ("%.16f" % t), data)) # 16 decimal places
            writer.writerow(row)
        """

# Set time step
O.stopAtTime = 250.0
O.dt = DT*PWaveTimeStep()
print('dt = ', O.dt)


steps = int(O.stopAtTime/O.dt)
vis_steps = int(steps/samples)

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
            exactAsphericalRot=True, rotAlgorithm=integrationType),
        VTKRecorder(fileName='temp-', recorders=['all'], iterPeriod=500, dead=True),
        PyRunner(command='addPlotData()', iterPeriod=vis_steps)
]

# Plot energies
plot.plots={'t':('total', 'kinTrans', 'kinRot', 'Elastic')}
plot.plot(subPlots=False)

def addPlotData(): 
    normElastEnergy = O.engines[2].lawDispatcher.functors[0].normElastEnergy() 

    plot.addData(t=O.time, total=O.energy['kinTrans'] + O.energy['kinRot'] + normElastEnergy, **O.energy, Elastic=normElastEnergy)

    # Save data
    """
    with open(fileE, 'a') as f:
        writer = csv.writer(f)
        data = [O.time, O.energy.total() + normElastEnergy, O.energy['kinTrans'], O.energy['kinRot']]   
        row = list(map(lambda t: ("%.16f" % t), data)) # 16 decimal places
        writer.writerow(row)
    """

    
#------------------------------------------------------------------------------
# Function: Start the user interface if YADE is ran with GUI (i.e. if "yadedaily -n clumpBouncingOnFacetBox.py" is not used)
try:
    from yade import qt
#    rndr = yade.qt.Renderer()
#    v=qt.View()
except: pass

O.run()