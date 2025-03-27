#---------------------------------------------------------------------
# Rods bouncing inside a box
# 2023 © Vasileios Angelidakis <v.angelidakis@qub.ac.uk>
# 2023 © Carlos Andrés del Valle <cdelv@unal.edu.co>
#---------------------------------------------------------------------
from yade import qt
from yade import plot
from yade import utils
import random
import numpy as np
import csv
import pandas as pd 

O.trackEnergy = True
paraview = False

rotAlgorithm = "Fincham1992"
rotAlgorithm = "Omelyan1998"
#rotAlgorithm = "delValle2023"

data_file = rotAlgorithm+".csv"

# PWaveTimeStep step fraction
DT = 0.7
samples = 20000
Vel = 4.0

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
O.bodies.updateClumpProperties(discretization=100)

# Ensure that all particles have the same kinetic energy
random.seed(11)
for b in O.bodies:
    if isinstance(b.shape,Clump):
        theta = random.random()
        phi = random.random()
        vv = Vector3(np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta))
        b.state.vel=Vel*vv
        b.state.angVel=Vel*vv


# Set time step
O.stopAtTime = 10000.0
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
        NewtonIntegrator(rotAlgorithm=rotAlgorithm, damping=0.0, gravity=[0,0,0], label='newton', kinSplit=True, exactAsphericalRot=True),
        PyRunner(command='addPlotData()', iterPeriod=vis_steps),
        VTKRecorder(fileName='frames/temp-', recorders=['all'], iterPeriod=vis_steps, dead=not paraview)
]

# Plot energies
plot.plots={'t':('total', 'kinTrans', 'kinRot', 'Elastic')}
plot.plot(subPlots=False)

energy = []
time = []

def addPlotData(): 
    normElastEnergy = O.engines[2].lawDispatcher.functors[0].normElastEnergy()
    total_energy = utils.kineticEnergy() + normElastEnergy 
    plot.addData(t=O.time, total=total_energy, **O.energy, Elastic=normElastEnergy)
    energy.append(total_energy)
    time.append(O.time)

    
try:
    from yade import qt
    rndr = yade.qt.Renderer()
    v=qt.View()
except: pass

O.run()
O.wait()


df = pd.DataFrame({'t': time, 'E': energy})
df.to_csv(data_file)