#---------------------------------------------------------------------
# This script is meant to be used to compare with the analytical 
# Solution in the Jupyter notebook for the case Tx = 0.
# 2023 © Carlos Andrés del Valle <cdelv@unal.edu.co>
# Institute for Multiscale Simulation
#---------------------------------------------------------------------
from yade import qt
from yade import plot
import csv

# Initial angular velocity
wx0 = 1.0
wy0 = 1.0
wz0 = 1.0

# Inertia tensor components
Ix = 2.0
Iy = 1.0
Iz = 1.0

# Simulation time
O.stopAtTime = 20.0
O.dt = 0.0001
n_samples = 100
n_steps = int(O.stopAtTime/O.dt)
vis_steps = int(n_steps/n_samples)

# File data
file = "Data/Yade_Analytical_Compare.csv"
with open(file, 'w') as f:
    f.write('t,q0,q1,q2,q3,wx,wy,wz\n')

# Create a sphere
Body = utils.sphere([0,0,0], 1.0)

# Make it non spherical and put it to spin
Body.aspherical = True
Body.state.inertia = Vector3(Ix, Iy, Iz)
Body.state.angVel = Vector3(wx0, wy0, wz0)
Body.state.angMom = Vector3(Ix*wx0, Iy*wy0, Iz*wz0)
Body.state.ori = Quaternion().Identity

# Add it to the simulation
O.bodies.append(Body)

# Create the engine
O.engines=[
    ForceResetter(),
    InsertionSortCollider(),
    InteractionLoop(),
    NewtonIntegrator(damping=0.0, gravity=[0,0,0], exactAsphericalRot=True, Omelyan98=False, Allen89=True),
    PyRunner(iterPeriod=vis_steps, command='addplot()')
]

# Set plot and data saving
def addplot():
    b = O.bodies[0]
    q = b.state.ori
    w = q.toRotationMatrix().transpose()*b.state.angVel # Move to body frame
    plot.addData(
            wx=w[0],
            wy=w[1],
            wz=w[2],
            q0=q[3], # Just so it looks the same as the jupyter
            q1=q[0], # In YADE Quaternion(q3, (q0, q1, q2))
            q2=q[1], # To make the equivalence with the Jupyter notebook
            q3=q[2], # Remember that when you print to the screen is displayed in Angle-Axis form
            t=O.time,
            tt=O.time
        )
    
    # Save data
    with open(file, 'a') as f:
        writer = csv.writer(f)
        data = [O.time, q[3], q[0], q[1], q[2], w[0], w[1], w[2]]   
        row = list(map(lambda t: ("%.16f" % t), data)) # 16 decimal places
        writer.writerow(row)


# Show plot
plot.plots = {'t': ('wx', 'wy', 'wz'), 'tt':('q0', 'q1', 'q2', 'q3')}
plot.plot(subPlots=True)

# Run the simulation
O.run()