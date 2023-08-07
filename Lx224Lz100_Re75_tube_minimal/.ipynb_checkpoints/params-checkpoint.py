import numpy as np
import sys

#######################
# Physical Parameters #
#######################
# Reynolds number
Re = 66.0

# Turbulent Reynolds numner
C_zt = (0.175/66)*Re**2
C_yt = 0.1*(Re/67)**2

# Domain size
Lx = 112.0
Lz = 50.0

# Beta
beta = np.pi/2

# Drag
alpha = 0.02
# Extra drag for q1
kexp = -0.25
kappa = 6.0/67.*(Re/67)**(kexp)

# Closures
a    = 0.8
b = 0.9
eta  = 2.e-4
eps1 = 0.164*2
eps2 = 0.0 

# Kinetic energy of initial condition:
E0 = 1e-4  # Large scale flow
Eq0 = 0.05 # turbulence
nICx = 8 # Mode number cutoff for noisy IC
nICz = 8 # Mode number cutoff for noisy IC

# IC file (if None, then blank restart. Otherwise specify a txt file which is an array for q only)
#IC_file=None
IC_file = 'gauss'
#IC_file='equil'
#IC_file = '../q0_Lz256.txt'

# Random force
rand_force = False
mu = 0.0
sigma = 1.0

########################
# Numerical Parameters #
########################
# Numerical resolution
Nx = int(Lx*2)
Nz = int(Lz*2)
dt_init=(Lz/Nz)**2/0.08/C_zt/8

# Set how long you want to run for:
sim_tmax = np.inf  # simulation time units
#sim_tmax = 3000  # simulation time units
real_tmax = (7+50/60.)*60*60 # 12*60= 12 mins. Real time is in seconds ... 12 hours = 43200
#real_tmax = (3+55/60.)*60*60 # 12*60= 12 mins. Real time is in seconds ... 12 hours = 43200
#real_tmax = (14/60.)*60*60 # 12*60= 12 mins. Real time is in seconds ... 12 hours = 43200
# real_tmax = np.inf
# Simulation will stop when the first of either is reached.

# Real (wall) time interval between snapshots exporting
#tsnap_wall = real_tmax/3.-15*60 # So that we save 3 snapshots 15 minutes before each third
#tsnap_wall = real_tmax/500.
tsnap_wall = np.inf

# Simulation time interval between snapshots exporting
tsnap_sim = 150.

# Time steps between outputting spectra and fluxes:
#tspec = 230
#tspec = 100

# Time steps between outputting time series info
tseries = 50.
#tseries = 10
