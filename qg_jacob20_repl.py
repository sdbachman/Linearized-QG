#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 10:24:38 2017

@author: jacob
"""

"""
Dedalus script for Finite Amplitude Calcs

To run, merge, and plot using 4 processes, for instance, you could use:
    $ mpiexec -n 4 python3 rayleigh_benard.py
    $ mpiexec -n 4 python3 merge.py snapshots
    $ mpiexec -n 4 python3 plot_2d_series.py snapshots/*.h5
    $ mpiexec -n 2 python3 qg_bathy_w_param.py
"""

import numpy as np
from mpi4py import MPI
CW = MPI.COMM_WORLD
import time
from pylab import *
from dedalus import public as de
from dedalus.extras import flow_tools
import scipy.integrate as integrate
import logging
logger = logging.getLogger(__name__)
import subprocess

# Parameters
#directoryname = '/scratch/jacob13/NLSIM/'
#ly_global = np.logspace(-5, -3, 192)*2*np.pi
#OD = False



nx = 512
ny = 512

mesh = [4,4]


f = 1e-4 # Coriolis parameter
#beta = 2e-11 # Beta Parameter
A4 = 1e-11
H1 = 500.0
H2 = 500.0

#Magnitude of shear flow
U1m = 0.0247
U2m = 0.0234
U3m = 0.0221
U4m = 0.0208
U5m = 0.0195
U6m = 0.0182
U7m = 0.0169
U8m = 0.0156
U9m = 0.0143
U10m = 0.0130
U11m = 0.0117
U12m = 0.0104
U13m = 0.0091
U14m = 0.0078
U15m = 0.0065
U16m = 0.0052
U17m = 0.0039
U18m = 0.0026
U19m = 0.0013
U20m = 0.000
#print(gp1)
#print(Hp)
# Deformation Radius
Ld = 4.5e3     # 3e4 for the three layer case, 4.5e3 for the 20-layer case
Ld_total = Ld * 20

dx = Ld_total / 8.0
dy = dx

#Set the strength of the bottom drag (see Thompson and Young 2006)
REPL1
beta = b_star * U1m / Ld_total / Ld_total
#beta = 1.2e-11

REPL2
#k_star = 0.2
rek = U1m * k_star / Ld_total

#beta = 5.555555555555557e-12/4.0
#rek = 1.6666666666666667e-08

f = open('./PARAMS.txt', 'w')
f.write('beta = ' + str(beta) + '\n')
f.write('rek = ' + str(rek) + '\n')
f.write('\n')
f.close()


#gridratio = 1.7 # gridratio = nx/(2*pi*Lx/Ld)See Thompson and Young  2006

#Lx, Ly = (nx/(gridratio/Ld),ny/(gridratio/Ld))
Lx, Ly = nx*dx, ny*dy

# Visc ratio, see Thompson and Young 2006
vr = 1e-13
L = Lx/(2*np.pi)
A8 = vr*(U1m*L**7)
#L = 10000.0*nx / (2*np.pi)
#A8 = vr*(0.25*L**7)
A8 = 0

#%% 3D PROBLEM
# Create basis and domain
start_init_time = time.time()
x_basis = de.Fourier('x', nx, interval=(0, Lx), dealias=3/2)
y_basis = de.Fourier('y', ny, interval=(0, Ly), dealias=3/2)

domain = de.Domain([x_basis, y_basis], grid_dtype=np.float64, mesh=None)
y = domain.grid(1)
x = domain.grid(0)

#########################################



#H = domain.new_field(name='H') # Total depth
#H['g'] = H0*(x/x+0*y) - eta_b

# Define Fields (these can later define spatially varying velocities)
# Define Fields (these can later define spatially varying velocities)
#U1 = domain.new_field(name='U1')
#U2 = domain.new_field(name='U2')

#U1['g'] = U1m
#U2['g'] = U2m


# set up IVP
problem = de.IVP(domain, variables=['q1', 'q2','q3','q4', 'q5','q6','q7', 'q8','q9','q10', 'q11','q12','q13', 'q14','q15','q16', 'q17','q18','q19', 'q20','psi1', 'psi2','psi3','psi4', 'psi5','psi6','psi7', 'psi8','psi9','psi10', 'psi11','psi12','psi13', 'psi14','psi15','psi16', 'psi17','psi18','psi19', 'psi20', 'u1', 'v1'])
#,'v2','v3','v4','v5','v6','v7','v8','v9','v10', 'v11','v12','v13','v14','v15','v16','v17','v18','v19','v20'])

slices = domain.dist.grid_layout.slices(scales=1)

problem.parameters['delta'] = H1/H2
problem.parameters['DeltaU'] = U1m - U2m
problem.parameters['beta'] = beta
problem.parameters['Ld'] = Ld
problem.parameters['U1m'] = U1m
problem.parameters['U2m'] = U2m
problem.parameters['U3m'] = U3m
problem.parameters['U4m'] = U4m
problem.parameters['U5m'] = U5m
problem.parameters['U6m'] = U6m
problem.parameters['U7m'] = U7m
problem.parameters['U8m'] = U8m
problem.parameters['U9m'] = U9m
problem.parameters['U10m'] = U10m
problem.parameters['U11m'] = U11m
problem.parameters['U12m'] = U12m
problem.parameters['U13m'] = U13m
problem.parameters['U14m'] = U14m
problem.parameters['U15m'] = U15m
problem.parameters['U16m'] = U16m
problem.parameters['U17m'] = U17m
problem.parameters['U18m'] = U18m
problem.parameters['U19m'] = U19m
problem.parameters['U20m'] = U20m
problem.parameters['A4'] = A4
problem.parameters['A8'] = A8
problem.parameters['rek'] = rek

#Following notation in Arbic and Flierl
# define substitutions
problem.substitutions['L4(A)'] = '(dx(dx(dx(dx(A)))) + 2*dx(dx(dy(dy(A)))) + dy(dy(dy(dy(A)))))' #Horizontal biharmonic diff
problem.substitutions['HV(A)'] = '-A8*(L4(dx(dx(dx(dx(A))))) + L4(2*dx(dx(dy(dy(A))))) + L4(dy(dy(dy(dy(A))))))'
problem.substitutions['Jac(A,B)'] = 'dx(A)*dy(B) - dy(A)*dx(B) ' # Jacobian
problem.substitutions['Lap(A)'] = 'dx(dx(A)) + dy(dy(A))'
  
problem.substitutions['F'] = '1/(Ld**2)'
#problem.substitutions['F2'] = 'delta/((delta+1)*Ld**2)'
problem.substitutions['Qy1'] = 'beta + F*DeltaU'
problem.substitutions['Qy2'] = 'beta'
problem.substitutions['Qy3'] = 'beta'
problem.substitutions['Qy4'] = 'beta'
problem.substitutions['Qy5'] = 'beta'
problem.substitutions['Qy6'] = 'beta'
problem.substitutions['Qy7'] = 'beta'
problem.substitutions['Qy8'] = 'beta'
problem.substitutions['Qy9'] = 'beta'
problem.substitutions['Qy10'] = 'beta'
problem.substitutions['Qy11'] = 'beta'
problem.substitutions['Qy12'] = 'beta'
problem.substitutions['Qy13'] = 'beta'
problem.substitutions['Qy14'] = 'beta'
problem.substitutions['Qy15'] = 'beta'
problem.substitutions['Qy16'] = 'beta'
problem.substitutions['Qy17'] = 'beta'
problem.substitutions['Qy18'] = 'beta'
problem.substitutions['Qy19'] = 'beta'
problem.substitutions['Qy20'] = 'beta - F*DeltaU'

# Calculate background PV gradient using this:
# https://pyqg.readthedocs.io/en/latest/equations/notation_layered.html

# define equations
problem.add_equation('dt(q1)  + dx(psi1)*Qy1 + U1m*dx(q1) = 0', condition='(nx!=0) or (ny !=0)')
problem.add_equation('dt(q2)  + dx(psi2)*Qy2 + U2m*dx(q2) = 0', condition='(nx!=0) or (ny !=0)')
problem.add_equation('dt(q3)  + dx(psi3)*Qy3 + U3m*dx(q3) = 0', condition='(nx!=0) or (ny !=0)')
problem.add_equation('dt(q4)  + dx(psi4)*Qy4 + U4m*dx(q4) = 0', condition='(nx!=0) or (ny !=0)')
problem.add_equation('dt(q5)  + dx(psi5)*Qy5 + U5m*dx(q5) = 0', condition='(nx!=0) or (ny !=0)')
problem.add_equation('dt(q6)  + dx(psi6)*Qy6 + U6m*dx(q6) = 0', condition='(nx!=0) or (ny !=0)')
problem.add_equation('dt(q7)  + dx(psi7)*Qy7 + U7m*dx(q7) = 0', condition='(nx!=0) or (ny !=0)')
problem.add_equation('dt(q8)  + dx(psi8)*Qy8 + U8m*dx(q8) = 0', condition='(nx!=0) or (ny !=0)')
problem.add_equation('dt(q9)  + dx(psi9)*Qy9 + U9m*dx(q9) = 0', condition='(nx!=0) or (ny !=0)')
problem.add_equation('dt(q10)  + dx(psi10)*Qy10 + U10m*dx(q10) = 0', condition='(nx!=0) or (ny !=0)')
problem.add_equation('dt(q11)  + dx(psi11)*Qy11 + U11m*dx(q11) = 0', condition='(nx!=0) or (ny !=0)')
problem.add_equation('dt(q12)  + dx(psi12)*Qy12 + U12m*dx(q12) = 0', condition='(nx!=0) or (ny !=0)')
problem.add_equation('dt(q13)  + dx(psi13)*Qy13 + U13m*dx(q13) = 0', condition='(nx!=0) or (ny !=0)')
problem.add_equation('dt(q14)  + dx(psi14)*Qy14 + U14m*dx(q14) = 0', condition='(nx!=0) or (ny !=0)')
problem.add_equation('dt(q15)  + dx(psi15)*Qy15 + U15m*dx(q15) = 0', condition='(nx!=0) or (ny !=0)')
problem.add_equation('dt(q16)  + dx(psi16)*Qy16 + U16m*dx(q16) = 0', condition='(nx!=0) or (ny !=0)')
problem.add_equation('dt(q17)  + dx(psi17)*Qy17 + U17m*dx(q17) = 0', condition='(nx!=0) or (ny !=0)')
problem.add_equation('dt(q18)  + dx(psi18)*Qy18 + U18m*dx(q18) = 0', condition='(nx!=0) or (ny !=0)')
problem.add_equation('dt(q19)  + dx(psi19)*Qy19 + U19m*dx(q19) = 0', condition='(nx!=0) or (ny !=0)')
problem.add_equation('dt(q20)  + dx(psi20)*Qy20 + U20m*dx(q20) = 0', condition='(nx!=0) or (ny != 0)')

problem.add_equation('psi1=0', condition='(nx==0) and (ny==0)')
problem.add_equation('psi2=0', condition='(nx==0) and (ny==0)')
problem.add_equation('psi3=0', condition='(nx==0) and (ny==0)')
problem.add_equation('psi4=0', condition='(nx==0) and (ny==0)')
problem.add_equation('psi5=0', condition='(nx==0) and (ny==0)')
problem.add_equation('psi6=0', condition='(nx==0) and (ny==0)')
problem.add_equation('psi7=0', condition='(nx==0) and (ny==0)')
problem.add_equation('psi8=0', condition='(nx==0) and (ny==0)')
problem.add_equation('psi9=0', condition='(nx==0) and (ny==0)')
problem.add_equation('psi10=0', condition='(nx==0) and (ny==0)')
problem.add_equation('psi11=0', condition='(nx==0) and (ny==0)')
problem.add_equation('psi12=0', condition='(nx==0) and (ny==0)')
problem.add_equation('psi13=0', condition='(nx==0) and (ny==0)')
problem.add_equation('psi14=0', condition='(nx==0) and (ny==0)')
problem.add_equation('psi15=0', condition='(nx==0) and (ny==0)')
problem.add_equation('psi16=0', condition='(nx==0) and (ny==0)')
problem.add_equation('psi17=0', condition='(nx==0) and (ny==0)')
problem.add_equation('psi18=0', condition='(nx==0) and (ny==0)')
problem.add_equation('psi19=0', condition='(nx==0) and (ny==0)')
problem.add_equation('psi20=0', condition='(nx==0) and (ny==0)')

problem.add_equation('-q1 + Lap(psi1) +F*(psi2-psi1)=0')
problem.add_equation('-q2 + Lap(psi2) +F*(psi1-2*psi2+psi3)=0' )
problem.add_equation('-q3 + Lap(psi3) +F*(psi2-2*psi3+psi4)=0' )
problem.add_equation('-q4 + Lap(psi4) +F*(psi3-2*psi4+psi5)=0' )
problem.add_equation('-q5 + Lap(psi5) +F*(psi4-2*psi5+psi6)=0' )
problem.add_equation('-q6 + Lap(psi6) +F*(psi5-2*psi6+psi7)=0' )
problem.add_equation('-q7 + Lap(psi7) +F*(psi6-2*psi7+psi8)=0' )
problem.add_equation('-q8 + Lap(psi8) +F*(psi7-2*psi8+psi9)=0' )
problem.add_equation('-q9 + Lap(psi9) +F*(psi8-2*psi9+psi10)=0' )
problem.add_equation('-q10 + Lap(psi10) +F*(psi9-2*psi10+psi11)=0' )
problem.add_equation('-q11 + Lap(psi11) +F*(psi10-2*psi11+psi12)=0' )
problem.add_equation('-q12 + Lap(psi12) +F*(psi11-2*psi12+psi13)=0' )
problem.add_equation('-q13 + Lap(psi13) +F*(psi12-2*psi13+psi14)=0' )
problem.add_equation('-q14 + Lap(psi14) +F*(psi13-2*psi14+psi15)=0' )
problem.add_equation('-q15 + Lap(psi15) +F*(psi14-2*psi15+psi16)=0' )
problem.add_equation('-q16 + Lap(psi16) +F*(psi15-2*psi16+psi17)=0' )
problem.add_equation('-q17 + Lap(psi17) +F*(psi16-2*psi17+psi18)=0' )
problem.add_equation('-q18 + Lap(psi18) +F*(psi17-2*psi18+psi19)=0' )
problem.add_equation('-q19 + Lap(psi19) +F*(psi18-2*psi19+psi20)=0' )
problem.add_equation('-q20 + Lap(psi20) +F*(psi19-psi20)=0' )

problem.add_equation('u1 + dy(psi1) = 0') #Necessary only for the CFL condition
problem.add_equation('v1 - dx(psi1) = 0')
#problem.add_equation('v2 - dx(psi2) = 0')
#problem.add_equation('v3 - dx(psi3) = 0')
#problem.add_equation('v4 - dx(psi4) = 0')
#problem.add_equation('v5 - dx(psi5) = 0')
#problem.add_equation('v6 - dx(psi6) = 0')
#problem.add_equation('v7 - dx(psi7) = 0')
#problem.add_equation('v8 - dx(psi8) = 0')
#problem.add_equation('v9 - dx(psi9) = 0')
#problem.add_equation('v10 - dx(psi10) = 0')
#problem.add_equation('v11 - dx(psi11) = 0')
#problem.add_equation('v12 - dx(psi12) = 0')
#problem.add_equation('v13 - dx(psi13) = 0')
#problem.add_equation('v14 - dx(psi14) = 0')
#problem.add_equation('v15 - dx(psi15) = 0')
#problem.add_equation('v16 - dx(psi16) = 0')
#problem.add_equation('v17 - dx(psi17) = 0')
#problem.add_equation('v18 - dx(psi18) = 0')
#problem.add_equation('v19 - dx(psi19) = 0')
#problem.add_equation('v20 - dx(psi20) = 0')
# Build solver
solver = problem.build_solver(de.timesteppers.RK443) 
logger.info('Solver built')

#%%
# define initial condtions
q1 = solver.state['q1']
q2 = solver.state['q2']
q3 = solver.state['q3']
q4 = solver.state['q4']
q5 = solver.state['q5']
q6 = solver.state['q6']
q7 = solver.state['q7']
q8 = solver.state['q8']
q9 = solver.state['q9']
q10 = solver.state['q10']
q11 = solver.state['q11']
q12 = solver.state['q12']
q13 = solver.state['q13']
q14 = solver.state['q14']
q15 = solver.state['q15']
q16 = solver.state['q16']
q17 = solver.state['q17']
q18 = solver.state['q18']
q19 = solver.state['q19']
q20 = solver.state['q20']

psi1 = solver.state['psi1']
psi2 = solver.state['psi2']
psi3 = solver.state['psi3']
psi4 = solver.state['psi4']
psi5 = solver.state['psi5']
psi6 = solver.state['psi6']
psi7 = solver.state['psi7']
psi8 = solver.state['psi8']
psi9 = solver.state['psi9']
psi10 = solver.state['psi10']
psi11 = solver.state['psi11']
psi12 = solver.state['psi12']
psi13 = solver.state['psi13']
psi14 = solver.state['psi14']
psi15 = solver.state['psi15']
psi16 = solver.state['psi16']
psi17 = solver.state['psi17']
psi18 = solver.state['psi18']
psi19 = solver.state['psi19']
psi20 = solver.state['psi20']

#v2 = solver.state['v2']
#v3 = solver.state['v3']
#v4 = solver.state['v4']
#v5 = solver.state['v5']
#v6 = solver.state['v6']
#v7 = solver.state['v7']
#v8 = solver.state['v8']
#v9 = solver.state['v9']
#v10 = solver.state['v10']
#v11 = solver.state['v11']
#v12 = solver.state['v12']
#v13 = solver.state['v13']
#v14 = solver.state['v14']
#v15 = solver.state['v15']
#v16 = solver.state['v16']
#v17 = solver.state['v17']
#v18 = solver.state['v18']
#v19 = solver.state['v19']
#v20 = solver.state['v20']
# Random perturbations, initialized globally for same results in parallel
gshape = domain.dist.grid_layout.global_shape(scales=1)
slices = domain.dist.grid_layout.slices(scales=1)
rand = np.random.RandomState(seed=23)
noise = rand.standard_normal(gshape)[slices]  

q1['g'] = 0
q1['g'] +=1e-7*noise
rand = np.random.RandomState(seed=24)
noise = rand.standard_normal(gshape)[slices] 
q2['g'] = 1e-7*noise
rand = np.random.RandomState(seed=25)
noise = rand.standard_normal(gshape)[slices]
q3['g'] = 1e-7*noise
rand = np.random.RandomState(seed=26)
noise = rand.standard_normal(gshape)[slices]
q4['g'] = 1e-7*noise
rand = np.random.RandomState(seed=27)
noise = rand.standard_normal(gshape)[slices]
q5['g'] = 1e-7*noise
rand = np.random.RandomState(seed=28)
noise = rand.standard_normal(gshape)[slices]
q6['g'] = 1e-7*noise
rand = np.random.RandomState(seed=29)
noise = rand.standard_normal(gshape)[slices]
q7['g'] = 1e-7*noise
rand = np.random.RandomState(seed=30)
noise = rand.standard_normal(gshape)[slices]
q8['g'] = 1e-7*noise
rand = np.random.RandomState(seed=31)
noise = rand.standard_normal(gshape)[slices]
q9['g'] = 1e-7*noise
rand = np.random.RandomState(seed=32)
noise = rand.standard_normal(gshape)[slices]
q10['g'] = 1e-7*noise
rand = np.random.RandomState(seed=33)
noise = rand.standard_normal(gshape)[slices]
q11['g'] = 1e-7*noise
rand = np.random.RandomState(seed=34)
noise = rand.standard_normal(gshape)[slices]
q12['g'] = 1e-7*noise
rand = np.random.RandomState(seed=35)
noise = rand.standard_normal(gshape)[slices]
q13['g'] = 1e-7*noise
rand = np.random.RandomState(seed=36)
noise = rand.standard_normal(gshape)[slices]
q14['g'] = 1e-7*noise
rand = np.random.RandomState(seed=37)
noise = rand.standard_normal(gshape)[slices]
q15['g'] = 1e-7*noise
rand = np.random.RandomState(seed=38)
noise = rand.standard_normal(gshape)[slices]
q16['g'] = 1e-7*noise
rand = np.random.RandomState(seed=39)
noise = rand.standard_normal(gshape)[slices]
q17['g'] = 1e-7*noise
rand = np.random.RandomState(seed=40)
noise = rand.standard_normal(gshape)[slices]
q18['g'] = 1e-7*noise
rand = np.random.RandomState(seed=41)
noise = rand.standard_normal(gshape)[slices]
q19['g'] = 1e-7*noise
rand = np.random.RandomState(seed=42)
noise = rand.standard_normal(gshape)[slices]
q20['g'] = 1e-7*noise
rand = np.random.RandomState(seed=43)
noise = rand.standard_normal(gshape)[slices]

psi1['g']=0
psi2['g']=0
psi3['g']=0
psi4['g']=0
psi5['g']=0
psi6['g']=0
psi7['g']=0
psi8['g']=0
psi9['g']=0
psi10['g']=0
psi11['g']=0
psi12['g']=0
psi13['g']=0
psi14['g']=0
psi15['g']=0
psi16['g']=0
psi17['g']=0
psi18['g']=0
psi19['g']=0
psi20['g']=0

#v2['g']=0
#v3['g']=0
#v4['g']=0
#v5['g']=0
#v6['g']=0
#v7['g']=0
#v8['g']=0
#v9['g']=0
#v10['g']=0
#v11['g']=0
#v12['g']=0
#v13['g']=0
#v14['g']=0
#v15['g']=0
#v16['g']=0
#v17['g']=0
#v18['g']=0
#v19['g']=0
#v20['g']=0
#%%
# With a 4x4 mesh I get 1000 timesteps per 5 minutes
# If I do a snapshot every 1000 timesteps that is 12 per hour
# I can run for 30 hours on izumi on the verylong queue = 360 snapshots

# I will do a snapshot every 500 timesteps = 24/hour and run for 20 hours

# Integration parameters
#solver.stop_sim_time = 1000*86400
#solver.stop_wall_time = 5000
#solver.stop_iteration = np.inf
solver.stop_wall_time = np.inf
#solver.stop_iteration = 5000
solver.stop_iteration = 20000
Ti = Ld/abs(U1m) #Eddy turnover timescale
#dt = Ti/200.0
dt = 3600.0
#Ti/10000.0 #First guess at time step

#solver.stop_sim_time = 3600*24*4166
solver.stop_sim_time = np.inf
 
# Analysis
#snap = solver.evaluator.add_file_handler('snapshots', sim_dt=3600*24*100, max_writes=24*1000, parallel=False)
snap = solver.evaluator.add_file_handler('snapshots', iter=100, max_writes=24*1000, parallel=False)

# Basic Diagnostics
#snap.add_task('q1', name = 'q1')
snap.add_task('psi1', name = 'psi1')
snap.add_task('psi2', name = 'psi2')
snap.add_task('psi3', name = 'psi3')
snap.add_task('psi4', name = 'psi4')
snap.add_task('psi5', name = 'psi5')
snap.add_task('psi6', name = 'psi6')
snap.add_task('psi7', name = 'psi7')
snap.add_task('psi8', name = 'psi8')
snap.add_task('psi9', name = 'psi9')
snap.add_task('psi10', name = 'psi10')
snap.add_task('psi11', name = 'psi11')
snap.add_task('psi12', name = 'psi12')
snap.add_task('psi13', name = 'psi13')
snap.add_task('psi14', name = 'psi14')
snap.add_task('psi15', name = 'psi15')
snap.add_task('psi16', name = 'psi16')
snap.add_task('psi17', name = 'psi17')
snap.add_task('psi18', name = 'psi18')
snap.add_task('psi19', name = 'psi19')
snap.add_task('psi20', name = 'psi20')
snap.add_task('q1', name = 'q1')
snap.add_task('q2', name = 'q2')
snap.add_task('q3', name = 'q3')
snap.add_task('q4', name = 'q4')
snap.add_task('q5', name = 'q5')
snap.add_task('q6', name = 'q6')
snap.add_task('q7', name = 'q7')
snap.add_task('q8', name = 'q8')
snap.add_task('q9', name = 'q9')
snap.add_task('q10', name = 'q10')
snap.add_task('q11', name = 'q11')
snap.add_task('q12', name = 'q12')
snap.add_task('q13', name = 'q13')
snap.add_task('q14', name = 'q14')
snap.add_task('q15', name = 'q15')
snap.add_task('q16', name = 'q16')
snap.add_task('q17', name = 'q17')
snap.add_task('q18', name = 'q18')
snap.add_task('q19', name = 'q19')
snap.add_task('q20', name = 'q20')
#snap.add_task('u1', name = 'u1')
#snap.add_task('v1', name = 'v1')
#snap.add_task('dx(psi2)', name = 'v2')
#snap.add_task('dx(psi3)', name = 'v3')
#snap.add_task('dx(psi4)', name = 'v4')
#snap.add_task('dx(psi5)', name = 'v5')
#snap.add_task('dx(psi6)', name = 'v6')
#snap.add_task('dx(psi7)', name = 'v7')
#snap.add_task('dx(psi8)', name = 'v8')
#snap.add_task('dx(psi9)', name = 'v9')
#snap.add_task('dx(psi10)', name = 'v10')
#snap.add_task('dx(psi11)', name = 'v11')
#snap.add_task('dx(psi12)', name = 'v12')
#snap.add_task('dx(psi13)', name = 'v13')
#snap.add_task('dx(psi14)', name = 'v14')
#snap.add_task('dx(psi15)', name = 'v15')
#snap.add_task('dx(psi16)', name = 'v16')
#snap.add_task('dx(psi17)', name = 'v17')
#snap.add_task('dx(psi18)', name = 'v18')
#snap.add_task('dx(psi19)', name = 'v19')
#snap.add_task('dx(psi20)', name = 'v20')

#snap.add_task('v2', name = 'v2')
#snap.add_task('v3', name = 'v3')
#snap.add_task('v4', name = 'v4')
#snap.add_task('v5', name = 'v5')
#snap.add_task('v6', name = 'v6')
#snap.add_task('v7', name = 'v7')
#snap.add_task('v8', name = 'v8')
#snap.add_task('v9', name = 'v9')
#snap.add_task('v10', name = 'v10')
#snap.add_task('v11', name = 'v11')
#snap.add_task('v12', name = 'v12')
#snap.add_task('v13', name = 'v13')
#snap.add_task('v14', name = 'v14')
#snap.add_task('v15', name = 'v15')
#snap.add_task('v16', name = 'v16')
#snap.add_task('v17', name = 'v17')
#snap.add_task('v18', name = 'v18')
#snap.add_task('v19', name = 'v19')
#snap.add_task('v20', name = 'v20')

#ssnap = solver.evaluator.add_file_handler('ssnapshots', iter=10, max_writes=24*1000, parallel=False)
#ssnap.add_task('v1', name = 'v1', layout='c')
#ssnap.add_task('u1', name = 'u1', layout='c')


# CFL
CFL = flow_tools.CFL(solver, initial_dt=dt, cadence=1, safety=2,
                     max_change=1.5, min_change=0, max_dt=86400.0) #max_dt=Ti/2.0)
CFL.add_velocities(('u1', 'v1'))



# Main loop
end_init_time = time.time()
logger.info('Initialization time: %f' %(end_init_time-start_init_time))

try:
  logger.info('Starting loop')
  start_run_time = time.time()


  while solver.ok:
    dt = CFL.compute_dt()                                                                                                                                                                                                                                                                                                                                                                                          
    solver.step(dt)

    if (solver.iteration-1) % 500 == 1:
      logger.info('Iteration: %i, Days: %1.1f, dt: %e' %(solver.iteration, solver.sim_time/86400, dt))
      qtemp = solver.state['q1']
      if qtemp['g'].size > 0:
        qm = np.max(np.abs(qtemp['g']))
        logger.info('q1 Val: %f' % qm)
        if np.isnan(qm):
          print('NaN encountered.')
except:
    logger.error('Exception raised, triggering end of main loop.')
    exit()
    raise
finally:
    end_run_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_run_time-start_run_time))
    logger.info('Run time: %f cpu-hr' %((end_run_time-start_run_time)/60/60*domain.dist.comm_cart.size))
