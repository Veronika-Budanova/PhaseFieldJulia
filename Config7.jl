# GRID
L = 2 #1e-3 
N = 101 #51 
fict = 3
h = L / (N - 1)

# PHYSIC PARAMETERS
#b = 100
#R = 0.12
b = 50

rho1A = 90
rho1B = 1200
rho2A = 900
rho2B = 120

lam11 = 4e-8
lam12 = 0
lam21 = lam12
lam22 = lam11

Apsi = 1e-4
Bpsi = Apsi

lamelB = 2e2
muelB = lamelB

theta = 1

M0 = 1

eta = 10
zeta = 0

g = 0.05

u_init = -0.1

bound = "wall"

# TIME PARAMETERS
#CFL = 0.125

t_max = 30
delta_t = 1e-6
step_max = t_max/delta_t

frame_dt = 0.1

# PISTON PARAMETERS
#gamma = 1.4

piston_width = 0.3 * L

x0_piston = 0.6 * L
#v0_piston = 0.0

#S = 1.0
#V0 = (x0_piston - piston_width/2) * S # Начальный объем каждого из газов
#p0 = 100.0

# DIRECTORY
direct = "CSframes/Config7"
