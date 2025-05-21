#
# https://scicomp.stackexchange.com/questions/36013/numerical-computation-of-lyapunov-exponent
#
import numpy as np
from scipy.integrate import solve_ivp

# parameters
a = 0.6
b = 0.3
c = 0.5


# Para el sistema de Barati 2010 las ecuaciones son
# \dot{x} = az
# \dot{y} = z(by² + cxz)
# \dot{z} = y² - 1 - xyz

# Su Jacobiana es 

# J[0] = az' 	
# J[1] = cz²x' + 2byzy' + (by²+2cxz)z'
# J[2] = -yzx' + (2y - xz)y' - xyz'

# ODE system
def func(t, v, sigma, rho, beta ) :
	x, y, z = v
	return [ sigma*(y-x), x*(rho-z)-y, x*y - beta*z]

# Jacobian matrix
def JM(t, v, sigma, rho, beta) :
	x, y, z = v
	return np.array([[ -sigma, sigma, 0], [rho-z, -1, -x], [y, x, -beta]])


# dimension of the system (to keep things general)
n = 3

# number of Lyapunov exponents
n_lyap = 3

# (dis)assemble state and tangent vectors (or their derivatives)
# into one vector for integrator:
assemble = lambda v,U: np.hstack((v,U.flatten()))
disassemble = lambda state: (state[:n], state[n:].reshape(n,n_lyap))

def func_with_lyaps(t, state, *pars):
	v, U = disassemble(state)
	dv = func(t, v, *pars)
	dU = JM(t, v, *pars) @ U
	return assemble(dv, dU)

# initial states:
v = np.random.random(n)
U = np.random.random((n_lyap,n))

v = np.array( [10,0,30] )
U = np.array( [ [1,0,0], [0,1,0], [0,0,1] ] )
# U = np.random.random((n_lyap,n))

lyaps = [] # local Lyapunov exponents

dt = 1
iters = 100

# 1 dt for transient time
sol = solve_ivp(
	func_with_lyaps,
	[0, dt],
	assemble(v,U),
	t_eval=[dt],
	args=(sigma, rho, beta),
	max_step=dt,
)
v,U = disassemble(sol.y.flatten())
U,R = np.linalg.qr(U)

exps = np.zeros( (3,) )
for _ in range(iters):
	sol = solve_ivp(
			func_with_lyaps,
			[0, dt],
			assemble(v,U),
			t_eval=[dt],
			args=(sigma, rho, beta),
			max_step=dt,
		)
	v,U = disassemble(sol.y.flatten())
	U,R = np.linalg.qr(U)
	w = np.log(abs(R.diagonal()))/dt
	exps += w
	t = dt + _ * dt
	av = exps/t
	print( t, av[0], av[1], av[2] )
	lyaps.append( w )

transient_steps = 50
# result:
print( "# ", *np.average(lyaps[transient_steps:],axis=0) )
