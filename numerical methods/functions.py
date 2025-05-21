def lorenz_equations(X, sigma, rho, beta):
    x, y, z = X
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    return [dx, dy, dz]

def moore_spiegel(X, T, R):
    x, y, z = X
    dx = y
    dy = z
    dz = -z-(T-R+R*x**2)*y-T*x
    return [dx, dy, dz]

def whang(X, m, n):
    x, y, z = X
    dx = y
    dy = -x+(y*z)+m
    dz = -x-(15*x*y)-(x*z)+n
    return [dx, dy, dz]

def barati(X, a, b, c):
    x, y, z = X
    dx = a*z
    dy = z*(b*(y**2) + c*x*z)
    dz = y**2 - 1 - x*y*z
    return [dx,dy,dz]

def wei_chen(X, a, b, c):
    x, y, z = X
    dx = a*(y-x)
    dy = -c*y - x*z
    dz = -b + x*y
    return [dx,dy,dz]

def xu(X, a, b, c, d, e, f, g):
    x, y, z = X
    dx = g*z
    dy = d*x**2 + e*y**2 - f
    dz = -a*x - b*x**2 + c*y**2
    return [dx,dy,dz] 

def rossler(X, a, b, c):
    x, y, z = X
    dx = -y - z
    dy = x + a*y
    dz = b + z*(x - c)
    return [dx, dy, dz]

def chen(X, a, b, c):
    x, y, z = X
    dx = a*(y - x)
    dy = (c - a)*x - x*z + c*y
    dz = x*y - b*z
    return [dx, dy, dz]
