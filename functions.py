def lorenz_equations(X, sigma, rho, beta):
    x, y, z = X
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    return [dx, dy, dz]

def moore_spiegel(X, T, R):
    #Moore-Spiegel
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

