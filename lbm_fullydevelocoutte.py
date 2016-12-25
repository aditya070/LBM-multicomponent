''' not acheived the objective'''
weights = [4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36]
cx = [0, 1, 0, -1, 0, 1, -1, -1, 1]
cy = [0, 0, 1, 0, -1, 1, 1, -1, -1]
q = 9
nx = 20
ny = 20
u0 = 0.05
omega = 0.56
ts = 2
rho0 = 1
feq = [0] * q
f = [[[0] * q] * nx] * ny
fs = [[[0] * q] * nx] * ny
ux = [[[0] * q] * nx] * ny
uy = [[[0] * q] * nx] * ny
rho = [[[0] * q] * nx] * ny

# initialise the distributive function
for x in range(0, nx):
    for y in range(0, ny):
        for k in range(0, q):
            feq[k] = weights[k] * rho0
            f[x][y][k] = feq[k]

# main for ts time step

for t in range(1, ts):
    for x in range(0, nx):
        for y in range(0, ny):
            dense = 0
            vx = 0
            vy = 0
    for k in range(0, q):
        dense = dense + f[x][y][k]
        vx = vx + cx[k] * f[x][y][k]
        vy = vy + cy[k] * f[x][y][k]

    if not dense==0:
        rho[x][y] = dense
        ux[x][y] = vx / dense
        uy[x][y] = vy / dense
    else:
        rho[x][y] = 0
        ux[x][y] = 0
        uy[x][y] = 0

    for k in range(0, q):
        feq[k] = weights[k] * rho[x][y] * (
        1 + 3 * (cx[k] * ux[x][y] + cy[k] * uy[x][y]) + 4.5 * (cx[k] * ux[x][y] + cy[k] * uy[x][y]) ** 2 - 1.5 * (ux[x][y] ** 2 + uy[x][y] ** 2))
        f[x][y][k] = f[x][y][k] * (1 - omega) + feq[k] * omega
        newx = (x + cx[k]) % nx
        newy = (y + cy[k]) % ny
        fs[newx][newy][k] = f[x][y][k]

    f = fs

for x in range(0, nx):
    for y in range(0, ny):
        if y == ny:
            rhow = f[0][x][y] + f[1][x][y] + f[3][x][y]+2 * (f[6][x][y] + f[2][x][y] + f[5][x][y])
            # uy[x][y]= ((f[0][x][y]+f[1][x][y]+f[3][x][y]+2*(f[2][x][y]+f[5][x][y]+f[6][x][y]))/rho0)-1
            # f[4][x][y]= f[2][x][y]-(2/3)*rho0*uy[x][y]
            f[4][x][y] = f[2][x][y]
            # f[7][x][y]= f[5][x][y]+0.5*(f[1][x][y]-f[3][x][y])-(1/6)*rho0*uy[x][y]-0.5*rho0*ux[x][y]
            f[7][x][y] = f[5][x][y] - 0.5 * rhow * u0 + 0.5 * (f[1][x][y] - f[3][x][y])
            # f[8][x][y] =f[6][x][y]-0.5*(f[1][x][y]-f[3][x][y])-(1/6)*rho0*uy[x][y]+0.5*rho0*ux[x][y]
            f[8][x][y] = f[6][x][y] + 0.5 * rhow * u0 - 0.5 * (f[1][x][y] - f[3][x][y])

        if y == 1:
            f[2][x][y] = f[4][x][y]
            f[5][x][y] = f[7][x][y]
            f[6][x][y] = f[8][x][y]

for x in range(0, nx):
    for y in range(0, ny):
        for k in range(0, q):
            print(f[x][y][k])
        print("\n")

            #  i = 1:1:ny
            # plot(i,ux(15,i))
            # hold on
