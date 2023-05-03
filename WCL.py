import pandas as pd
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from matplotlib import cm

H = 1
l = 0.3
L = 0.7
X = 1
T = 1

def calc_sigma (ts, xs, u):
	return ts/xs*np.max(u)

def calc_initial_t (t):
    # return 1
	if t < l or t >= L:
		return 0
	if t >= l and t < (l+L)/2:
		return (2*H/(L-l))*(t-l)
	if t >= (l+L)/2 and t < L:
		return -(2*H/(L-l))*(t-L)

def calc_initial_x (x):
    # return 2*x - x**2 + 1
	if x < l or x >= L:
		return 0
	if x >= l and x < L:
		return H

def left_angel(ud, udl, xs, ts):
    return ud - ts*ud*(ud-udl)/xs

def calc_t(u,ur,s):
    return u-2*s*(ur**2/2-u**2/2)/3

def calc_tt(u,ut,utl,s):
    return (u+ut)/2-s*(ut**2/2-utl**2/2)/3

def scheme_WCL(ud, udl, udll, udr, udrr, xs, ts, s):
    # omega = (4*s**2+1)*(4-s**2)/5
    omega = 4*s**2-s**4
    # omega = 3
    
    utll = calc_t(udll,udl,s)
    utl = calc_t(udl,ud,s)
    ut = calc_t(ud,udr,s)
    utr = calc_t(udr,udrr,s)
    
    uttl = calc_tt(udl,utl,utll,s)
    uttr = calc_tt(udr,utr,ut,s)
    
    F = -2*udrr**2/2+7*udr**2/2-7*udl**2/2+2*udll**2/2
    G = uttr**2/2-uttl**2/2
    H = udrr-4*udr+6*ud-4*udl+6*udll
    return ud - s*(F/24+3*G/8+omega*H/24)

ts = 0.006
xs = 0.06

Nx = int(X/xs)+1
Nt = int(T/ts)+1

# Инициализация начальных данных сетки
u = np.full((Nt,Nx), 0, dtype="float64")
u[0,:] = np.linspace(0,X,Nx)
u[:,0] = np.linspace(0,T,Nt)
for i in range (0, Nx):
    u[0,i] = calc_initial_x(u[0,i])
for i in range (0, Nt):
    u[i,0] = calc_initial_t(u[i,0])
print(u)

WCL = u.copy()

for i in range (1, Nt):
    for j in range (1,Nx):
        if ((j != Nx-1) and (j != Nx-2) and (j != 1)):
            WCL[i,j] = scheme_WCL(WCL[i-1,j],WCL[i-1,j-1],WCL[i-1,j-2],WCL[i-1,j+1],WCL[i-1,j+2],xs,ts,ts/xs)
        else:
            WCL[i,j] = left_angel(WCL[i-1,j],WCL[i-1,j-1],xs,ts)
   
print(WCL)

# Сохранение данных в таблицу excel
df = pd.DataFrame(data=WCL)
df.to_excel("WCL.xlsx")

# Создание 3D графика
x = np.linspace(0, X, Nx)
t = np.linspace(0, T, Nt)

xx, yy = np.meshgrid(x, t)

fig, ax = plt.subplots(subplot_kw = {"projection" : "3d" })

surf = ax.plot_surface(xx, yy, np.array(WCL), cmap = cm.turbo, linewidth = 0, antialiased = True)

fig.colorbar(surf, shrink = 0.5, aspect = 5)

plt.show()