import pandas as pd
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from matplotlib import cm

H = 1
l = 3
L = 7
X = 10
T = 10

def calc_sigma (ts, xs, u):
	return ts/xs*np.max(u)

def calc_initial_t (t):
	if t < l or t >= L:
		return 0
	if t >= l and t < (l+L)/2:
		return (2*H/(L-l))*(t-l)
	if t >= (l+L)/2 and t < L:
		return -(2*H/(L-l))*(t-L)

def calc_initial_x (x):
	if x < l or x >= L:
		return 0
	if x >= l and x < L:
		return H

def left_angel(ud, udl, xs, ts):
    return ts*(ud**2 - udl**2)/(4*xs) + ud

def sceme_L(ud, udl, udr, xs, ts):
    return (udl**2-udr**2)*ts/(4*xs) + 0.5*(udr + udl)

ts = 0.1
xs = 0.1

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

L = u.copy()

for i in range (1, Nt):
	for j in range (1,Nx):
		if (j != Nx-1):
			L[i,j] = sceme_L(L[i-1,j],L[i-1,j-1],L[i-1,j+1],xs,ts)
		else:
			L[i,j] = left_angel(L[i-1,j],L[i-1,j-1],xs,ts)
   
print(L)

# Сохранение данных в таблицу excel
df = pd.DataFrame(data=L)
df.to_excel("L.xlsx")

# Создание 3D графика
x = np.linspace(0, X, Nx)
t = np.linspace(0, T, Nt)

xx, yy = np.meshgrid(x, t)

fig, ax = plt.subplots(subplot_kw = {"projection" : "3d" })

surf = ax.plot_surface(xx, yy, np.array(L), cmap = cm.turbo, linewidth = 0, antialiased = True)

fig.colorbar(surf, shrink = 0.5, aspect = 5)

plt.show()