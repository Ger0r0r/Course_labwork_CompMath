import pandas as pd
import numpy as np
import sympy as sp
import shutil
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from datetime import datetime
import time
from PIL import Image

def create_image(Data, Y, X):
    raw_image = np.full((Y,X,4), 255)
    H = np.max(Data)
    L = np.min(Data)
    
    raw_image[:,:,1] = (Data[:,:]-L)*255/(H-L)
    raw_image[:,:,0] = 0
    raw_image[:,:,2] = 0
    raw_image[:,:,3] = 255
    return raw_image

H = 1
l = 0.3
L = 0.7
X = 3
T = 1
Kurant = 0.1


def calc_sigma (ts, xs, u):
	return ts/xs*np.max(u)

def calc_initial_t (t):
	if t <l or t >=L:
		return 0
	if t >= l and t < (l+L)/2:
		return (2*(H)/(L-l))*(t-l)
	if t >= (l+L)/2 and t < L:
		return -(2*(H)/(L-l))*(t-L)

def calc_initial_x (x):
	# return 0
	if x < l or x >= L:
		return 0
	if x >= l and x < L:
		return H

def left_angel(ud, udl, xs, ts):
	return ts * (ud**2 - udl**2)/(2*xs) + ud


def scheme_L(udl, udr, xs, ts):
    return 0.5*(udr + udl) - 0.5*(udr**2-udl**2)*ts/(xs)


def create_GIF(ts, Kurant):

	dir_name = "./temp/"
	os.mkdir(dir_name)

	xs = ts * int(1 / Kurant)

	Nx = int(X/xs)+1
	Nframes = int(T/ts)+1

	# CALCULATION
	u = np.full((Nframes,Nx), 0, dtype="float64")
	u[0,:] = np.linspace(0,X,Nx)
	u[:,0] = np.linspace(0,T,Nframes)
	for i in range (0, Nx):
		u[0,i] = calc_initial_x(u[0,i])
	for i in range (0, Nframes):
		u[i,0] = calc_initial_t(u[i,0])

	L = u.copy()

	for i in range(1, Nframes):
		for j in range (1,Nx):
			if (j != Nx-1):
				L[i,j] = scheme_L(L[i-1,j-1],L[i-1,j+1],xs,ts)
			else:
				L[i,j] = scheme_L(L[i-1,j-1],L[i-1,j],xs,ts)

	# MAKE FRAMES
	x = np.linspace(0, X, Nx)

	for frame in range(0, Nframes):
		fig = plt.figure(figsize = [12,7])
		fig = plt.plot(x, L[frame, :], '-b', ms = 0.1)
		fig = plt.xlabel('x')
		fig = plt.xlim([0, X])
		fig = plt.ylabel('u(x)')
		fig = plt.title(f"t = {round(frame * ts, 4)}")
		fig = plt.savefig(dir_name + f"{frame}.png")

	# GIF CREATION
	im = []
	for frame in range(0, Nframes):
		img = Image.open(dir_name + f"{frame}.png")
		img.load()
		im.append(img) 

	im[0].save(f"./GIF/L_{xs}.gif", save_all = "True", append_images = im[1:], duration = 100)

	shutil.rmtree(dir_name)



def do_Laks(ts, Kurant):

	start_time = datetime.now()

	xs = ts * int(1 / Kurant)

	print("L prec: ", xs)
	Nx = int(X/xs)+1
	Nt = int(T/ts)+1

	u = np.full((Nt,Nx), 0, dtype="float64")
	u[0,:] = np.linspace(0,X,Nx)
	u[:,0] = np.linspace(0,T,Nt)
	for i in range (0, Nx):
		u[0,i] = calc_initial_x(u[0,i])
	for i in range (0, Nt):
		u[i,0] = calc_initial_t(u[i,0])
	# print(u)

	L = u.copy()

	for i in range (1, Nt):
		for j in range (1,Nx):
			if (j != Nx-1):
				L[i,j] = scheme_L(L[i-1,j-1],L[i-1,j+1],xs,ts)
			else:
				L[i,j] = scheme_L(L[i-1,j-1],L[i-1,j],xs,ts)
	# print(L)

	calc_time = datetime.now() - start_time
	print("calc: ", calc_time)

	df = pd.DataFrame(data=L)
	df.to_csv(f"./tables/L-{xs}.csv")
	# ./Course_labwork_CompMath/tables/

	save_time = datetime.now() - start_time - calc_time
	print("save: ", save_time)

	x = np.linspace(0, X, Nx)
	t = np.linspace(0, T, Nt)

	raw = create_image(L,Nt,Nx)
	im = Image.fromarray(raw.astype(np.uint8))
	size = (1000,1000)
	im = im.resize(size)
	im = im.transpose(Image.Transpose.FLIP_TOP_BOTTOM)
	im.save(f"./test/L-{xs}.png")

	total_time = datetime.now() - start_time - save_time - calc_time
	print("show: ", total_time)

create_GIF(0.001, 0.1)

# t = [0.01, 0.001, 0.0001]
# Kurant = 0.2

# for ts in t:
# 	do_Laks(ts, Kurant)