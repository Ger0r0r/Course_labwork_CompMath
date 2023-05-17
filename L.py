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

H = 1
l = 0.3
L = 0.7
X = 3
T = 1
Kurant = 0.1

def f(u):
	return (u**2)/2

def create_image(Data, Y, X):
    raw_image = np.full((Y,X,4), 255)
    H = np.max(Data)
    L = np.min(Data)
    
    raw_image[:,:,1] = (Data[:,:]-L)*255/(H-L)
    raw_image[:,:,0] = 0
    raw_image[:,:,2] = 0
    raw_image[:,:,3] = 255
    return raw_image


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
	if x < l or x >= L:
		return 0
	if x >= l and x < L:
		return H

def left_angel(ud, udl, xs, ts):
	return ts* (ud**2 - udl**2)/(xs*2) + ud


def scheme_L(udl, udr, Kurant):
    return 0.5*(udr + udl) - 0.5*(f(udr)-f(udl))*Kurant

def scheme_LF(ud, udl, udr, Kurant):
	uar = 0.5*(udr+ud) - (udr**2-ud**2)*Kurant/(4) 
	ual = 0.5*(ud+udl) - (ud**2-udl**2)*Kurant/(4)
	return ud - (uar**2-ual**2)*Kurant/(2)

def scheme_CIR(ud, udl, udr, Kurant):
    if ud < 0:
        return ud - (udr**2 - ud**2)*Kurant/(2)
    else:
        return ud - (ud**2 - udl**2)*Kurant/(2)

def calc_t(u,ur,s):
    return u-2*s*(f(ur)-f(u))/3

def calc_tt(u,ut,utl,s):
    return (u+ut)/2-s*(f(ut)-f(utl))/3

def scheme_WCL(ud, udl, udll, udr, udrr, xs, ts, Kurant, omega = 1.66):
    # omega = 4*Kurant**2-Kurant**4
    
    utll = calc_t(udll,udl,Kurant)
    utl = calc_t(udl,ud,Kurant)
    ut = calc_t(ud,udr,Kurant)
    utr = calc_t(udr,udrr,Kurant)
    
    uttl = calc_tt(udl,utl,utll,Kurant)
    uttr = calc_tt(udr,utr,ut,Kurant)
    
    F = 2*f(udrr)- 7*f(udr) + 7*f(udl) - 2*f(udll)
    G = f(uttr)-f(uttl)
    H = udrr-4*udr+6*ud-4*udl+udll
    return ud - Kurant*2*(F/24+3*G/8+omega*H/24)

def calculation (xs, ts, Kurant, scheme = "L", omega = 1.66):

	Nx = int(X/xs)+1
	Nt = int(T/ts)+1

	u = np.full((Nt,Nx), 0, dtype="float128")
	u[0,:] = np.linspace(0,X,Nx)
	u[:,0] = np.linspace(0,T,Nt)
	for i in range (0, Nx):
		u[0,i] = calc_initial_x(u[0,i])
	for i in range (0, Nt):
		u[i,0] = calc_initial_t(u[i,0])

	L = u.copy()

	if (scheme == "L"):
		for i in range(1, Nt):
			for j in range (1,Nx):
				if (j != Nx-1):
					L[i,j] = scheme_L(L[i-1,j-1],L[i-1,j+1], Kurant)
				else:
					L[i,j] = scheme_L(L[i-1,j-1],L[i-1,j], Kurant)
	
	if (scheme == "LF"):
		for i in range (1, Nt):
			for j in range (1,Nx):
				if (j != Nx-1):
					L[i,j] = scheme_LF(L[i-1,j],L[i-1,j-1],L[i-1,j+1], Kurant)
				else:
					# experiment
					L[i, j] = scheme_LF(L[i-1,j],L[i-1,j-1],L[i-1,j], Kurant)
					# L[i,j] = left_angel(L[i-1,j],L[i-1,j-1],xs,ts)

	if (scheme == "CIR"):
		for i in range (1, Nt):
			for j in range (1,Nx):
				if (j != Nx-1):
					L[i,j] = scheme_CIR(L[i-1,j],L[i-1,j-1], L[i-1,j+1],Kurant)
				else:
				    # experiment
					L[i,j] = scheme_CIR(L[i-1,j], L[i-1,j-1], L[i-1, j],Kurant)

	if (scheme == "WCL"):
		for i in range (1, Nt):
			for j in range (1,Nx):
				if ((j != Nx-1) and (j != Nx-2) and (j != 1)):
					L[i,j] = scheme_WCL(L[i-1,j], L[i-1,j-1],L[i-1,j-2],L[i-1,j+1],L[i-1,j+2],xs,ts,ts/xs, omega)
				else:
					# experiment
					L[i,j] = scheme_WCL(L[i-1,j], L[i-1,j-1],L[i-1,j-2],L[i-1,j],L[i-1,j],xs,ts,ts/xs, omega)
					# L[i,j] = left_angel(L[i-1,j],L[i-1,j-1],xs,ts)
	return L



def create_GIF(xs, Kurant, scheme = "L"):

	dir_name = "./temp/"
	os.mkdir(dir_name)

	ts = xs * Kurant

	Nx = int(X/xs)+1
	Nframes = int(T/ts)+1

	if (Nframes > 100):
		stepes = int(Nframes / 100)
	else:
		stepes = 1

	L = calculation(xs, ts, Kurant, scheme)

	# MAKE FRAMES
	x = np.linspace(0, X, Nx)

	for frame in range(0, Nframes,stepes):
		fig = plt.figure(figsize = [12,7])
		fig = plt.plot(x, L[frame, :], '-', ms = 0.1)
		fig = plt.xlabel('x')
		fig = plt.xlim([0, X])
		fig = plt.ylim ([0, 1.5])
		fig = plt.ylabel('u(x)')
		fig = plt.title(f"t = {round(frame * ts, 4)}")
		fig = plt.savefig(dir_name + f"{frame}.png")

	# GIF CREATION
	im = []
	for frame in range(0, Nframes, stepes):
		img = Image.open(dir_name + f"{frame}.png")
		img.load()
		im.append(img) 

	im[0].save(f"./GIF/" + scheme + f"_{xs}_{Kurant}.gif", save_all = "True", append_images = im[1:], duration = 60)

	shutil.rmtree(dir_name)



def do_scheme(ts, Kurant, scheme = "L"):

	start_time = datetime.now()

	xs = ts * int(1 / Kurant)

	print("X prec: ", xs)
	Nx = int(X/xs)+1
	Nt = int(T/ts)+1

	L = calculation (ts, Kurant, scheme)

	calc_time = datetime.now() - start_time
	print("calc: ", calc_time)

	df = pd.DataFrame(data=L)
	df.to_csv(f"./tables/L-{xs}.csv")

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


x = [0.005]
schemes = ["L", "LF", "CIR", "WCL"]
Kurant = [0.7]

for xs in x:
	for c in Kurant:
		for sc in schemes:
			create_GIF(xs, c, sc)
