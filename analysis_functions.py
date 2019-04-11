# Duncan Barlow, Odin project, Warwick University, 01/19
import sdf_helper as sh
import matplotlib.pyplot as plt
import numpy as np
import glob
import sys, os
from matplotlib.widgets import Slider, RadioButtons
		


def basic(dat, one_sdf):
	
	fac = 1.0
	if dat.Logical_flags.use_rz:
		fac = 2*np.pi
	
	vol = dat.Fluid_Volume.data * fac
	setattr(one_sdf, "Fluid_Volume", vol)
	
	mass = dat.Fluid_Rho.data[:,:] * vol[:,:]
	setattr(one_sdf, "Fluid_Mass", mass)
	
	radius = np.sqrt(dat.Grid_Grid_mid.data[0]**2 + dat.Grid_Grid_mid.data[1]**2)
	
	com = np.sum(np.sum(mass * one_sdf.radius)) / np.sum(np.sum(mass))
	setattr(one_sdf, "com", com)
	
	return one_sdf



def laser(dat, one_sdf):
	
	one_sdf = basic(dat, one_sdf)
	
	laser_dep = dat.Fluid_Energy_deposited_laser.data
	tot_laser_dep = np.sum(np.sum(one_sdf.Fluid_Mass * laser_dep))
	setattr(one_sdf, "tot_laser_dep", tot_laser_dep)
	
	return one_sdf



def main():
	"""
	"""


if __name__ == "__main__":
	main()
