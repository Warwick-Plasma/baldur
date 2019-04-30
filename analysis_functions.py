# Duncan Barlow, Odin project, Warwick University, 01/19
import sdf_helper as sh
import matplotlib.pyplot as plt
import numpy as np
import glob
import sys, os
from matplotlib.widgets import Slider, RadioButtons



class new_variable:
  """
  """
  def __init__(self, *args, **kwargs):
    self.data = 1 #np.zeros(RunCounter)
    self.label = 'None'
    self.grid = 1



def basic(dat):
	
	fac = 1.0
	if dat.Logical_flags.use_rz:
		fac = 2.0 * np.pi

	grid_mid = dat.Grid_Grid_mid.data
	xc = grid_mid[0]
	yc = grid_mid[1]
	
	var_name = "Radius_mid"
	setattr(dat, var_name, dat.Grid_Grid_mid)
	radius = np.sqrt(xc**2 + yc**2)
	dat.Radius_mid.data = radius
	
	var_list = ["None"]
	
	var_name = "Fluid_Volume_rz"
	var_list.append(var_name)
	setattr(dat, var_name, new_variable())
	vol = dat.Fluid_Volume.data * fac
	dat.Fluid_Volume_rz.data = vol
	dat.Fluid_Volume_rz.grid = grid_mid
	
	var_name = "Fluid_Mass"
	var_list.append(var_name)
	setattr(dat, var_name, new_variable())
	mass = dat.Fluid_Rho.data[:,:] * vol[:,:]
	dat.Fluid_Mass.data = mass
	dat.Fluid_Mass.grid = grid_mid
	
	var_name = "Centre_Of_Mass"
	var_list.append(var_name)
	setattr(dat, var_name, new_variable())
	dat.Centre_Of_Mass.data = np.sum(np.sum(mass * radius)) / np.sum(np.sum(mass))
	
	var_list = var_list[1:]
	setattr(dat, "variables_time", var_list)
	
	return dat



def laser(dat):
	
	dat = basic(dat)
	
	var_list = dat.variables_time
	
	var_name = "Total_Energy_Laser"
	var_list.append(var_name)
	setattr(dat, var_name, new_variable())
	laser_dep = dat.Fluid_Energy_deposited_laser.data
	tot_laser_dep = np.sum(np.sum(dat.Fluid_Mass.data * laser_dep))
	dat.Total_Energy_Laser.data = tot_laser_dep
	
	
	return dat



def main():
	"""
	"""


if __name__ == "__main__":
	main()
