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
    self.data = kwargs.get('data', 1)
    self.grid =  kwargs.get('grid', 1)
    self.units = kwargs.get('units', 'Not set')
    self.name = kwargs.get('name', 'Not set')
    self.unit_conversion = kwargs.get('unit_conversion', 1)
    self.units_new = kwargs.get('units_new', 'Not set')



def basic(dat):
	
	fac = 1.0
	if dat.Logical_flags.use_rz:
		fac = 2.0 * np.pi

	grid_mid = dat.Grid_Grid_mid.data
	xc = grid_mid[0]
	yc = grid_mid[1]
	
	# Grids
	var_list = dat.grids
	
	var_name = "Radius_mid"
	var_list.append(var_name)
	radius = np.sqrt(xc**2 + yc**2)
	setattr(dat, var_name, new_variable(data = radius,
	                                    units = dat.Grid_Grid_mid.units_new,
	                                    unit_conversion = dat.Grid_Grid_mid.unit_conversion,
	                                    name = "Radius"))
	
	setattr(dat, "grids", var_list)
	
	# Variables that change in time and space
	var_list = dat.variables
	
	var_name = "Fluid_Volume_rz"
	var_list.append(var_name)
	vol = dat.Fluid_Volume.data * fac
	setattr(dat, var_name, new_variable(data = vol,
	                                    grid = dat.Grid_Grid,
	                                    units_new = "m$^3$",
	                                    unit_conversion = 1,
	                                    name = "Volume"))
	
	var_name = "Cell_Mass"
	var_list.append(var_name)
	mass = dat.Fluid_Rho.data[:,:] * vol[:,:]
	setattr(dat, var_name, new_variable(data = mass,
	                                    grid = dat.Grid_Grid,
	                                    units_new = "kg/m$^3$",
	                                    unit_conversion = 1,
	                                    name = "Mass"))
  
	setattr(dat, "variables", var_list)
	
	# variables that only change in time
	var_list = dat.variables_time
	
	var_name = "Centre_Of_Mass"
	var_list.append(var_name)
	com = np.sum(np.sum(mass * radius)) / np.sum(np.sum(mass))
	setattr(dat, var_name, new_variable(data = com,
	                                    units = dat.Grid_Grid_mid.units_new,
	                                    unit_conversion = dat.Grid_Grid_mid.unit_conversion,
	                                    name = "Centre of Mass"))
	
	setattr(dat, "variables_time", var_list)
	
	return dat



def laser(dat, *args, **kwargs):
  call_basic = kwargs.get('call_basic', True)
  
  if call_basic:
    dat = basic(dat)
  
  # variables that only change in time
  var_list = dat.variables_time
  
  var_name = "Total_Energy_Laser_deposited"
  var_list.append(var_name)
  laser_dep = dat.Fluid_Energy_deposited_laser.data
  tot_laser_dep = np.sum(np.sum(dat.Cell_Mass.data * laser_dep))
  setattr(dat, var_name, new_variable(data = tot_laser_dep,
                                      units_new = 'J',
                                      unit_conversion = 1,
                                      name = "Total laser energy deposited"))
  
  setattr(dat, "variables_time", var_list)
  
  return dat



def main():
	"""
	"""


if __name__ == "__main__":
	main()
