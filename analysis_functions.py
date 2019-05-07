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
	
	# Grids, all grids need to be 3D arrays so we stack radius
	var_list = dat.grids
	
	var_name = "Radius_mid"
	var_list.append(var_name)
	radius = np.sqrt(xc**2 + yc**2)
	theta = np.arctan2(yc, xc)
	setattr(dat, var_name, new_variable(data = np.array([radius, theta]),
	                                    units_new = dat.Grid_Grid_mid.units_new,
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
	
	var_name = "Times"
	var_list.append(var_name)
	setattr(dat, var_name, new_variable(data = dat.Header["time"],
	                                    units_new = "ns",
	                                    unit_conversion = 1.0e9,
	                                    name = "Time"))
	
	var_name = "Centre_Of_Mass"
	var_list.append(var_name)
	com = np.sum(np.sum(mass * radius)) / np.sum(np.sum(mass))
	setattr(dat, var_name, new_variable(data = com,
	                                    units_new = dat.Grid_Grid_mid.units_new,
	                                    unit_conversion = dat.Grid_Grid_mid.unit_conversion,
	                                    name = "Centre of Mass"))
	
	setattr(dat, "variables_time", var_list)
	
	return dat



def laser(dat, *args, **kwargs):
  call_basic = kwargs.get('call_basic', True)
  laser_change = kwargs.get('laser_change', False)
  sdf_num = kwargs.get('sdf_num', 0)
  istart = kwargs.get('istart', 0)
  pathname = kwargs.get('pathname', os.path.abspath(os.getcwd()))
  
  if call_basic:
    dat = basic(dat)
  
  laser_dep = dat.Fluid_Energy_deposited_laser.data
  
  if laser_change:
    # Variables that change in time and space
    var_list = dat.variables
	  
    var_name = "Laser_Energy_per_step"
    var_list.append(var_name)
    if sdf_num == istart:
      lap_dep_step = laser_dep
    elif sdf_num >= istart:
      SDFName=pathname+'/'+str(sdf_num-1).zfill(4)+'.sdf'
      dat2 = sh.getdata(SDFName,verbose=False)
      lap_dep_step = laser_dep - dat2.Fluid_Energy_deposited_laser.data
    else:
      print('Error with laser change calculation')
      print('sdf_num = ', sdf_num, ' and the minimum = ', istart)
    setattr(dat, var_name, new_variable(data = lap_dep_step,
                                        grid = dat.Grid_Grid,
                                        units_new = "J/kg",
                                        unit_conversion = 1,
                                        name = "Laser Energy Deposited"))
  
    setattr(dat, "variables", var_list)
  
  # variables that only change in time
  var_list = dat.variables_time
  
  var_name = "Total_Energy_Laser_deposited"
  var_list.append(var_name)
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
