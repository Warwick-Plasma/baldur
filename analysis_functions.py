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
    self.all_time_data = kwargs.get('all_time_data', 1)
    self.grid =  kwargs.get('grid', 1)
    self.units = kwargs.get('units', 'Not set')
    self.name = kwargs.get('name', 'Not set')
    self.unit_conversion = kwargs.get('unit_conversion', 1)
    self.units_new = kwargs.get('units_new', 'Not set')



def basic(dat):
	
  fac = 1.0
  if dat.Logical_flags.use_rz:
    fac = 2.0 * np.pi
  
  if dat.Logical_flags.use_rz:
    v1 = dat.Velocity_VTheta.data
    v2 = dat.Velocity_Vr.data
    v3 = dat.Velocity_Vz.data
  else:
    v1 = dat.Velocity_Vx.data
    v2 = dat.Velocity_Vy.data
    v3 = dat.Velocity_Vz.data

  grid_mid = dat.Grid_Grid_mid.data
  xc = grid_mid[0]
  yc = grid_mid[1]
  grid = dat.Grid_Grid.data
  x = grid[0]
  y = grid[1]

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
                                    
  var_name = "Area"
  var_list.append(var_name)
  len_a = np.sqrt((x[:-1,:-1] - x[:-1,1:])**2 + (y[:-1,:-1] - y[:-1,1:])**2)
  len_b = np.sqrt((x[:-1,1:] - x[1:,1:])**2 + (y[:-1,1:] - y[1:,1:])**2)
  len_c = np.sqrt((x[1:,1:] - x[1:,:-1])**2 + (y[1:,1:] - y[1:,:-1])**2)
  len_d = np.sqrt((x[1:,:-1] - x[:-1,:-1])**2 + (y[1:,:-1] - y[:-1,:-1])**2)
  angle_a = np.arctan2(np.abs(y[:-1,:-1] - y[1:,:-1]), np.abs(x[:-1,:-1] - x[1:,:-1])) + np.arctan2(np.abs(y[:-1,1:] - y[:-1,:-1]),np.abs(x[:-1,1:] - x[:-1,:-1]))
  angle_c = np.arctan2(np.abs(y[:-1,1:] - y[1:,1:]),np.abs(x[:-1,1:] - x[1:,1:])) + np.arctan2(np.abs(y[1:,1:] - y[1:,:-1]),np.abs(x[1:,1:] - x[1:,:-1]))
  area = 0.5 * len_a * len_d * np.sin(angle_a) + 0.5 * len_b * len_c * np.sin(angle_c)
  setattr(dat, var_name, new_variable(data = area,
                                      grid = dat.Grid_Grid,
                                      units_new = "m$^2$",
                                      unit_conversion = 1,
                                      name = "Area"))
  
  var_name = "Cell_Mass"
  var_list.append(var_name)
  mass = dat.Fluid_Rho.data[:,:] * vol[:,:]
  setattr(dat, var_name, new_variable(data = mass,
                                      grid = dat.Grid_Grid,
                                      units_new = "kg/m$^3$",
                                      unit_conversion = 1,
                                      name = "Mass"))	
  var_name = "Fluid_Speed"
  var_list.append(var_name)
  speed = np.sqrt(v1**2 + v2**2 + v3**2)
  setattr(dat, var_name, new_variable(data = speed,
                                      grid = dat.Grid_Grid,
                                      units_new = "m/s$^1$",
                                      unit_conversion = 1,
                                      name = "Speed of Cell"))
  var_name = "Rho_r"
  var_list.append(var_name)
  dr = np.zeros(np.shape(radius))
  for i in range(len(radius)-1):
      dr[i, :] = radius[i+1] - radius[i]
  rho = dat.Fluid_Rho.data[:,:]
  rhor = rho * dr
  rhor_cumulative = np.cumsum(rhor, axis=0)
  setattr(dat, var_name, new_variable(data = rhor_cumulative,
                                      grid = dat.Grid_Grid,
                                      units_new = "g/cm$^2$",
                                      unit_conversion = 0.1,
                                      name = "Areal Density"))

  setattr(dat, "variables", var_list)
  
  # variables that only change in time
  var_list = dat.variables_time
  
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
  
  var_name = "Laser_Energy_Total_Deposited"
  var_list.append(var_name)
  tot_laser_dep = np.sum(np.sum(dat.Cell_Mass.data * laser_dep))
  setattr(dat, var_name, new_variable(data = tot_laser_dep,
                                      units_new = 'J',
                                      unit_conversion = 1,
                                      name = "Total laser energy deposited"))
  
  setattr(dat, "variables_time", var_list)
  
  return dat



def adiabat(dat, *args, **kwargs):
  call_basic = kwargs.get('call_basic', True)
  
  # Volume must be times by 2*pi in RZ
  fac = 1.0
  if dat.Logical_flags.use_rz:
    fac = 2*np.pi

  small_number = 1e-100

  rho = dat.Fluid_Rho.data
  pressure = dat.Fluid_Pressure.data
  
  # Variables that change in time and space
  var_list = dat.variables

  var_name = "Fluid_Adiabat"
  var_list.append(var_name)
  # The conversion to electron degeneracy pressure is only true for DT
  deg_pressure = 2.17e12 * (rho / 1000)**(5.0/3.0) / 10 + small_number
  adiabat = pressure / deg_pressure
  max_val = 100.0
  adiabat = np.where(adiabat < max_val, adiabat, max_val)
  setattr(dat, var_name, new_variable(data = adiabat,
                                      grid = dat.Grid_Grid,
                                      units_new = "unitless",
                                      unit_conversion = 1,
                                      name = "Adiabat"))

  var_name = "Fluid_Inverse_Pressure_Length_Scale"
  var_list.append(var_name)
  # As used by Craxton et al 2015 the inverse pressure scale length
  # makes the discontinous shock clear. Requires similar spatial and
  # temporal resolution
  dx = dat.Radius_mid.data[0][:-1,:-1] - dat.Radius_mid.data[0][1:,:-1]
  dlnp = np.log(pressure[:-1,:-1] + small_number) - np.log(pressure[1:,:-1] + small_number)
  pressure_ls = np.zeros(dat.Radius_mid.data[0].shape)
  pressure_ls[1:,1:] = np.abs(dlnp / dx)
  setattr(dat, var_name, new_variable(data = pressure_ls,
                                      grid = dat.Grid_Grid,
                                      units_new = "unitless",
                                      unit_conversion = 1,
                                      name = "Inverse Pressure Length Scale"))
  
  setattr(dat, "variables", var_list)

  return dat



def energy(dat, *args, **kwargs):
  call_basic = kwargs.get('call_basic', True)

  # Volume must be times by 2*pi in RZ
  fac = 1.0
  if dat.Logical_flags.use_rz:
    fac = 2*np.pi

  mass = dat.Cell_Mass.data
  corner_mass = dat.Test_Corner_Mass.data * fac
  vel_sqr = dat.Fluid_Speed.data**2
  Ei = dat.Fluid_Energy_ion.data * mass
  Ee = dat.Fluid_Energy_electron.data * mass

  # Variables that change in time and space
  var_list = dat.variables

  var_name = "Fluid_Energy_Kinetic"
  var_list.append(var_name)
  KE = 0.5 * (corner_mass[::2,::2] * vel_sqr[:-1,:-1]\
     + corner_mass[1::2,1::2] * vel_sqr[1:,1:]\
     + corner_mass[1::2,::2] * vel_sqr[1:,:-1]\
     + corner_mass[::2,1::2] * vel_sqr[:-1,1:])
  setattr(dat, var_name, new_variable(data = KE,
                                      grid = dat.Grid_Grid,
                                      units_new = "J",
                                      unit_conversion = 1,
                                      name = "Kinetic Energy"))

  setattr(dat, "variables", var_list)
  
  # variables that only change in time
  var_list = dat.variables_time
  
  var_name = "Total_Kinetic_Energy"
  var_list.append(var_name)
  tot_KE = np.sum(np.sum(KE))
  setattr(dat, var_name, new_variable(data = tot_KE,
                                      units_new = 'J',
                                      unit_conversion = 1,
                                      name = "Total Kinetic Energy"))
                                      
  var_name = "Total_Internal_Energy"
  var_list.append(var_name)
  tot_Ei = np.sum(np.sum(Ei))
  tot_Ee = np.sum(np.sum(Ee))
  tot_IE = tot_Ei + tot_Ee
  setattr(dat, var_name, new_variable(data = tot_IE,
                                      units_new = 'J',
                                      unit_conversion = 1,
                                      name = "Total Internal Energy"))
  
  var_name = "Total_Energy"
  var_list.append(var_name)
  tot_LE = dat.Laser_Energy_Total_Deposited.data
  tot_energy = tot_IE + tot_KE - tot_LE
  setattr(dat, var_name, new_variable(data = tot_energy,
                                      units_new = 'J',
                                      unit_conversion = 1,
                                      name = "Total Energy"))
  
  setattr(dat, "variables_time", var_list)
  
  return dat



def time_variables(dat, *args, **kwargs):
  """
  """

  # variables that only change in time
  var_list = dat.variables_time
  
  var_name = "Laser_Power_Total_Deposited"
  var_list.append(var_name)
  dt = dat.Times.all_time_data[1:] - dat.Times.all_time_data[:-1]
  tot_laser_pwr_dep = (dat.Laser_Energy_Total_Deposited.all_time_data[1:]
      - dat.Laser_Energy_Total_Deposited.all_time_data[:-1]) / dt
  tot_laser_pwr_dep = np.insert(tot_laser_pwr_dep, 0, 0)
  setattr(dat, var_name, new_variable(data = 0.0,
                                      all_time_data = tot_laser_pwr_dep,
                                      units_new = 'W',
                                      unit_conversion = 1,
                                      name = "Total laser power deposited"))
  
  setattr(dat, "variables_time", var_list)



def main():
	"""
	"""


if __name__ == "__main__":
	main()
