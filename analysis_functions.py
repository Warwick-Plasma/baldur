# Duncan Barlow, Odin project, Warwick University, 01/19
import sdf_helper as sh
import matplotlib.pyplot as plt
import numpy as np
import glob
import sys, os
from matplotlib.widgets import Slider, RadioButtons


global small_number
small_number = 1e-100


class new_variable:
  """
  """
  def __init__(self, *args, **kwargs):
    self.data = kwargs.get('data', 1)
    self.index = kwargs.get('indices', 1)
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
  
  nmat = dat.Integer_flags.nmat
  amu = 1.66053904e-27

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
  rho = dat.Fluid_Rho.data
  rhor = rho * dr
  rhor_cumulative = np.cumsum(rhor, axis=0)
  setattr(dat, var_name, new_variable(data = rhor_cumulative,
                                      grid = dat.Grid_Grid,
                                      units_new = "g/cm$^2$",
                                      unit_conversion = 0.1,
                                      name = "Areal Density"))
                                      
  var_name = "Fluid_Number_density_ion"
  var_list.append(var_name)
  ni_density = 0.0
  for imat in range(1,nmat+1):
    mat_name = getattr(getattr(dat, "material_string_flags_"+str(imat).zfill(3)), "data")['name']
    a_bar = getattr(getattr(dat, "material_real_flags_"+str(imat).zfill(3)), "a_bar")
    mat_den = getattr(getattr(dat, "Fluid_Rho_"+mat_name), "data")
    ni_density = ni_density + mat_den / a_bar / amu
    
  setattr(dat, var_name, new_variable(data = ni_density,
                                      grid = dat.Grid_Grid,
                                      units_new = "#/m$^3$",
                                      unit_conversion = 1,
                                      name = "Number density of ions"))
  
  var_name = "Fluid_Number_density_electron"
  var_list.append(var_name)
  Z = dat.Fluid_Charge_State.data
  ne_density = Z * ni_density
  setattr(dat, var_name, new_variable(data = ne_density,
                                      grid = dat.Grid_Grid,
                                      units_new = "#/m$^3$",
                                      unit_conversion = 1,
                                      name = "Number density of electrons"))

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
  
  radius = dat.Radius_mid.data[0]
  
  laser_wavelength = 351.0e-9
  laser_k = 2 * np.pi / laser_wavelength
  n_crit = 8.8e14 / laser_wavelength**2
  laser_dir = -1
  
  laser_dep = dat.Fluid_Energy_deposited_laser.data
  
  var_name = "Critical_Density"
  setattr(dat, var_name, new_variable(data = n_crit,
                                      units_new = "#/m^3",
                                      unit_conversion = 1,
                                      name = "Critical density"))

  ne_density = dat.Fluid_Number_density_electron.data
  crit_crossing = ne_density - n_crit
  quart_crit_crossing = ne_density - n_crit / 4.0
  nx, ny = ne_density.shape
  crit_surf_ind = [0] * ny
  crit_rad = np.zeros(ny)
  quart_crit_surf_ind = [0] * ny
  quart_crit_rad = np.zeros(ny)
  # 1 critical surface is chosen based on direction of laser propagation
  for iy in range(0,ny):
    zero_crossings = np.where(np.diff(np.sign(crit_crossing[:,iy])))[0]
    zero_crossings = np.append(0, zero_crossings)
    crit_surf_ind[iy] = int(zero_crossings[laser_dir])
    crit_rad[iy] = radius[crit_surf_ind[iy],iy]
    
    zero_crossings = np.where(np.diff(np.sign(quart_crit_crossing[:,iy])))[0]
    zero_crossings = np.append(0, zero_crossings)
    quart_crit_surf_ind[iy] = int(zero_crossings[laser_dir])
    quart_crit_rad[iy] = radius[quart_crit_surf_ind[iy],iy]
  
  var_list = dat.track_surfaces
  
  var_name = "Critical_Surface"
  var_list.append(var_name)
  setattr(dat, var_name, new_variable(data = crit_rad,
                                      index = crit_surf_ind,
                                      units_new = dat.Grid_Grid.units_new,
                                      unit_conversion = dat.Grid_Grid.unit_conversion,
                                      name = "Location of critical surface"))
    
  var_name = "Critical_Surface_quarter"
  var_list.append(var_name)
  setattr(dat, var_name, new_variable(data = quart_crit_rad,
                                      index = quart_crit_surf_ind,
                                      units_new = dat.Grid_Grid.units_new,
                                      unit_conversion = dat.Grid_Grid.unit_conversion,
                                      name = "Location of quarter critical surface"))
  
  setattr(dat, "track_surfaces", var_list)
  
  # Variables that change in time and space 
  var_list = dat.variables

  var_name = "Laser_Energy_per_step"
  var_list.append(var_name)
  if sdf_num == istart:
    las_dep_step = laser_dep
  elif sdf_num >= istart:
    SDFName=pathname+'/'+str(sdf_num-1).zfill(4)+'.sdf'
    dat2 = sh.getdata(SDFName,verbose=False)
    las_dep_step = laser_dep - dat2.Fluid_Energy_deposited_laser.data
  else:
    print('Error with laser change calculation')
    print('sdf_num = ', sdf_num, ' and the minimum = ', istart)
  setattr(dat, var_name, new_variable(data = las_dep_step,
                                      grid = dat.Grid_Grid,
                                      units_new = "J/kg",
                                      unit_conversion = 1,
                                      name = "Laser Energy Deposited"))
  
  var_name = "Fluid_Number_density_electron_per_critical"
  var_list.append(var_name)
  ne_per_crit = ne_density / n_crit
  setattr(dat, var_name, new_variable(data = ne_per_crit,
                                      grid = dat.Grid_Grid,
                                      units_new = '$n_{crit}$',
                                      unit_conversion = 1,
                                      name = "Electron number density"))
  
  var_name = "Fluid_Density_scale_length"
  var_list.append(var_name)
  grad_ne_density = gradient_function(ne_density, dat.Grid_Grid_mid.data)
  density_scale_length = np.zeros(dat.Grid_Grid_mid.data[0].shape)
  density_scale_length[1:-1,1:-1] = abs(ne_density[1:-1,1:-1] / (grad_ne_density[1:-1,1:-1] + small_number))
  # Remeber unit conversions are applied after!!
  max_val = 1e-2
  density_scale_length = np.where(density_scale_length < max_val, density_scale_length, 0.0)
  setattr(dat, var_name, new_variable(data = density_scale_length,
                                      grid = dat.Grid_Grid,
                                      units_new = dat.Grid_Grid.units_new,
                                      unit_conversion = dat.Grid_Grid.unit_conversion,
                                      name = "Density Scale length $l_n$"))
  
  
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
  pressure_ls = gradient_function(np.log(pressure + small_number), dat.Grid_Grid_mid.data)
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

def gradient_function(param, grid):
  
  xc = grid[0]
  yc = grid[1]
  
  dx = ((xc[:-2,1:-1] - xc[2:,1:-1])**2 +
        (yc[:-2,1:-1] - yc[2:,1:-1])**2)**0.5
  dy = ((xc[1:-1,:-2] - xc[1:-1,2:])**2 +
        (yc[1:-1,:-2] - yc[1:-1,2:])**2)**0.5
  dlnpx = param[:-2,1:-1] - param[2:,1:-1]
  dlnpy = param[1:-1,:-2] - param[1:-1,2:]
  grad_param = np.zeros(grid[0].shape)
  grad_param[1:-1,1:-1] = np.abs(dlnpx / dx) + np.abs(dlnpy / dy)
  
  return grad_param

def main():
	"""
	"""


if __name__ == "__main__":
	main()
