# Duncan Barlow, Odin project, Warwick University, 01/19
import sdf_helper as sh
import matplotlib.pyplot as plt
import numpy as np
import glob
import sys, os
import csv
from matplotlib.widgets import Slider, RadioButtons
plt.switch_backend('TkAgg')

# small number stored to avoid divide by zero errors
global small_number
small_number = 1e-100


class new_variable:
  """This class is designed to mimic an SDF variable such that the plotting
  routine can do either with the same method. It is used for creating new
  variables after analysis
  """
  def __init__(self, *args, **kwargs):
    self.data = kwargs.get('data', 1)
    self.index = kwargs.get('index', 1)
    self.all_time_data = kwargs.get('all_time_data', 1)
    self.grid =  kwargs.get('grid', 1)
    self.units = kwargs.get('units', 'Not set')
    self.name = kwargs.get('name', 'Not set')
    self.unit_conversion = kwargs.get('unit_conversion', 1)
    self.units_new = kwargs.get('units_new', 'Not set')



def basic(dat):
  """Basic analysis functions that would be useful for any simulation
  """
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

  # Add to the list of grids
  var_list = dat.grids

  var_name = "Radius_mid"
  var_list.append(var_name)
  radius = np.sqrt(xc**2 + yc**2)
  theta = np.arctan2(yc, xc)
  setattr(dat, var_name,
          new_variable(data = np.array([radius, theta]),
                       units_new = dat.Grid_Grid_mid.units_new,
                       unit_conversion = dat.Grid_Grid_mid.unit_conversion,
                       name = "Radius"))

  setattr(dat, "grids", var_list)

  # Variables that change in time and space
  var_list = dat.variables

  var_name = "None"
  var_list.insert(0, var_name)
  blank = dat.Fluid_Rho.data * 0.0 + 1.0
  setattr(dat, var_name, new_variable(data = blank,
                                      grid = dat.Grid_Grid,
                                      units_new = " ",
                                      unit_conversion = 1,
                                      name = "No Variable"))

  var_name = "Fluid_Pressure_Gbar"
  var_list.append(var_name)
  setattr(dat, var_name, new_variable(data = dat.Fluid_Pressure.data,
                                      grid = dat.Grid_Grid,
                                      units_new = "Gbar",
                                      unit_conversion = 1.0e-14,
                                      name = "Pressure"))

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
  # I believe this method has been copied from Odin source code
  len_a = np.sqrt((x[:-1,:-1] - x[:-1,1:])**2 + (y[:-1,:-1] - y[:-1,1:])**2)
  len_b = np.sqrt((x[:-1,1:] - x[1:,1:])**2 + (y[:-1,1:] - y[1:,1:])**2)
  len_c = np.sqrt((x[1:,1:] - x[1:,:-1])**2 + (y[1:,1:] - y[1:,:-1])**2)
  len_d = np.sqrt((x[1:,:-1] - x[:-1,:-1])**2 + (y[1:,:-1] - y[:-1,:-1])**2)
  angle_a = np.arctan2(np.abs(y[:-1,:-1] - y[1:,:-1]), np.abs(x[:-1,:-1] \
      - x[1:,:-1])) + np.arctan2(np.abs(y[:-1,1:] - y[:-1,:-1]),
                                 np.abs(x[:-1,1:] - x[:-1,:-1]))
  angle_c = np.arctan2(np.abs(y[:-1,1:] - y[1:,1:]),
                       np.abs(x[:-1,1:] - x[1:,1:])) \
          + np.arctan2(np.abs(y[1:,1:] - y[1:,:-1]),
                       np.abs(x[1:,1:] - x[1:,:-1]))
  area = 0.5 * len_a * len_d * np.sin(angle_a) \
       + 0.5 * len_b * len_c * np.sin(angle_c)
  setattr(dat, var_name, new_variable(data = area,
                                      grid = dat.Grid_Grid,
                                      units_new = "m$^2$",
                                      unit_conversion = 1,
                                      name = "Area"))

  var_name = "Cell_Mass"
  var_list.append(var_name)
  mass = dat.Fluid_Rho.data * vol
  setattr(dat, var_name, new_variable(data = mass,
                                      grid = dat.Grid_Grid,
                                      units_new = "kg",
                                      unit_conversion = 1,
                                      name = "Mass"))

  var_name = "Relative_Mass_Change"
  var_list.append(var_name)
  delta_mass = mass * 0.0
  delta_mass[1:,:] = mass[:-1,:] - mass[1:,:]
  rel_mass_change = delta_mass / (mass + small_number)
  max_val = 100.0
  rel_mass_change = np.where(rel_mass_change < max_val, rel_mass_change, -1.0)
  setattr(dat, var_name, new_variable(data = rel_mass_change,
                                      grid = dat.Grid_Grid,
                                      units_new = "Unitless",
                                      unit_conversion = 1,
                                      name = "Relative Mass Change"))

  var_name = "Fluid_Speed"
  var_list.append(var_name)
  speed = np.sqrt(v1**2 + v2**2 + v3**2)
  setattr(dat, var_name, new_variable(data = speed,
                                      grid = dat.Grid_Grid,
                                      units_new = "m/s",
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

  var_name = "Rho_r_reversed"
  var_list.append(var_name)
  rhor_cumulative2 = np.flip(np.cumsum(np.flip(rhor), axis=0))
  setattr(dat, var_name, new_variable(data = rhor_cumulative2,
                                      grid = dat.Grid_Grid,
                                      units_new = "g/cm$^2$",
                                      unit_conversion = 0.1,
                                      name = "Areal Density"))

  var_name = "Fluid_Number_density_ion"
  var_list.append(var_name)
  ni_density = 0.0
  if nmat == 1:
    mat_name = getattr(getattr(dat, "material_string_flags"), "data")['name']
    a_bar = getattr(getattr(dat, "material_real_flags"), "a_bar")
    mat_den = getattr(getattr(dat, "Fluid_Rho"), "data")
    ni_density = ni_density + mat_den / a_bar / amu
  else:
    for imat in range(1,nmat+1):
      mat_name = getattr(getattr(dat, "material_string_flags_" \
               + str(imat).zfill(3)), "data")['name']
      a_bar = getattr(getattr(dat, "material_real_flags_" \
            + str(imat).zfill(3)), "a_bar")
      mat_den = getattr(getattr(dat, "Fluid_Rho_"+mat_name), "data")
      ni_density = ni_density + mat_den / a_bar / amu

  setattr(dat, var_name, new_variable(data = ni_density,
                                      grid = dat.Grid_Grid,
                                      units_new = "#/m$^3$",
                                      unit_conversion = 1,
                                      name = "Number density of ions"))

  if hasattr(dat, "Fluid_Charge_State"):
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
  setattr(dat, var_name,
          new_variable(data = com,
                       units_new = dat.Grid_Grid_mid.units_new,
                       unit_conversion = dat.Grid_Grid_mid.unit_conversion,
                       name = "Centre of Mass"))

  setattr(dat, "variables_time", var_list)

  return dat



def laser(dat, *args, **kwargs):
  """Analysis of variables that are linked to the laser
  """
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
  crit_surf_ind = np.zeros(ny, dtype=int)
  crit_rad = np.zeros(ny)
  quart_crit_surf_ind = np.zeros(ny, dtype=int)
  quart_crit_rad = np.zeros(ny)
  # critical surface is chosen based on direction of laser propagation
  # This needs upgrading to use the output laser direction.
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
  setattr(dat, var_name,
          new_variable(data = crit_rad,
                       index = crit_surf_ind,
                       units_new = dat.Grid_Grid.units_new,
                       unit_conversion = dat.Grid_Grid.unit_conversion,
                       name = "Location of critical surface"))

  var_name = "Critical_Surface_quarter"
  var_list.append(var_name)
  setattr(dat, var_name,
          new_variable(data = quart_crit_rad,
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
  density_scale_length[1:-1,1:-1] = abs(ne_density[1:-1,1:-1] \
                                  / (grad_ne_density[1:-1,1:-1] + small_number))
  # Remeber unit conversions are applied after!!
  max_val = 1e-2
  density_scale_length = np.where(density_scale_length < max_val,
                                  density_scale_length, 0.0)
  setattr(dat, var_name,
          new_variable(data = density_scale_length,
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
                                      units_new = 'kJ',
                                      unit_conversion = 1.0e-3,
                                      name = "Total laser energy deposited"))

  setattr(dat, "variables_time", var_list)

  return dat



def hot_electron(dat, *args, **kwargs):
  """
  """
  laser_change = kwargs.get('laser_change', False)
  sdf_num = kwargs.get('sdf_num', 0)
  istart = kwargs.get('istart', 0)
  pathname = kwargs.get('pathname', os.path.abspath(os.getcwd()))

  electron_dep = dat.Fluid_Energy_deposited_hot_electron.data

  # Variables that change in time and space
  var_list = dat.variables

  var_name = "Electron_energy_per_step"
  var_list.append(var_name)
  if sdf_num == istart:
    electron_dep_step = electron_dep
    dt = dat.Header.get('time')
  elif sdf_num >= istart:
    SDFName=pathname+'/'+str(sdf_num-1).zfill(4)+'.sdf'
    dat2 = sh.getdata(SDFName,verbose=False)
    electron_dep_step = electron_dep \
        - dat2.Fluid_Energy_deposited_hot_electron.data
    dt = dat.Header.get('time') - dat2.Header.get('time')
  else:
    print('Error with electron change calculation')
    print('sdf_num = ', sdf_num, ' and the minimum = ', istart)
  setattr(dat, var_name, new_variable(data = electron_dep_step,
                                      grid = dat.Grid_Grid,
                                      units_new = "J/kg",
                                      unit_conversion = 1,
                                      name = "Hot Electron Energy Deposited"))

  var_name = "Electron_Energy_per_cell"
  var_list.append(var_name)
  electron_dep_cell = electron_dep_step * dat.Cell_Mass.data
  setattr(dat, var_name, new_variable(data = electron_dep_cell,
                                      grid = dat.Grid_Grid,
                                      units_new = "J",
                                      unit_conversion = 1,
                                      name = "Hot Electron Energy per Cell"))

  var_name = "Electron_Power_per_volume"
  var_list.append(var_name)
  if dt < small_number:
    dt = 1.0
  electron_pwr_per_vol = electron_dep_step * dat.Cell_Mass.data \
      / dat.Fluid_Volume_rz.data / dt
  setattr(dat, var_name, new_variable(data = electron_pwr_per_vol,
                                      grid = dat.Grid_Grid,
                                      units_new = "W/m$^3$",
                                      unit_conversion = 1,
                                      name = "Hot Electron Power Per Volume"))

  setattr(dat, "variables", var_list)

  # variables that only change in time
  var_list = dat.variables_time

  var_name = "Hot_Electron_Power_Total_Deposited"
  var_list.append(var_name)
  tot_ele_dep = np.sum(np.sum(electron_dep_step * dat.Cell_Mass.data)) / dt
  setattr(dat, var_name,
          new_variable(data = tot_ele_dep,
                       units_new = 'TW',
                       unit_conversion = 1.0e-12,
                       name = "Total Hot Electron Power Deposited"))

  var_name = "Hot_Electron_Energy_Total_Deposited"
  var_list.append(var_name)
  tot_electron_dep = np.sum(np.sum(dat.Cell_Mass.data * electron_dep))
  setattr(dat, var_name, new_variable(data = tot_electron_dep,
                                      units_new = 'kJ',
                                      unit_conversion = 1.0e-3,
                                      name = "Total Hot Electron Energy Deposited"))

  setattr(dat, "variables_time", var_list)

  return dat



def adiabat(dat, *args, **kwargs):
  """Calculate paramaters pertaining to shocks: the fluid adiabat and inverse
  pressure length scale.
  """
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
  # See "The Physics of Inertial Fusion" by Atzeni and Meyer-ter-Vehn
  # page number pending (sorry!)
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
  pressure_ls = gradient_function(np.log(pressure + small_number),
                                  dat.Grid_Grid_mid.data)
  setattr(dat, var_name, new_variable(data = pressure_ls,
                                      grid = dat.Grid_Grid,
                                      units_new = "unitless",
                                      unit_conversion = 1,
                                      name = "Inverse Pressure Length Scale"))

  setattr(dat, "variables", var_list)

  return dat



def energy(dat, *args, **kwargs):
  """Energy calculations that are not already accounted for by Odin
  """
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
  if "Laser_Energy_Total_Deposited" in dat.variables_time:
    tot_LE = dat.Laser_Energy_Total_Deposited.data
  else:
    tot_LE = 0.0
  if "Hot_Electron_Energy_Total_Deposited" in dat.variables_time:
    tot_HEE = dat.Hot_Electron_Energy_Total_Deposited.data
  else:
    tot_HEE = 0.0
  tot_energy = tot_IE + tot_KE - tot_LE - tot_HEE
  setattr(dat, var_name, new_variable(data = tot_energy,
                                      units_new = 'J',
                                      unit_conversion = 1,
                                      name = "Total Energy"))

  setattr(dat, "variables_time", var_list)

  return dat



def time_variables(dat, *args, **kwargs):
  """This routine calculates variables that have no spatial dependance but
  do change in time.
  """
  # fraction of total sphere
  frac_sphere = 1.0
  if dat.Logical_flags.use_rz:
    frac_sphere = (dat.Real_flags.y_max - dat.Real_flags.y_min) / 2.0

  laser_profile_filename = "laser_profile.csv"
  if os.path.exists(laser_profile_filename):

    var_name = "Input_laser_profile"
    with open(laser_profile_filename) as csv_file:
      csv_reader = csv.reader(csv_file,delimiter=',')
      laser_profile_time = []
      laser_profile_pwr = []
      for row in csv_reader:
        laser_profile_time.append(float(row[0]))
        laser_profile_pwr.append(float(row[1]))
    laser_profile_time = np.asarray(laser_profile_time)
    laser_profile_pwr = np.asarray(laser_profile_pwr) * frac_sphere
    setattr(dat, var_name, new_variable(data = 0.0,
                                        all_time_data = laser_profile_pwr,
                                        units_new = 'TW',
                                        unit_conversion = 1.0e-12,
                                        name = "Original Laser Profile"))
    setattr(getattr(dat,var_name), "times", laser_profile_time)
    setattr(getattr(dat,var_name), "times_units", 'ns')
    setattr(getattr(dat,var_name), "times_conversion", 1.0e9)

    var_name = "Input_laser_profile_energy"
    profile_length = len(laser_profile_time)
    laser_profile_energy = np.zeros(profile_length)
    for i in range(1, profile_length):
      laser_profile_energy[i] = laser_profile_energy[i-1] + \
      (laser_profile_time[i] - laser_profile_time[i-1]) * \
      (0.5 * abs(laser_profile_pwr[i] - laser_profile_pwr[i-1]) + \
      min(laser_profile_pwr[i], laser_profile_pwr[i-1]))
    setattr(dat, var_name, new_variable(data = 0.0,
                                        all_time_data = laser_profile_energy,
                                        units_new = 'kJ',
                                        unit_conversion = 1.0e-3,
                                        name = "Original Laser Profile Energy"))
    setattr(getattr(dat,var_name), "times", laser_profile_time)
    setattr(getattr(dat,var_name), "times_units", 'ns')
    setattr(getattr(dat,var_name), "times_conversion", 1.0e9)
  else:
    print(laser_profile_filename, " DOES NOT exist")

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
                                      units_new = 'TW',
                                      unit_conversion = 1.0e-12,
                                      name = "Total laser power deposited"))

  setattr(dat, "variables_time", var_list)

def gradient_function(param, grid):
  """ This function calculates the 2D gradient of the [param] relative
  to the [grid]
  """
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
