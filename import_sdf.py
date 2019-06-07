# Duncan Barlow, Odin project, Warwick University, 01/19
import sdf_helper as sh
import matplotlib.pyplot as plt
import numpy as np
import glob
import sys, os
from matplotlib.widgets import Slider, RadioButtons
import analysis_functions as afunc



def add_label(dat):
  
  # default labels
  for var_name in dat.variables:
    var = getattr(dat, var_name)
    setattr(var, "unit_conversion", 1)
    units = getattr(var, "units")
    setattr(var, "units_new", units)
  
  for var_name in dat.grids:
    var = getattr(dat, var_name)
    setattr(var, "unit_conversion", 1)
    units = getattr(var, "units")
    setattr(var, "units_new", units[0])
  
  #print("Warning: User defined units will overwrite original")
  #print("and the conversion might only be true from SI.")
  
  # User defined
  var = getattr(dat, "Fluid_Temperature_electron")
  name = 'Electron Temperature'
  setattr(var, "name", name)
  unit_conversion = 1.0 / 11604.5 / 1000.0 # from Kelvin
  setattr(var, "unit_conversion", unit_conversion)
  units_new = "KeV"
  setattr(var, "units_new", units_new)
  
  var = getattr(dat, "Fluid_Temperature_ion")
  name = 'Ion Temperature'
  setattr(var, "name", name)
  unit_conversion = 1.0 / 11604.5 / 1000.0 # from Kelvin
  setattr(var, "unit_conversion", unit_conversion)
  units_new = "KeV"
  setattr(var, "units_new", units_new)
  
  var = getattr(dat, "Grid_Grid_mid")
  name = 'Distance, cell centred'
  setattr(var, "name", name)
  unit_conversion = 1.0e6 # from m
  setattr(var, "unit_conversion", unit_conversion)
  units_new = "$\mu m$"
  setattr(var, "units_new", units_new)
  
  var = getattr(dat, "Grid_Grid")
  name = 'Distance, node centred'
  setattr(var, "name", name)
  unit_conversion = 1.0e6 # from m
  setattr(var, "unit_conversion", unit_conversion)
  units_new = "$\mu m$"
  setattr(var, "units_new", units_new)
  
  var = getattr(dat, "Fluid_Rho")
  name = 'Density'
  setattr(var, "name", name)
  unit_conversion = 1.0 / 1000.0 # from kg to g
  setattr(var, "unit_conversion", unit_conversion)
  units_new = "g/cm$^3$"
  setattr(var, "units_new", units_new)
  
  var = getattr(dat, "Fluid_Pressure_ion")
  name = 'Ion Pressure'
  setattr(var, "name", name)
  unit_conversion = 1.0e-11 # from Pascal to Mbar
  setattr(var, "unit_conversion", unit_conversion)
  units_new = "Mbar"
  setattr(var, "units_new", units_new)
  
  var = getattr(dat, "Fluid_Pressure_electron")
  name = 'Electron Pressure'
  setattr(var, "name", name)
  unit_conversion = 1.0e-11 # from Pascal to Mbar
  setattr(var, "unit_conversion", unit_conversion)
  units_new = "Mbar"
  setattr(var, "units_new", units_new)
  
  var = getattr(dat, "Fluid_Pressure")
  name = 'Total Pressure'
  setattr(var, "name", name)
  unit_conversion = 1.0e-11 # from Pascal to Mbar
  setattr(var, "unit_conversion", unit_conversion)
  units_new = "Mbar"
  setattr(var, "units_new", units_new)
  
  return dat



def preallocate_dat(dat, iend, cs):
  
  for var_name in dat.grids:
    var = getattr(dat, var_name)
    data = getattr(var, "data")
    len_x = np.shape(data[0])[0]
    array = np.zeros((2, iend, len_x))
    array[0,0,:] = data[0][:,cs]
    array[1,0,:] = data[1][:,cs]
    setattr(var, "all_time_data", array)
    setattr(var, "all_time_data_polar", array)
  
  for var_name in dat.variables:
    var = getattr(dat, var_name)
    data = getattr(var, "data")
    len_x = np.shape(data)[0]
    array = np.zeros((iend, len_x))
    array[0,:] = data[:,cs]
    setattr(var, "all_time_data", array)
  
  for var_name in dat.variables_time:
    var = getattr(dat, var_name)
    data = getattr(var, "data")
    array = np.zeros(iend)
    array[0] = data
    setattr(var, "all_time_data", array)
  
  return dat



def use_sdf(sdf_num, pathname, *args, **kwargs):
  istart = kwargs.get('istart', 0)
  use_analysis = kwargs.get('use_analysis', False)

  SDFName=pathname+'/'+str(sdf_num).zfill(4)+'.sdf'
  dat = sh.getdata(SDFName,verbose=False)
  
  dat_names = list(dat.__dict__.keys())
  variable_type = type(dat.Fluid_Rho)
  grid_type = type(dat.Grid_Grid)
  
  dat_variable_names = []
  dat_grid_names = []
  dat_variable_time_names = []
  for n in range(0, len(dat_names)):
    var = getattr(dat, dat_names[n])
    if type(var) == variable_type:
      dat_variable_names.append(dat_names[n])
    elif type(var) == grid_type:
      dat_grid_names.append(dat_names[n])
  
  setattr(dat, "grids", dat_grid_names)
  setattr(dat, "variables", dat_variable_names)
  setattr(dat, "variables_time", dat_variable_time_names)
		
  # Save time variable
  var_list = dat.variables_time
  var_name = "Times"
  var_list.append(var_name)
  setattr(dat, var_name, afunc.new_variable(data = dat.Header["time"],
                                            units_new = "ns",
                                            unit_conversion = 1.0e9,
                                            name = "Time"))
  setattr(dat, "variables_time", var_list)
  
  dat = add_label(dat)
  
  if use_analysis:
    dat = afunc.basic(dat)
    dat = afunc.laser(dat, call_basic = False, laser_change = True,
        sdf_num = sdf_num, istart = istart, pathname = pathname)
    dat = afunc.adiabat(dat, call_basic = False)
    dat = afunc.energy(dat, call_basic = False)
    # IFAR
  
  return dat



def get_data_all(dat1, istart, iend, pathname, use_analysis, cs):
  """
  """
  irange = iend - istart + 1
  dat1 = preallocate_dat(dat1, irange, cs)
  
  for n in range(istart+1, irange):     
    dat = use_sdf(n, pathname, use_analysis = use_analysis)
    
    # grid for this data is either radius or X depending on rz t/f
    for var_name in dat.grids:
    	array = getattr(getattr(dat1, var_name), "all_time_data")
    	data = getattr(getattr(dat, var_name), "data")
    	array[0,n,:] = data[0][:,cs]
    	array[1,n,:] = data[1][:,cs]
    	setattr(getattr(dat1, var_name), "all_time_data", array)
    
    for var_name in dat.variables:
      array = getattr(getattr(dat1, var_name), "all_time_data")
      data = getattr(getattr(dat, var_name), "data")
      array[n,:] = data[:,cs]
      setattr(getattr(dat1, var_name), "all_time_data", array)
    
    for var_name in dat.variables_time:
      array = getattr(getattr(dat1, var_name), "all_time_data")
      data = getattr(getattr(dat, var_name), "data")
      array[n] = data
      setattr(getattr(dat1, var_name), "all_time_data", array)
  
  if use_analysis:
    afunc.time_variables(dat1)
  
  return dat1



def main():
        """
        """


if __name__ == "__main__":
        main()
