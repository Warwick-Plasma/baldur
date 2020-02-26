# Duncan Barlow, Odin project, Warwick University, 01/19
import sdf_helper as sh
import sdf
import matplotlib.pyplot as plt
import numpy as np
import glob
import sys, os
from matplotlib.widgets import Slider, RadioButtons
import analysis_functions as afunc
import time
plt.switch_backend('TkAgg')


def add_label(dat):
  """When plot the data can have different units, this routine makes sure all
  the Odin output variables have the [unit_converion] and [units_new] labels
  that are required by odin_plot.py for plots. It also allows the user to
  input their own unit changes for default variable, unit changes for post-
  processed variables is done in analysis_function.py.
  """
  # default labels for varaibles
  for var_name in dat.variables:
    var = getattr(dat, var_name)
    setattr(var, "unit_conversion", 1)
    units = getattr(var, "units")
    setattr(var, "units_new", units)

  # default labels for grids
  for var_name in dat.grids:
    var = getattr(dat, var_name)
    setattr(var, "unit_conversion", 1)
    units = getattr(var, "units")
    setattr(var, "units_new", units[0])

  # Epoch grids require change from a length N and a length M array
  # to a meshgrid of size N by M.
  if (dat.Header['code_name'] == 'Epoch2d'):
    for var_name in dat.grids:
      var = getattr(dat, var_name)
      data_x = getattr(var, "data")[0]
      data_y = getattr(var, "data")[1]
      if (len(np.shape(data_x)) == 1):
        data_X, data_Y = np.meshgrid(data_x, data_y, indexing='ij')
        data = np.array([data_X, data_Y])
        setattr(var, "data", data)

  # User defined change in units
  if (dat.Header['code_name'] == 'Odin2D'):
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
  """Preallocates empty arrays to [dat] file variables such that a full
  history of the varaible can be added later. The history is taken at
  cross-section [cs].
  """
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
  """The basic routine for taking data from an sdf file and putting it in
  [dat]. The varaibles are put into lists depending on type so they can be
  easily plotted later, the lists are cleaned of certain variables that
  cannot easily be plotted. Finally the data is labelled with new units
  and passed to the analysis functions for post-processing.
  """
  istart = kwargs.get('istart', 0)
  use_analysis = kwargs.get('use_analysis', False)

  SDFName=pathname+'/'+str(sdf_num).zfill(4)+'.sdf'
  dat = sh.getdata(SDFName,verbose=False)

  # Get all variables
  dat_names = list(dat.__dict__.keys())
  variable_type = sdf.BlockPlainVariable
  grid_type = sdf.BlockLagrangianMesh
  beam_type = sdf.BlockStitchedPath

  # Create lists of certain types of variable for plotting
  dat_grid_names = []
  dat_variable_names = []
  dat_variable_time_names = []
  dat_track_surfaces = []
  dat_beam_names = []
  dat_burst_names = []
  for n in range(0, len(dat_names)):
    var = getattr(dat, dat_names[n])
    if type(var) == variable_type:
      dat_variable_names.append(dat_names[n])
    elif type(var) == grid_type:
      dat_grid_names.append(dat_names[n])
    elif type(var) == beam_type:
      dat_beam_names.append(dat_names[n])

  # Clean grid list of things that require special effort to plot
  bad_var_list = []
  for var in dat_grid_names:
    if ('Ray' in var):
      bad_var_list.append(var)
    if ('Electrons_Electron' in var):
      bad_var_list.append(var)
  for var in bad_var_list:
    dat_grid_names.remove(var)

  # Clean variable list of things that require special effort to plot
  bad_var_list = []
  for var in dat_variable_names:
    if ('Ray' in var):
      bad_var_list.append(var)
    if ('Electrons_Electron' in var):
      bad_var_list.append(var)
  for var in bad_var_list:
    dat_variable_names.remove(var)

  # Clean beams list of individual rays
  bad_var_list = []
  for var in dat_beam_names:
    if ('_' in var):
      bad_var_list.append(var)
    else:
      if ('Burst' in var):
        dat_burst_names.append(var)
        bad_var_list.append(var)
  for var in bad_var_list:
    dat_beam_names.remove(var)

  # Clean Epoch list of CPU and particle distributions that I cannot plot yet
  if (dat.Header['code_name'] == 'Epoch2d'):
    dat_grid_names = ['Grid_Grid', 'Grid_Grid_mid']
    bad_var_list = []
    for var in dat_variable_names:
      if ('CPU' in var):
        bad_var_list.append(var)
      if ('dist' in var):
        bad_var_list.append(var)
    for var in bad_var_list:
      dat_variable_names.remove(var)

  # Add lists to dat
  setattr(dat, "grids", dat_grid_names)
  setattr(dat, "track_surfaces", dat_track_surfaces)
  setattr(dat, "variables", dat_variable_names)
  setattr(dat, "variables_time", dat_variable_time_names)
  setattr(dat, "beams", dat_beam_names)
  setattr(dat, "bursts", dat_burst_names)

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

  # Analysis functions are applied to SDF data
  if use_analysis:
    dat = afunc.basic(dat)
    dat = afunc.laser(dat, call_basic = False, laser_change = True,
        sdf_num = sdf_num, istart = istart, pathname = pathname)
    dat = afunc.adiabat(dat, call_basic = False)
    dat = afunc.energy(dat, call_basic = False)

  return dat



def get_data_all(dat1, istart, iend, pathname, use_analysis, cs):
  """This routine is a wrapper for [use_sdf] which it applies to all time sdf
  files in current directory and compile cross-sections [cs] into a single
  object [dat1]. [istart] and [iend] dictate the bounds of the analysis.
  """
  irange = iend - istart + 1
  dat1 = preallocate_dat(dat1, irange, cs)

  # loop over sdf files in range [istart] to [iend].
  for n in range(1, irange):
    print_string = 'Processing file {:4d}'.format(istart+n) + \
        ' of {:4d}'.format(iend)
    sys.stdout.write('\r' + print_string)
    sys.stdout.flush()

    dat = use_sdf(istart+n, pathname, use_analysis = use_analysis)

    # Assemble arrays with position X and Y data stored for all time
    for var_name in dat.grids:
    	array = getattr(getattr(dat1, var_name), "all_time_data")
    	data = getattr(getattr(dat, var_name), "data")
    	array[0,n,:] = data[0][:,cs]
    	array[1,n,:] = data[1][:,cs]
    	setattr(getattr(dat1, var_name), "all_time_data", array)

    # Assemble arrays with variable data stored for all time
    for var_name in dat.variables:
      array = getattr(getattr(dat1, var_name), "all_time_data")
      data = getattr(getattr(dat, var_name), "data")
      array[n,:] = data[:,cs]
      setattr(getattr(dat1, var_name), "all_time_data", array)

    # Assemble lists of numbers that quantify each sdf file ie. total KE
    for var_name in dat.variables_time:
      array = getattr(getattr(dat1, var_name), "all_time_data")
      data = getattr(getattr(dat, var_name), "data")
      array[n] = data
      setattr(getattr(dat1, var_name), "all_time_data", array)

  # special analysis functions applied to time series data
  if use_analysis:
    afunc.time_variables(dat1)

  return dat1



def main():
        """
        """


if __name__ == "__main__":
        main()
