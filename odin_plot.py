# Duncan Barlow, Odin project, Warwick University, 01/19
import sdf_helper as sh
import matplotlib.pyplot as plt
import numpy as np
import glob
import import_sdf as isdf
import sys, os
import csv
from matplotlib import cm
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.widgets import Slider, RadioButtons
plt.switch_backend('TkAgg')

# This sets a global fontsize
global fs
fs = 15

# Globals used to stop dvide by zero errors
global small_num
global big_num
small_num = 1e-100
big_num = 1e100



def plot_thermodynamic_path(fig, ax, num, x_data, y_data, x_label, y_label, line_label):
  """
  """
  ax.lines[num].set_visible(True)

  ax.xaxis.get_offset_text().set_size(fs)
  ax.yaxis.get_offset_text().set_size(fs)

  ax.set_ylabel(y_label, fontsize = fs)
  ax.set_xlabel(x_label, fontsize = fs)

  if num == 1:
    ax.lines[num].set_linestyle('-')

  ax.lines[num].set_ydata(y_data)
  ax.lines[num].set_xdata(x_data)

  ax.tick_params(axis='x', labelsize = fs)
  ax.tick_params(axis='y', labelsize = fs)

  ax.set_xlim(np.min(x_data[:-1]), np.max(x_data[:-1]))
  ax.set_ylim(np.min(y_data[:-1]), 1.3 * np.max(y_data[:-1]))

  try:
    line_colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
  except:
    line_colours = plt.rcParams['axes.color_cycle']
 
  ax.lines[num].set_color(line_colours[num])
  ax.lines[num].set_color(line_colours[num])
  ax.set_yscale('log')
  ax.set_xscale('log')

  ax.lines[num].set_label(line_label)



def save_as_csv(file_name, x_data, y_data, *args, **kwargs):
  x_label = kwargs.get('xlabel', 'x label (units)')
  y_label = kwargs.get('ylabel', 'y label (units)')

  data = [None] * (len(x_data)+1)
  with open(file_name+".csv", "w", newline='') as file:

    writer = csv.writer(file)
    data[0] = x_label
    data[1:] = x_data
    writer.writerow(data)
    data[0] = y_label
    data[1:] = y_data
    writer.writerow(data)

def plot_laser_profile(*args, **kwargs):
  """Simple, Self-suffient routine to plot a .csv file as a laser profile
  """
  name1 = kwargs.get('name1', 'laser_profile.csv')
  name2 = kwargs.get('name2', 'laser_profile.csv')

  with open(name1) as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    times = []
    powers = []
    for row in readCSV:
      time = float(row[0]) * 1.0e9
      times.append(time)
      power = float(row[1]) * 1.0e-12
      powers.append(power)
  plt.figure()
  ax = plt.axes()
  ax.plot(times, powers, linewidth = 2)

  with open(name2) as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    times = []
    powers = []
    for row in readCSV:
      time = float(row[0]) * 1.0e9
      times.append(time)
      power = float(row[1]) * 1.0e-12
      powers.append(power)
  ax.plot(times, powers, color = 'tab:red', linestyle = "--", linewidth = 2)

  ax.set_xlabel('Time (ns)', fontsize = fs)
  ax.set_ylabel('Power (TW)', fontsize = fs)
  ax.tick_params(axis='x', labelsize = fs)
  ax.tick_params(axis='y', labelsize = fs)

  plt.show()



def time_history(fig, ax1, cax1, *args, **kwargs):
  """A pcolormesh plot of space against time with a variable shown in colour
  """
  parameters = kwargs.get('parameters', plot_parameters())
  data_struct = kwargs.get('data', data_structure())

  var = getattr(data_struct.data, parameters.var_name)
  unit_conv = getattr(var, "unit_conversion")
  units = getattr(var, "units_new")
  name = getattr(var, "name")
  grid = getattr(var, "grid")
  grid_name = getattr(grid, "name")
  grid_data = getattr(grid, "all_time_data")
  grid_conv = getattr(grid, "unit_conversion")
  grid_units = getattr(grid, "units_new")

  c_data = getattr(var, "all_time_data") * unit_conv
  y_data, c_data = two_dim_grid(data_struct.data, c_data)
  times = data_struct.data.Times
  x_data, y_data1 = np.meshgrid(times.all_time_data \
                  * times.unit_conversion, y_data[0,:], indexing='ij')

  x_label = times.name + ' (' + getattr(times, 'units_new') + ')'
  y_label = 'Radius (' + grid_units + ')'
  c_label = name + " (" + units + ")"

  if (parameters.grid_choice == 'default'):
    y_data = y_data
  elif (parameters.grid_choice == 'initial'):
    y_data = y_data1
    y_label = 'Initial ' + y_label
  elif (parameters.grid_choice == 'cell number'):
    pos = np.linspace(0, np.shape(c_data)[1]-1, np.shape(c_data)[1])
    _, y_data = np.meshgrid(times.all_time_data, pos, indexing='ij')
    y_label = 'Cell Number'

  cbar_range = np.max(c_data) - np.min(c_data) + small_num
  cbar_max = np.min(c_data) + np.exp(np.log(cbar_range) \
           + parameters.cbar_colour_scale)

  if parameters.reset_axis:
    zoomed_axis1 = np.array([np.min(x_data[:-1,:]), np.max(x_data[:-1,:]),
                             np.min(y_data[:-1,:]), np.max(y_data[:-1,:])])
  else:
    zoomed_axis1 = np.array([ax1.get_xlim()[0], ax1.get_xlim()[1],
                             ax1.get_ylim()[0], ax1.get_ylim()[1]])

  ax1.clear() # This is nessasary for speed
  cax1.clear()

  cmesh = ax1.pcolormesh(x_data, y_data, c_data, linewidth=0.1, shading='auto')
  cmesh.set_clim(np.min(c_data), cbar_max)

  cbar = fig.colorbar(cmesh, cax=cax1)

  ax1.set_xlabel(x_label, fontsize = fs)
  ax1.set_ylabel(y_label, fontsize = fs)
  cbar.set_label(c_label, fontsize = fs)

  ax1.tick_params(axis='x', labelsize = fs)
  ax1.tick_params(axis='y', labelsize = fs)
  cbar.ax.tick_params(labelsize=fs)

  cbar.ax.yaxis.get_offset_text().set(size = fs)
  ax1.xaxis.get_offset_text().set_size(fs)
  ax1.yaxis.get_offset_text().set_size(fs)

  ax1.set_xlim(zoomed_axis1[:2])
  ax1.set_ylim(zoomed_axis1[2:])
  cbar.draw_all()

  plt.show()



def two_dim_grid(dat, data):
  """ This routine is called when making a time vs radius plot for the time
  history. It creates the spatial and time grid.
  """
  x_mid = dat.Grid_Grid_mid.all_time_data[0] * dat.Grid_Grid_mid.unit_conversion
  y_mid = dat.Grid_Grid_mid.all_time_data[1] * dat.Grid_Grid_mid.unit_conversion
  grid = np.sqrt(x_mid**2 + y_mid**2)

  x_edge = dat.Grid_Grid.all_time_data[0] * dat.Grid_Grid_mid.unit_conversion
  y_edge = dat.Grid_Grid.all_time_data[1] * dat.Grid_Grid_mid.unit_conversion
  grid_edge = np.sqrt(x_edge**2 + y_edge**2)

  if np.shape(data) != np.shape(grid):
    if np.shape(data) == np.shape(grid_edge):
      grid = grid_edge
    else:
      print("Unknown geometry variable")
      print("Creating uniform grid")
      pos = np.linspace(0, np.shape(data)[1]-1, np.shape(data)[1])
      times, grid = np.meshgrid(dat.Times.all_time_data, pos, indexing='ij')

  return grid, data



def time_history_lineout(fig, ax, ax1, *args, **kwargs):
  """ A 1D line of time vs amplitude for the tkiter selected variable. Only
  called from time history menu.
  """
  parameters = kwargs.get('parameters', plot_parameters())
  data_struct = kwargs.get('data', data_structure())

  if parameters.use_analysis:

    ax1.lines[0].set_visible(True)

    var = getattr(data_struct.data, parameters.var_name2)
    unit_conv = getattr(var, "unit_conversion")
    units = getattr(var, "units_new")
    name = getattr(var, "name")
    y_data1 = getattr(var, "all_time_data") * unit_conv
    y_data = y_data1

    times = data_struct.data.Times
    x_data = times.all_time_data * times.unit_conversion
    ax.lines[1].set_xdata(x_data)
    ax1.lines[0].set_xdata(x_data)
    ax1.lines[0].set_ydata(y_data1)

    x_label = times.name + " (" + times.units_new + ")"
    y_label = name + " (" + units + ")"

    ax.set_xlabel(x_label, fontsize = fs)
    ax1.set_ylabel(y_label, color='tab:red', fontsize = fs)

    if hasattr(data_struct.data, "Laser_Power_Total_Deposited"):
      ax.lines[0].set_visible(True)

      var = data_struct.data.Laser_Power_Total_Deposited
      y_data = var.all_time_data * var.unit_conversion
      name = var.name
      units = var.units_new

      ax.lines[0].set_ydata(y_data)
      ax.lines[0].set_xdata(x_data)
      y_label = name + " (" + units + ")"
      ax.set_ylabel(y_label, fontsize = fs)

    ax.xaxis.get_offset_text().set_size(fs)
    ax.yaxis.get_offset_text().set_size(fs)
    ax.tick_params(axis='x', labelsize = fs)
    ax.tick_params(axis='y', labelsize = fs)
    ax1.tick_params(axis='y', labelsize = fs)
    ax1.yaxis.get_offset_text().set_size(fs)

    ax.set_xlim(np.min(x_data[:-1]), np.max(x_data[:-1]))
    ax.set_ylim(np.min(y_data[:-1]), 1.3 * np.max(y_data[:-1]))
    ax1.set_xlim(np.min(x_data[:-1]), np.max(x_data[:-1]))
    ax1.set_ylim(np.min(y_data1[:-1]), 1.3 * np.max(y_data1[:-1]))

    if hasattr(data_struct.data, "Input_laser_profile"):
      ax.lines[1].set_visible(True)
      var = getattr(data_struct.data, "Input_laser_profile")
      x_data = var.times * var.times_conversion
      ax.lines[1].set_xdata(x_data)
      y_data = var.all_time_data * var.unit_conversion
      ax.lines[1].set_ydata(y_data)
      if (parameters.var_name2 == "Laser_Energy_Total_Deposited"):
        ax1.lines[1].set_visible(True)
        var = getattr(data_struct.data, "Input_laser_profile_energy")
        x_data = var.times * var.times_conversion
        ax1.lines[1].set_xdata(x_data)
        y_data = var.all_time_data * var.unit_conversion
        ax1.lines[1].set_ydata(y_data)
      else:
        ax1.lines[1].set_visible(False)
    else:
      ax.lines[1].set_visible(False)

    plt.show()



def check_analysis(use_analysis):
  """Checks the terminal input for whether analysis of sdf data is
  requested.
  """
  if use_analysis == True:
    print("starting analysis")
  elif use_analysis == False:
    print("set: <use_analysis = True>, for analysis")
  else:
    print("set: <use_analysis = True> or <False>")
    print("it requires certain dump_masks")
    sys.exit()



def data_and_plot(sdf_num, fig, ax1, cax1, fig2, ax2, ax3, parameters):
  """This routine is called from menu.py and is a wrapper for calls to data
  collection the plotting routines. The dat file is created with all the
  data from the sdf file.
  """
  print_string = 'Processing file {:4d}'.format(sdf_num) \
               + ' of {:4d}'.format(parameters.iend) + '   '
  sys.stdout.write('\r' + print_string)
  sys.stdout.flush()
  data_struct = data_structure()
  data_struct.data = [None] * parameters.num_dir

  for num in range(0, parameters.num_dir):
    if parameters.movie:
      parameters.sdf_num = sdf_num
      parameters.sdf_num2 = sdf_num
    else:
      if num == 0:
        sdf_num = parameters.sdf_num
      else:
        sdf_num = parameters.sdf_num2
    if parameters.apply_comparison[num]:
      if os.path.isdir(parameters.entry_comparison[num]):
        #print("Data loaded from: " + parameters.entry_comparison[num] + " to " + num)
        data_struct.data[num] = isdf.use_sdf(sdf_num,
                                             parameters.entry_comparison[num],
                                             use_analysis = parameters.use_analysis,
                                             istart = parameters.istart)
      else:
        parameters.apply_comparison[num] = False
        print()
        print("Warning: " + parameters.entry_comparison[num] + " is not a directory")

  snapshot(fig, ax1, cax1, parameters = parameters, data = data_struct)
  lineout(fig2, ax2, ax3, parameters = parameters, data = data_struct)

  return parameters, data_struct



# © Copyright 2002 - 2012 John Hunter, Darren Dale, Eric Firing, Michael
# Droettboom and the Matplotlib development team; 2012 - 2018 The Matplotlib
# development team.
# https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/multicolored_line.html
def plot_colourline(fig1, ax1, x, y, c, cnorm):
  """Modified from matplotlib examples. This plotting routine uses
  LineCollection to plot a line with different segments in different
  colours.
  """
  # Create a set of line segments so that we can color them individually
  # This creates the points as a N x 1 x 2 array so that we can stack points
  # together easily to get the segments. The segments array for line collection
  # needs to be (numlines) x (points per line) x 2 (for x and y)
  points = np.array([x, y]).T.reshape(-1, 1, 2)
  segments = np.concatenate([points[:-1], points[1:]], axis=1)

  # Create a continuous norm to map from data points to colors
  lc = LineCollection(segments, cmap='Reds', norm=cnorm)
  # Set the values used for colormapping
  lc.set_array(c)
  lc.set_linewidth(2)
  ax1.add_collection(lc)



def plot_rays(name, name_var, skip, dat, fig1, ax1, use_polar, grid_conv):
  """Adds rays to the 2d plot. Currently does not have a colourbar. Colour
  of the line segments is determined by the [name_var] that is passed in and
  which type of ray is determined by [name].
  """
  beam = getattr(dat, name)
  beam_energy = getattr(dat, name + '_' + name_var)
  nrays = len(beam.data)
  for iray in range(0, nrays, skip):
    print_string = 'Processing ray {:4d}'.format(iray+1) \
                 + ' of {:4d}'.format(nrays) + '   '
    sys.stdout.write('\r' + print_string)
    sys.stdout.flush()

    x_ray = beam.data[iray].data[0] * grid_conv
    y_ray = beam.data[iray].data[1] * grid_conv
    c_ray = beam_energy.data[iray].data
    cmax = np.max(beam_energy.data[iray].data)
    cmin = np.min(beam_energy.data[iray].data)
    cnorm = plt.Normalize(cmin, cmax)

    if use_polar: x_ray, y_ray, _, _ = polar_coordinates(x_ray, y_ray, " ")

    plot_colourline(fig1, ax1, x_ray, y_ray, c_ray, cnorm)
  smap = cm.ScalarMappable(norm=cnorm, cmap='Reds')
  smap.set_array([])
  #fig1.colorbar(smap)



class data_structure:
  """A class to save the Odin data
  """
  def __init__(self):
    self.data = [None]



class plot_parameters:
  """A class defined to save all the values generated by the tkinter controls.
  """
  def __init__(self):
    self.first_call = True
    self.sdf_num = 0
    self.use_analysis = False
    self.pathname = 'None'
    self.istart = 0
    self.iend = 0
    self.grid_boolean = False
    self.use_polar = False
    self.var_name = ['None'] * 2
    self.reset_axis = False
    self.view_anisotropies = False
    self.use_log = False
    self.surface_name = 'None'
    self.apply_scale_max = False
    self.scale_max = 1.0
    self.apply_scale_min = False
    self.scale_min = 0.0
    self.plot_rays_on = False
    self.entry_comparison = os.path.abspath(os.getcwd())
    self.apply_comparison = False
    self.cross_section = 1
    self.show_legend = False
    self.line_label = 'None'
    self.sdf_num2 = 0
    self.plot_all_rays = False
    self.select_ray = False
    self.y_dir_cross_section = False
    self.dir_list = [' '] * 2
    self.num_dir = 2
    self.movie = False

    # Time history params
    self.var_name2 = 'None'
    self.grid_choice = 'None'
    self.cbar_colour_scale = 1.0


def open_var_2d(dat, var_name, parameters):
  """ This routine uses the [var_name] to determine which variable to plot
  from [dat] and find the correct grid to plot it on. It also applies the
  change to polar coordinates (radius vs theta), whether the variable is
  logged, and whether we are doing a mean subtraction from the input controls
  in [parameters].
  """
  var = getattr(dat, var_name)
  var_grid = getattr(var, 'grid')

  grid_conv = getattr(var_grid, 'unit_conversion')
  x_data = getattr(var_grid, 'data')[0] * grid_conv
  y_data = getattr(var_grid, 'data')[1] * grid_conv

  if dat.Logical_flags.use_rz:
    x_label = 'R (' + getattr(var_grid, 'units_new') + ')'
    y_label = 'Z (' + getattr(var_grid, 'units_new') + ')'
  else:
    x_label = 'X (' + getattr(var_grid, 'units_new') + ')'
    y_label = 'Y (' + getattr(var_grid, 'units_new') + ')'

  if parameters.use_polar:
    units = getattr(var_grid, 'units_new')
    x_data, y_data, x_label, y_label = polar_coordinates(x_data, y_data, units)

  c_data = getattr(var, 'data') * getattr(var, 'unit_conversion')
  c_label = getattr(var, "name") + " (" + getattr(var, "units_new") + ")"
  if parameters.view_anisotropies:
    c_data, c_label = mean_subtract(c_data, c_label)

  if parameters.use_log:
    c_data = abs(c_data) + small_num
    c_data = np.log10(c_data)
    c_label = 'log10(' + c_label + ')'

  return x_data, y_data, c_data, x_label, y_label, c_label



def snapshot(fig, ax1, cax1, *args, **kwargs):
  """This function plots [var_name] from the data set [dat] on [ax1].
  """
  parameters = kwargs.get('parameters', plot_parameters())
  data_struct = kwargs.get('data', data_structure())
  if parameters.var_name[1] == "None":
    var_name = parameters.var_name[0]
  else:
    var_name = parameters.var_name[1]
  if parameters.grid_boolean == False:
    grid_colour = 'None'
  else:
    grid_colour = 'k'

  x_data, y_data, c_data, x_label, y_label, c_label = \
      open_var_2d(data_struct.data[0], var_name, parameters)

  if parameters.apply_comparison[1]:
    x_data1, y_data1, c_data1, _, _, _ \
        = open_var_2d(data_struct.data[1], var_name, parameters)

    x_size = max(np.shape(x_data)[0], np.shape(x_data1)[0])
    y_size = np.shape(x_data)[1] + np.shape(x_data1)[1]

    new_x_data = np.zeros((x_size, y_size))
    new_x_data[x_size-np.shape(x_data1)[0]:,:np.shape(x_data1)[1]] \
        = np.flip(-x_data1,1)
    new_x_data[x_size-np.shape(x_data)[0]:,np.shape(x_data1)[1]:] \
        = x_data

    new_y_data = np.zeros((x_size, y_size))
    new_y_data[x_size-np.shape(x_data1)[0]:,:np.shape(x_data1)[1]] \
        = np.flip(y_data1,1)
    new_y_data[x_size-np.shape(x_data)[0]:,np.shape(x_data1)[1]:] \
        = y_data

    # c_data is 1 smaller than x and y as it is cell centred so we minus 1 from
    # x and 2 from y but we need an extra blank space to combine the data sets
    # so we +1 in the y direction
    new_c_data = np.zeros((x_size-1, y_size-1))
    new_c_data[x_size-1-np.shape(c_data1)[0]:,:np.shape(c_data1)[1]] \
        = np.flip(c_data1,1)
    new_c_data[x_size-1-np.shape(c_data)[0]:,np.shape(c_data1)[1]+1:] = c_data

    x_data = new_x_data
    y_data = new_y_data
    c_data = new_c_data

    t_label = time_label(data_struct.data[0], parameters.line_labels[0])
    t_label1 = time_label(data_struct.data[1], parameters.line_labels[1])
    t_label = t_label1 + " and " + t_label
  else:
    t_label = time_label(data_struct.data[0], " ")

  if parameters.reset_axis:
    zoomed_axis1 = np.array([np.min(x_data[:-1,:]), np.max(x_data[:-1,:]),
                             np.min(y_data[:-1,:]), np.max(y_data[:-1,:])])
  else:
    zoomed_axis1 = np.array([ax1.get_xlim()[0], ax1.get_xlim()[1],
                             ax1.get_ylim()[0], ax1.get_ylim()[1]])

  ax1.clear() # This is nessasary for speed
  cax1.clear()

  if parameters.use_log:
    cmin = np.mean(c_data) - 2.0
    cmax = np.max(c_data)
  else:
    cmin = np.min(c_data)
    cmax = np.max(c_data)

  if parameters.apply_scale_max:
    cmax = parameters.scale_max
  if parameters.apply_scale_min:
    cmin = parameters.scale_min

  cmesh = ax1.pcolormesh(x_data, y_data, c_data, linewidth=0.1)
  cmesh.set_edgecolor(grid_colour)
  cmesh.set_clim(cmin, cmax)
  cbar = fig.colorbar(cmesh, cax=cax1)

  if parameters.plot_rays_on:
    wrapper_plot_light_rays(data_struct.data[0], parameters, fig, ax1)
    wrapper_plot_electron_rays(data_struct.data[0], parameters, fig, ax1)

  if not (parameters.surface_name == 'None'):
    cs = parameters.cross_section
    surface = getattr(data_struct.data[0], parameters.surface_name)
    ax1.plot(x_data[surface.index[cs],:],y_data[surface.index[cs],:], 'w:',
        linewidth = 2)

  ax1.set_xlabel(x_label, fontsize = fs)
  ax1.set_ylabel(y_label, fontsize = fs)
  cbar.set_label(c_label, fontsize = fs)
  ax1.set_title(t_label, fontsize = fs)

  ax1.tick_params(axis='x', labelsize = fs)
  ax1.tick_params(axis='y', labelsize = fs)
  cbar.ax.tick_params(labelsize=fs)

  cbar.ax.yaxis.get_offset_text().set(size = fs)
  ax1.xaxis.get_offset_text().set_size(fs)
  ax1.yaxis.get_offset_text().set_size(fs)

  new_xlim = zoomed_axis1[:2]
  ax1.set_xlim(new_xlim)
  new_ylim = zoomed_axis1[2:]
  ax1.set_ylim(new_ylim)
  cbar.draw_all()
  # Remove plot axis, great for gifs
  #ax1.axis('off')

  plt.show()
  plt.pause(0.1)



def wrapper_plot_light_rays(dat, parameters, fig, ax1):
  var = getattr(dat, parameters.var_name[0])
  var_grid = getattr(var, 'grid')
  grid_conv = getattr(var_grid, 'unit_conversion')
  select_ray = parameters.select_ray
  if hasattr(dat, 'Beam1'):
    skip = 1
    plot_rays('Beam1', 'Energy', skip, dat, fig, ax1, parameters.use_polar,
              grid_conv)
  else:
    print(" ")
    print("No light ray data found")



def wrapper_plot_electron_rays(dat, parameters, fig, ax1):
  var = getattr(dat, parameters.var_name[0])
  var_grid = getattr(var, 'grid')
  grid_conv = getattr(var_grid, 'unit_conversion')
  select_ray = parameters.select_ray
  if hasattr(dat, 'Burst1'):
    if parameters.plot_all_rays:
      num_burs = len(dat.bursts)
      for iname in range(0, num_burs):
        skip = 1
        plot_rays(dat.bursts[iname], 'Energy_Deposited', skip, dat, fig, ax1,
                  parameters.use_polar, grid_conv)
    else:
      iname = parameters.select_ray
      skip = 1
      plot_rays(dat.bursts[iname], 'Energy_Deposited', skip, dat, fig, ax1,
                parameters.use_polar, grid_conv)
  else:
    print(" ")
    print("No electron path data found")



def time_label(dat, data_name):
  time = getattr(dat, "Times")
  t_data = getattr(time, "data") * getattr(time, 'unit_conversion')
  t_label = data_name + ' ' + getattr(time, "name") \
          + ' = {0:5.3f}'.format(t_data) + getattr(time, "units_new")
  return t_label



def mean_subtract(cc, cl):
  """Mean subtracts the input data [cc] and updates the label [cl].
  """
  c_data = (cc - np.mean(cc, 1, keepdims = True)) \
         / np.maximum(np.mean(cc, 1, keepdims = True), 1e-17)
  c_label = cl + "[As fraction of average]"
  return c_data, c_label



def polar_coordinates(xc, yc, units):
  """Chnge to polar coordinates, radius vs theta.
  """
  x_label = 'Distance from origin (' + units + ')'
  y_label = "Angle from X-axis (radians)"

  x_data = np.sqrt(xc**2 + yc**2)
  y_data = np.arctan2(yc, xc)
  # This correction makes singularity at radius = 0, look more intuitive.
  y_data[0,:] = y_data[1,:]

  return x_data, y_data, x_label, y_label



def mass(*args, **kwargs):
  """Self-sufficient routine that determines the mass of each cell in the
  starting grid and plots it.
  """
  dat=sh.getdata(0, verbose=False)

  fac = 1.0
  if dat.Logical_flags.use_rz:
    fac = 2*np.pi

  vol=dat.Fluid_Volume.data * fac
  mass=rho[:,:]*vol[:,:]

  # The linout is abitrarily taken at halfway through the domain.
  half = round(np.shape(dat.Fluid_Volume.data)[1] / 2)
  cross_section = kwargs.get('cross_section', half)

  print("Total mass is: ", np.sum(np.sum(mass)))

  X=x[:,:]
  Y=y[:,:]

  fig1=plt.figure()
  plt.pcolormesh(X,Y,mass,edgecolor='none')
  cbar = plt.colorbar()
  cbar.set_label('Mass (kg)')
  #plt.gca().set_aspect('equal', adjustable='box')

  fig2=plt.figure()
  plt.plot(xc[:,cross_section],mass[:,cross_section])
  plt.plot(xc[:,cross_section],mass[:,cross_section],'*')
  plt.xlabel('Radius (m)')
  plt.ylabel('Mass (kg)')

  plt.show()



def empty_lineout(fig, ax):
  """Initilise empty 1D plot and lines to be populated with data later.
  """
  ax1 = ax.twinx()
  ax1.tick_params(axis='y', labelcolor = 'tab:red')

  line_style = ('-','--','-.',':') * 10

  for il in range(0, 10):
    ax_line, = ax.plot(1, lw = 2.5, color='black', linestyle = line_style[il])
    ax_line.set_visible(False)
    line_name = 'line' + str(il + 1)
    ax1_line, = ax1.plot(1, lw = 2, color = 'tab:red',
        linestyle = line_style[il])
    ax1_line.set_visible(False)

  ax_surf = ax.axvline(-big_num, lw = 1, color = 'tab:blue', linestyle = '--')
  setattr(ax, 'surf_tracker', ax_surf)
  ax_surf.set_visible(False)

  return ax1



def lineout(fig, ax, ax1, *args, **kwargs):
  """1D plot of [var_name] from data set [dat] on axis [ax1]. Axis [ax] is
  currently resevred for default variable which in Odin is Fluid_Rho. The [cs]
  provides information about which slice through the data to take. The x and y
  axis maintain scale as the time slider is moved on [ax] but on [ax1] the y
  axis updates with the data.
  """
  parameters = kwargs.get('parameters', plot_parameters())
  data_struct = kwargs.get('data', data_structure())
  var_default = parameters.var_name[0]
  var_name = parameters.var_name[1]
  cs = parameters.cross_section

  if parameters.grid_boolean == False:
    grid_style = 'None'
  else:
    grid_style = 'x'

  try:
    line_colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
  except:
    line_colours = plt.rcParams['axes.color_cycle']

  for num in range(0, parameters.num_dir):
    ax.lines[num].set_label(' ')
    ax1.lines[num].set_label(' ')

    ax.lines[num].set_marker(grid_style)
    ax1.lines[num].set_marker(grid_style)
    if var_default == 'None':
      ax.lines[num].set_visible(False)
    else:
      ax.lines[num].set_visible(parameters.apply_comparison[num])
      ax.lines[num].set_label(parameters.line_labels[num])

    if var_name == 'None':
      ax1.lines[num].set_visible(False)
      ax1.set_visible(False)
      ax.lines[num].set_color(line_colours[num])
    else:
      ax1.lines[num].set_visible(parameters.apply_comparison[num])
      ax1.set_visible(True)
      ax.lines[num].set_color('black')

    if parameters.apply_comparison[num]:
      x_data_alt, y_data_alt, x_label_alt, y_label_alt \
                                = open_var_1d(data_struct.data[num], var_default,
                                              parameters.cross_section,
                                              parameters.use_log,
                                              parameters.y_dir_cross_section,
                                              parameters.use_polar)
      _, y_data1_alt, _, y_label1_alt = open_var_1d(data_struct.data[num], var_name,
                                                    parameters.cross_section,
                                                    parameters.use_log,
                                                    parameters.y_dir_cross_section,
                                                    parameters.use_polar)

      ax.lines[num].set_xdata(x_data_alt)
      ax.lines[num].set_ydata(y_data_alt)
      ax1.lines[num].set_xdata(x_data_alt)
      ax1.lines[num].set_ydata(y_data1_alt)

    if num == 0:
      x_data = x_data_alt
      y_data = y_data_alt
      x_label = x_label_alt
      y_label = y_label_alt
      y_data1 = y_data1_alt
      y_label1 = y_label1_alt

  if parameters.apply_comparison[1]:
    t_label = time_label(data_struct.data[0], parameters.line_labels[0])
    t_label1 = time_label(data_struct.data[1], parameters.line_labels[1])
    t_label = t_label + " and " + t_label1
  else:
    ax.lines[1].set_xdata(1)
    ax.lines[1].set_ydata(1)
    ax1.lines[1].set_xdata(1)
    ax1.lines[1].set_ydata(1)

    t_label = time_label(data_struct.data[0], " ")

  ax.set_xlabel(x_label, fontsize = fs)
  ax.set_ylabel(y_label, fontsize = fs)
  ax1.set_ylabel(y_label1, color='tab:red', fontsize = fs)
  ax.set_title(t_label, fontsize = fs)

  ax.xaxis.get_offset_text().set_size(fs)
  ax.yaxis.get_offset_text().set_size(fs)

  if parameters.use_log:
    ymin = np.mean(y_data[:-1]) - 2.0
    ymax = np.max(y_data[:-1]) + 0.3
    ymin1 = np.mean(y_data1[:-1]) - 2.0
    ymax1 =  np.max(y_data1[:-1]) + 0.3
  else:
    ymin = np.min(y_data[:-1])
    ymax = 1.3 * np.max(y_data[:-1])
    ymin1 = np.min(y_data1[:-1])
    ymax1 = 1.3 * np.max(y_data1[:-1])

  # This section updates the axis with information from the data if
  # [reset_axis] is true or from the prevouis plot if false.
  if parameters.reset_axis:
    zoomed_axis = np.array([np.min(x_data[:-1]), np.max(x_data[:-1]),
                            ymin, ymax])
    zoomed_axis1 = np.array([np.min(x_data[:-1]), np.max(x_data[:-1]),
                             ymin1, ymax1])
  else:
    zoomed_axis = np.array([ax.get_xlim()[0], ax.get_xlim()[1],
                             ax.get_ylim()[0], ax.get_ylim()[1]])
    zoomed_axis1 = np.array([ax1.get_xlim()[0], ax1.get_xlim()[1],
                             ymin1, ymax1])

  ax_surf = getattr(ax, 'surf_tracker')
  # Track a particular point in the data as time is updated
  if (parameters.surface_name == 'None' or parameters.y_dir_cross_section):
    ax_surf.set_visible(False)
    surface_location = 'None'
    surface_move = 0.0
  else:
    ax_surf.set_visible(True)
    old_surface_location = getattr(ax, "loc_cell_track")
    surface = getattr(data_struct.data[0], parameters.surface_name)
    surface_location = x_data[surface.index[cs]]
    if old_surface_location == 'None':
      surface_move = 0.0
    else:
      surface_move = surface_location - old_surface_location
    ax_surf.set_xdata(surface_location)

  setattr(ax, "loc_cell_track", surface_location)

  ax.tick_params(axis='x', labelsize = fs)
  ax.tick_params(axis='y', labelsize = fs)
  ax1.tick_params(axis='y', labelsize = fs)
  ax1.yaxis.get_offset_text().set_size(fs)

  ax.set_xlim(zoomed_axis[:2] + surface_move)
  ax.set_ylim(zoomed_axis[2:])
  ax1.set_xlim(zoomed_axis1[:2] + surface_move)
  ax1.set_ylim(zoomed_axis1[2:])

  if parameters.show_legend:
    lines, labels = ax.get_legend_handles_labels()
    lines1, labels1 = ax1.get_legend_handles_labels()
    all_lines = lines[0:parameters.num_dir]
    all_labels = labels[0:parameters.num_dir]
    ax.legend(all_lines, all_labels, loc = 'upper right', fontsize = fs-5)
  else:
    try:
      ax.get_legend().remove()
    except AttributeError:
      pass

  plt.show()
  plt.pause(0.1)



def open_var_1d(dat, var_name, cs, use_log, y_dir_cross_section, use_polar):
  """Aligns variable [var_name] from [dat] with the appropriate grid for
  plotting. Use [cs] to find which slice through the data is taken.
  """
  var = getattr(dat, var_name)
  unit_conv = getattr(var, "unit_conversion")
  units = getattr(var, "units_new")
  name = getattr(var, "name")
  grid = getattr(var, "grid")
  grid_name = getattr(grid, "name")
  grid_data = getattr(grid, "data")
  grid_conv = getattr(grid, "unit_conversion")
  grid_units = getattr(grid, "units_new")

  if not y_dir_cross_section:
    pos1 = dat.Grid_Grid_mid.data[0][:,cs] * dat.Grid_Grid_mid.unit_conversion
    pos2 = dat.Grid_Grid_mid.data[1][:,cs] * dat.Grid_Grid_mid.unit_conversion
    y_data = getattr(var, "data")[:,cs] * unit_conv
    if use_polar:
      x_data = np.sqrt(pos1**2 + pos2**2)
      x_label = "Radius " + " (" + grid_units + ")"
      y_label = name + " (" + units + ")"
    else:
      x_data = pos1
      x_label = "Distance in x" + " (" + grid_units + ")"
      y_label = name + " (" + units + ")"

    y_data = one_dim_grid(np.array(grid_data)[:,:,cs], grid_conv, x_data, y_data)
  else:
    pos1 = dat.Grid_Grid_mid.data[0][cs,:] * dat.Grid_Grid_mid.unit_conversion
    pos2 = dat.Grid_Grid_mid.data[1][cs,:] * dat.Grid_Grid_mid.unit_conversion
    y_data = getattr(var, "data")[cs,:] * unit_conv
    if use_polar:
      x_data = np.arctan(pos2 / pos1)
      x_label = "Angle from X-axis" + " (radians)"
      y_label = name + " (" + units + ")"
    else:
      x_data = pos2
      x_label = "Distance in Y" + " (" + grid_units + ")"
      y_label = name + " (" + units + ")"

    y_data = one_dim_grid(np.array(grid_data)[:,cs,:], grid_conv, x_data, y_data)

  if use_log:
    y_data = abs(y_data) + small_num
    y_data = np.log10(y_data)
    y_label = 'log10(' + y_label + ')'

  return x_data, y_data, x_label, y_label



def one_dim_grid(grid, grid_conv, x_data, y_data):
  """This code is used to find the correct grid for a cross section of 2D data
  as a midpoint variable can no longer be plotted against corners as in 2D but
  must be plot against midpoints.
  """
  edge = np.sqrt(grid[0,:]**2 + grid[1,:]**2) * grid_conv
  XP = (edge[:-1] + edge[1:]) * 0.5

  if np.shape(y_data) != np.shape(x_data):
    if np.shape(y_data) == np.shape(edge):
      XP = edge
    elif np.shape(y_data) == np.shape(XP):
      XP = XP
    else:
      print("Unknown geometry variable")
    print("Warning: Linear Interpolation!")
    y_data = np.interp(x_data, XP, y_data)

  return y_data



def main():
  """
  """



if __name__ == "__main__":
  main()
