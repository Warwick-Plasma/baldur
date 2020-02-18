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
global small_num
global big_num
fs = 15
small_num = 1e-100
big_num = 1e100



def plot_laser_profile(name):
  with open(name) as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    times = []
    powers = []
    for row in readCSV:
      time = float(row[0]) / 1000.0
      times.append(time)
      power = float(row[1])
      powers.append(power)
      print(time, power)
  plt.figure()
  plt.xlabel('Time (ns)')
  plt.ylabel('Power (TW)')
  plt.plot(times, powers)
  plt.show()



def time_history(dat, fig, ax1, cax1, *args, **kwargs):
  """A pcolormesh plot of space against time with a variable shown in colour
  """
  
  var_name = kwargs.get('var_name', "Fluid_Rho")
  cbar_upscale = kwargs.get('cbar_upscale', -10.0)
  reset_axis = kwargs.get('reset_axis', True)
  grid_choice = kwargs.get('grid', 'default')
  
  var = getattr(dat, var_name)
  unit_conv = getattr(var, "unit_conversion")
  units = getattr(var, "units_new")
  name = getattr(var, "name")
  grid = getattr(var, "grid")
  grid_name = getattr(grid, "name")
  grid_data = getattr(grid, "all_time_data")
  grid_conv = getattr(grid, "unit_conversion")
  grid_units = getattr(grid, "units_new")
  
  c_data = getattr(var, "all_time_data") * unit_conv
  y_data, c_data = two_dim_grid(dat, c_data)
  x_data, y_data1 = np.meshgrid(dat.Times.all_time_data * dat.Times.unit_conversion, y_data[0,:], indexing='ij')
  
  x_label = dat.Times.name + ' (' + getattr(dat.Times, 'units_new') + ')'
  y_label = 'Radius (' + grid_units + ')'
  c_label = name + " (" + units + ")"
  
  if (grid_choice == 'default'):
    y_data = y_data
  elif (grid_choice == 'initial'):
    y_data = y_data1
    y_label = 'Initial ' + y_label
  elif (grid_choice == 'cell number'):
    pos = np.linspace(0, np.shape(c_data)[1]-1, np.shape(c_data)[1])
    _, y_data = np.meshgrid(dat.Times.all_time_data, pos, indexing='ij')
    y_label = 'Cell Number'

  cbar_range = np.max(c_data) - np.min(c_data) + small_num
  cbar_max = np.min(c_data) + np.exp(np.log(cbar_range) + cbar_upscale)
  
  if reset_axis:
    zoomed_axis1 = np.array([np.min(x_data[:-1,:]), np.max(x_data[:-1,:]), 
                             np.min(y_data[:-1,:]), np.max(y_data[:-1,:])])
  else:
    zoomed_axis1 = np.array([ax1.get_xlim()[0], ax1.get_xlim()[1], 
                             ax1.get_ylim()[0], ax1.get_ylim()[1]])
  
  ax1.clear() # This is nessasary for speed
  cax1.clear()
  
  cmesh = ax1.pcolormesh(x_data, y_data, c_data, linewidth=0.1)
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



def time_history_lineout(dat, fig, ax, ax1, *args, **kwargs):
  use_analysis = kwargs.get('use_analysis', False)
  var_name = kwargs.get('var_name', "Laser_Energy_Total_Deposited")

  l1 = getattr(ax, 'line1')
  l2 = getattr(ax, 'line2')
  l3 = getattr(ax1, 'line1')

  if use_analysis:
    var = dat.Laser_Power_Total_Deposited
    y_data = var.all_time_data * var.unit_conversion
    name = var.name
    units = var.units_new

    l1.set_ydata(y_data)
    y_label = name + " (" + units + ")"
    ax.set_ylabel(y_label, fontsize = fs)
    
    ax.xaxis.get_offset_text().set_size(fs)
    ax.yaxis.get_offset_text().set_size(fs)

    var = getattr(dat, var_name)
    unit_conv = getattr(var, "unit_conversion")
    units = getattr(var, "units_new")
    name = getattr(var, "name")
    y_data1 = getattr(var, "all_time_data") * unit_conv

    x_data = dat.Times.all_time_data * dat.Times.unit_conversion
    l1.set_xdata(x_data)
    l2.set_xdata(0)
    l3.set_xdata(x_data)
    l3.set_ydata(y_data1)

    x_label = dat.Times.name + " (" + dat.Times.units_new + ")"
    y_label = name + " (" + units + ")"

    ax.set_xlabel(x_label, fontsize = fs)
    ax1.set_ylabel(y_label, color='tab:red', fontsize = fs)
    
    ax.tick_params(axis='x', labelsize = fs)
    ax.tick_params(axis='y', labelsize = fs)
    ax1.tick_params(axis='y', labelsize = fs)
    ax1.yaxis.get_offset_text().set_size(fs)

    ax.set_xlim(np.min(x_data[:-1]), np.max(x_data[:-1]))
    ax.set_ylim(np.min(y_data[:-1]), 1.3 * np.max(y_data[:-1]))
    ax1.set_xlim(np.min(x_data[:-1]), np.max(x_data[:-1]))
    ax1.set_ylim(np.min(y_data1[:-1]), 1.3 * np.max(y_data1[:-1]))

    plt.show()



def check_analysis(use_analysis):
  if use_analysis == True:
    print("starting analysis")
  elif use_analysis == False:
    print("set: <use_analysis = True>, for analysis")
  else:
    print("set: <use_analysis = True> or <False>")
    print("it requires certain dump_masks")
    sys.exit()



def data_and_plot(sdf_num, fig, ax1, cax1, fig2, ax2, ax3, parameters):
  """ This routine is called from Tkinter and calls all the plotting routines.
  The dat file is created with all the data from the sdf file indicated in
  parameters.
  """
  print_string = 'Processing file {:4d}'.format(sdf_num) + ' of {:4d}'.format(parameters.iend) + '   '
  sys.stdout.write('\r' + print_string)
  sys.stdout.flush()
  
  if parameters.apply_comparison:
    if os.path.isdir(parameters.entry_comparison):
      parameters.dat1 = isdf.use_sdf(sdf_num, parameters.entry_comparison, use_analysis = parameters.use_analysis, istart = parameters.istart)
    else:
      parameters.apply_comparison = False
      print()
      print("Warning: " + parameters.entry_comparison + " is not a directory")
   
  
  dat = isdf.use_sdf(sdf_num, parameters.pathname, use_analysis = parameters.use_analysis, istart = parameters.istart)
  
  snapshot(dat, fig, ax1, cax1, parameters.var_name, parameters = parameters)
  
  lineout(dat, parameters.cross_section, fig2, ax2, ax3, parameters.var_name, parameters = parameters)



def plot_colourline(fig1, ax1, x, y, c, cnorm):
  # Create a set of line segments so that we can color them individually
  # This creates the points as a N x 1 x 2 array so that we can stack points
  # together easily to get the segments. The segments array for line collection
  # needs to be (numlines) x (points per line) x 2 (for x and y)
  points = np.array([x, y]).T.reshape(-1, 1, 2)
  segments = np.concatenate([points[:-1], points[1:]], axis=1)

  # Create a continuous norm to map from data points to colors
  lc = LineCollection(segments, cmap='viridis', norm=cnorm)
  # Set the values used for colormapping
  lc.set_array(c)
  lc.set_linewidth(2)
  ax1.add_collection(lc)



def plot_rays(name, name_var, skip, dat, fig1, ax1, use_polar, grid_conv):
  
  beam = getattr(dat, name)
  beam_energy = getattr(dat, name + '_' + name_var)
  nrays = len(beam.data)
  cmax = max(max(beam_energy.data, key=lambda x: max(x.data)).data)
  cmin = min(min(beam_energy.data, key=lambda x: min(x.data)).data)
  cnorm = plt.Normalize(cmin, cmax)
  for iray in range(0, nrays, skip):
    print_string = 'Processing ray {:4d}'.format(iray+1) + ' of {:4d}'.format(nrays) + '   '
    sys.stdout.write('\r' + print_string)
    sys.stdout.flush()
    
    x_ray = beam.data[iray].data[0] * grid_conv
    y_ray = beam.data[iray].data[1] * grid_conv
    c_ray = beam_energy.data[iray].data
    
    if use_polar: x_ray, y_ray, y_label = polar_coordinates(x_ray, y_ray)
    
    plot_colourline(fig1, ax1, x_ray, y_ray, c_ray, cnorm)
  smap = cm.ScalarMappable(norm=cnorm, cmap='viridis')
  smap.set_array([])
  #fig1.colorbar(smap)



class plot_parameters:
  def __init__(self):
    self.sdf_num = 0
    self.use_analysis = False
    self.pathname = 'None'
    self.istart = 0
    self.iend = 0
    self.grid_boolean = False
    self.use_polar = False
    self.var_name = 'None'
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
    self.dat1 = None
    self.cross_section = 1



def open_var_2d(dat, var_name, parameters):
  
  var = getattr(dat, var_name)
  var_grid = getattr(var, 'grid')
  
  grid_conv = getattr(var_grid, 'unit_conversion')
  x_data = getattr(var_grid, 'data')[0] * grid_conv
  y_data = getattr(var_grid, 'data')[1] * grid_conv
  x_label = 'R (' + getattr(var_grid, 'units_new') + ')'
  y_label = 'Z (' + getattr(var_grid, 'units_new') + ')'
  if parameters.use_polar:
    x_data, y_data, y_label = polar_coordinates(x_data, y_data)
  
  c_data = getattr(var, 'data') * getattr(var, 'unit_conversion')
  c_label = getattr(var, "name") + " (" + getattr(var, "units_new") + ")"
  if parameters.view_anisotropies:
    c_data, c_label = mean_subtract(c_data, c_label)
  
  if parameters.use_log:
    c_data = abs(c_data) + small_num
    c_data = np.log10(c_data)
    c_label = 'log10(' + c_label + ')'
  
  return x_data, y_data, c_data, x_label, y_label, c_label



def snapshot(dat, fig, ax1, cax1, var_name, *args, **kwargs):
  """
  """
  
  parameters = kwargs.get('parameters', plot_parameters())
  if parameters.grid_boolean == False:
    grid_colour = 'None'
  else:
    grid_colour = 'k'
  
  x_data, y_data, c_data, x_label, y_label, c_label = open_var_2d(dat, var_name, parameters)
  
  if parameters.apply_comparison:
    x_data1 = np.zeros(np.shape(x_data)) # this might allow plotting of different sized arrays
    y_data1 = np.zeros(np.shape(y_data))
    c_data1 = np.zeros((np.shape(c_data)[0],np.shape(c_data)[1]+1))
    x_data1[0:,0:], y_data1[0:,0:], c_data1[:,:-1], _, _, _ = open_var_2d(parameters.dat1, var_name, parameters)
    
    x_data = np.hstack((x_data, np.flip(x_data1,1)))
    y_data = np.hstack((y_data, np.flip(-y_data1,1)))
    c_data = np.hstack((c_data, np.flip(c_data1,1)))
  
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
    var = getattr(dat, var_name)
    var_grid = getattr(var, 'grid')
    grid_conv = getattr(var_grid, 'unit_conversion')
    if hasattr(dat, 'Beam1'):
      skip = 1
      plot_rays('Beam1', 'Energy', skip, dat, fig, ax1, parameters.use_polar, grid_conv)
    if hasattr(dat, 'Burst1'):
      num_burs = len(dat.bursts)
      for iname in range(0, num_burs):
        skip = 1
        plot_rays(dat.bursts[iname], 'Energy_Deposited', skip, dat, fig, ax1, parameters.use_polar, grid_conv)
  
  ax1.set_xlabel(x_label, fontsize = fs)
  ax1.set_ylabel(y_label, fontsize = fs)
  cbar.set_label(c_label, fontsize = fs)
  
  ax1.tick_params(axis='x', labelsize = fs)
  ax1.tick_params(axis='y', labelsize = fs)
  cbar.ax.tick_params(labelsize=fs)
  
  cbar.ax.yaxis.get_offset_text().set(size = fs)
  ax1.xaxis.get_offset_text().set_size(fs)
  ax1.yaxis.get_offset_text().set_size(fs)
  
  time = getattr(dat, "Times")
  t_data = getattr(time, "data") * getattr(time, 'unit_conversion')
  t_label = getattr(time, "name") + ' = {0:5.3f}'.format(t_data) + getattr(time, "units_new")
  ax1.set_title(t_label, fontsize = fs)

  new_xlim = zoomed_axis1[:2]
  ax1.set_xlim(new_xlim)
  new_ylim = zoomed_axis1[2:]
  ax1.set_ylim(new_ylim)
  cbar.draw_all()
  
  plt.show()



def mean_subtract(cc, cl):
  c_data = (cc - np.mean(cc, 1, keepdims = True)) / np.maximum(np.mean(cc, 1, keepdims = True), 1e-17)
  c_label = cl + "[As fraction of average]"
  return c_data, c_label



def polar_coordinates(xc, yc):
  y_label = "Radians"
  
  x_data = np.sqrt(xc**2 + yc**2)
  y_data = np.arctan2(yc, xc) / np.pi
  y_data[0,:] = y_data[1,:]
  
  return x_data, y_data, y_label



def mass(*args, **kwargs):
  dat=sh.getdata(0, verbose=False)

  fac = 1.0
  if dat.Logical_flags.use_rz:
          fac = 2*np.pi

  vol=dat.Fluid_Volume.data * fac
  mass=rho[:,:]*vol[:,:]
  
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
  
  ax1 = ax.twinx()
  
  ax_l1, = ax.plot(1, lw = 2.5, color='black')
  setattr(ax, 'line1', ax_l1)
  ax_l2, = ax.plot(1, lw = 2.5, color='black', linestyle = '--')
  setattr(ax, 'line2', ax_l2)
  ax_l3 = ax.axvline(0, lw = 1, color = 'tab:blue', linestyle = '--')
  setattr(ax, 'line3', ax_l3)
  
  ax1_l1, = ax1.plot(1, lw = 2, color = 'tab:red')
  setattr(ax1, 'line1', ax1_l1)
  ax1_l2, = ax1.plot(1, lw = 2, color='tab:red', linestyle = '--')
  setattr(ax1, 'line2', ax1_l2)
  ax1.tick_params(axis='y', labelcolor = 'tab:red')
  return ax1



def lineout(dat, cs, fig, ax, ax1, var_name, *args, **kwargs):
  """
  """
  
  parameters = kwargs.get('parameters', plot_parameters())
  if parameters.grid_boolean == False:
    grid_style = 'None'
  else:
    grid_style = 'x'
  
  ax_l1 = getattr(ax, 'line1')
  ax_l2 = getattr(ax, 'line2')
  ax_l3 = getattr(ax, 'line3')
  ax1_l1 = getattr(ax1, 'line1')
  ax1_l2 = getattr(ax1, 'line2')
  
  if (dat.Header['code_name'] == 'Odin2D'):
    var_default = "Fluid_Rho"
  else:
    var_default = dat.variables[0]
    
  x_data, y_data, x_label, y_label = open_var_1d(dat, var_default, cs, parameters.use_log)
  _, y_data1, _, y_label1 = open_var_1d(dat, var_name, cs, parameters.use_log)
  
  ax_l1.set_xdata(x_data)
  ax_l1.set_ydata(y_data)
  ax1_l1.set_xdata(x_data)
  ax1_l1.set_ydata(y_data1)
  
  if parameters.apply_comparison:
    x_data_comp, y_data_comp, _, _ = open_var_1d(parameters.dat1, var_default, cs, parameters.use_log)
    _, y_data1_comp, _, _ = open_var_1d(parameters.dat1, var_name, cs, parameters.use_log)
  
    ax_l2.set_xdata(x_data_comp)
    ax_l2.set_ydata(y_data_comp)
    ax1_l2.set_xdata(x_data_comp)
    ax1_l2.set_ydata(y_data1_comp)
  else:
    ax_l2.set_xdata(1)
    ax_l2.set_ydata(1)
    ax1_l2.set_xdata(1)
    ax1_l2.set_ydata(1)
  
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
  
  if parameters.surface_name == 'None':
    ax_l3.set_xdata(-big_num)
    surface_location = 'None'
    surface_move = 0.0
  else:
    old_surface_location = getattr(ax, "loc_cell_track")
    surface = getattr(dat, parameters.surface_name)
    surface_location = surface.data[cs] * surface.unit_conversion
    if old_surface_location == 'None':
      surface_move = 0.0
    else:
      surface_move = surface_location - old_surface_location
    ax_l3.set_xdata(surface_location)
    
  setattr(ax, "loc_cell_track", surface_location)
  
  ax_l1.set_marker(grid_style)

  ax.set_xlabel(x_label, fontsize = fs)
  ax.set_ylabel(y_label, fontsize = fs)
  ax1.set_ylabel(y_label1, color='tab:red', fontsize = fs)

  ax.tick_params(axis='x', labelsize = fs)
  ax.tick_params(axis='y', labelsize = fs)
  ax1.tick_params(axis='y', labelsize = fs)
  ax1.yaxis.get_offset_text().set_size(fs)
  
  ax.set_xlim(zoomed_axis[:2] + surface_move)
  ax.set_ylim(zoomed_axis[2:])
  ax1.set_xlim(zoomed_axis1[:2] + surface_move)
  ax1.set_ylim(zoomed_axis1[2:])
  
  ax.set_title(dat.Times.name
      + ' = {0:5.3f}'.format(dat.Times.data
      * dat.Times.unit_conversion), fontsize = fs)

  plt.show()



def open_var_1d(dat, var_name, cs, use_log):
  
  var = getattr(dat, var_name)
  unit_conv = getattr(var, "unit_conversion")
  units = getattr(var, "units_new")
  name = getattr(var, "name")
  grid = getattr(var, "grid")
  grid_name = getattr(grid, "name")
  grid_data = getattr(grid, "data")
  grid_conv = getattr(grid, "unit_conversion")
  grid_units = getattr(grid, "units_new")
  
  pos1 = dat.Grid_Grid_mid.data[0][:,cs] * dat.Grid_Grid_mid.unit_conversion
  pos2 = dat.Grid_Grid_mid.data[1][:,cs] * dat.Grid_Grid_mid.unit_conversion
  x_data = np.sqrt(pos1**2 + pos2**2)
  y_data = getattr(var, "data")[:,cs] * unit_conv
  
  y_data = one_dim_grid(np.array(grid_data)[:,:,cs], grid_conv, x_data, y_data)
  
  x_label = grid_name + " (" + grid_units + ")"
  y_label = name + " (" + units + ")"
  
  if use_log:
    y_data = abs(y_data) + small_num
    y_data = np.log10(y_data)
    y_label = 'log10(' + y_label + ')'
  
  return x_data, y_data, x_label, y_label



def one_dim_grid(grid, grid_conv, x_data, y_data):
  
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
  mass(0)
  adiabat()
  which_sdf = 0
  total_energy(which_sdf)
  snapshot(start = 0, var_name = "Fluid_Rho")
  lineout(start = 0, var_name = "Fluid_Rho")



if __name__ == "__main__":
  main()
