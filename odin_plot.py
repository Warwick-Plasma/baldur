# Duncan Barlow, Odin project, Warwick University, 01/19
import sdf_helper as sh
import matplotlib.pyplot as plt
import numpy as np
import glob
import import_sdf as isdf
import sys, os
from matplotlib.widgets import Slider, RadioButtons

# This sets a global fontsize
global fs
fs = 12



def time_history(dat, fig, ax1, *args, **kwargs):
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
  x_data, y_data1 = np.meshgrid(dat.Times.all_time_data, y_data[0,:], indexing='ij')
  
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

  cbar_range = np.max(c_data) - np.min(c_data)
  cbar_max = np.min(c_data) + np.exp(np.log(cbar_range) + cbar_upscale)
  
  if reset_axis:
    zoomed_axis1 = np.array([np.min(x_data[:-1,:]), np.max(x_data[:-1,:]), 
                             np.min(y_data[:-1,:]), np.max(y_data[:-1,:])])
  else:
    zoomed_axis1 = np.array([ax1.get_xlim()[0], ax1.get_xlim()[1], 
                             ax1.get_ylim()[0], ax1.get_ylim()[1]])
  ax1.clear() # This is nessasary for speed
  
  cmesh = ax1.pcolormesh(x_data, y_data, c_data, linewidth=0.1)
  
  cbar = getattr(ax1, 'cbar')
  if cbar == 'None':
    cbar = fig.colorbar(cmesh)
    setattr(ax1, 'cbar', cbar)
  cmesh.set_clim(np.min(c_data), cbar_max)
  cbar.set_clim(np.min(c_data), cbar_max)
  
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
  l2 = getattr(ax1, 'line1')
  l3 = getattr(ax1, 'line2')

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
    y_data1 = getattr(var, "all_time_data")

    x_data = dat.Times.all_time_data
    l1.set_xdata(x_data)
    l2.set_xdata(x_data)
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



def snapshot(dat, fig, ax1, *args, **kwargs):
  """
  """
  
  var_name = kwargs.get('var_name', "Fluid_Rho")
  grid_boolean = kwargs.get('grid_boolean', False)
  if grid_boolean == False:
    grid_colour = 'None'
  else:
    grid_colour = 'k'
  use_polar = kwargs.get('use_polar', False)
  reset_axis = kwargs.get('reset_axis', True)
  view_anisotropies = kwargs.get('view_anisotropies', False)
  
  var = getattr(dat, var_name)
  var_grid = getattr(var, 'grid')
        
  x_data = getattr(var_grid, 'data')[0] * getattr(var_grid, 'unit_conversion')
  y_data = getattr(var_grid, 'data')[1] * getattr(var_grid, 'unit_conversion')
  x_label = 'R (' + getattr(var_grid, 'units_new') + ')'
  y_label = 'Z (' + getattr(var_grid, 'units_new') + ')'
  if use_polar: x_data, y_data, y_label = polar_coordinates(x_data, y_data)
  
  c_data = getattr(var, 'data') * getattr(var, 'unit_conversion')
  c_label = getattr(var, "name") + " (" + getattr(var, "units_new") + ")"
  if view_anisotropies: c_data, c_label = mean_subtract(c_data, c_label)
  
  if reset_axis:
    zoomed_axis1 = np.array([np.min(x_data[:-1,:]), np.max(x_data[:-1,:]), 
                             np.min(y_data[:-1,:]), np.max(y_data[:-1,:])])
  else:
    zoomed_axis1 = np.array([ax1.get_xlim()[0], ax1.get_xlim()[1], 
                             ax1.get_ylim()[0], ax1.get_ylim()[1]])
  ax1.clear() # This is nessasary for speed
  
  cbar = getattr(ax1, 'cbar')
  cmesh = ax1.pcolormesh(x_data, y_data, c_data, linewidth=0.1)
  cmesh.set_edgecolor(grid_colour)
  if cbar == 'None':
    cbar = fig.colorbar(cmesh)
    setattr(ax1, 'cbar', cbar)
  cmesh.set_clim(np.min(c_data), np.max(c_data))
  cbar.set_clim(np.min(c_data), np.max(c_data))
  
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
  
  ax1.set_xlim(zoomed_axis1[:2])
  ax1.set_ylim(zoomed_axis1[2:])
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
  
  l1, = ax.plot(1, lw = 2.5, color='black')
  setattr(ax, 'line1', l1)
  
  l2, = ax1.plot(0, lw = 2, color = 'tab:blue')
  setattr(ax1, 'line1', l2)
  l3, = ax1.plot(0, lw = 2, color = 'tab:red')
  setattr(ax1, 'line2', l3)
  ax1.tick_params(axis='y', labelcolor = 'tab:red')
  return ax1



def lineout(dat, cs, fig, ax, ax1, var_name, *args, **kwargs):
  """
  """
  reset_axis = kwargs.get('reset_axis', True)
  grid_boolean = kwargs.get('grid_boolean', False)
  if grid_boolean == False:
    grid_style = 'None'
  else:
    grid_style = 'x'
  
  l1 = getattr(ax, 'line1')
  l2 = getattr(ax1, 'line1')
  l3 = getattr(ax1, 'line2')
  
  var = dat.Fluid_Rho
  y_data = var.data[:,cs] * var.unit_conversion
  name = var.name
  units = var.units_new
  
  l1.set_ydata(y_data)
  y_label = name + ' (' + units + ')'
  ax.set_ylabel(y_label, fontsize = fs)
    
  ax.xaxis.get_offset_text().set_size(fs)
  ax.yaxis.get_offset_text().set_size(fs)
  
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
  y_data1 = getattr(var, "data")[:,cs] * unit_conv
  
  y_data1 = one_dim_grid(np.array(grid_data)[:,:,cs], grid_conv, x_data, y_data1)
  
  l1.set_xdata(x_data)
  l2.set_xdata(x_data)
  l3.set_xdata(x_data)
  l3.set_ydata(y_data1)
  
  l1.set_marker(grid_style)
  
  x_label = grid_name + " (" + grid_units + ")"
  y_label = name + " (" + units + ")"

  ax.set_xlabel(x_label, fontsize = fs)
  ax1.set_ylabel(y_label, color='tab:red', fontsize = fs)

  ax.tick_params(axis='x', labelsize = fs)
  ax.tick_params(axis='y', labelsize = fs)
  ax1.tick_params(axis='y', labelsize = fs)
  ax1.yaxis.get_offset_text().set_size(fs)
  
  if reset_axis:
    zoomed_axis = np.array([np.min(x_data[:-1]), np.max(x_data[:-1]), 
                             np.min(y_data[:-1]), 1.3 * np.max(y_data[:-1])])
    zoomed_axis1 = np.array([np.min(x_data[:-1]), np.max(x_data[:-1]), 
                             np.min(y_data1[:-1]), 1.3 * np.max(y_data1[:-1])])
  else:
    zoomed_axis = np.array([ax.get_xlim()[0], ax.get_xlim()[1], 
                             ax.get_ylim()[0], ax.get_ylim()[1]])
    zoomed_axis1 = np.array([ax1.get_xlim()[0], ax1.get_xlim()[1], 
                             np.min(y_data1[:-1]), 1.3 * np.max(y_data1[:-1])])
  ax.set_xlim(zoomed_axis[:2])
  ax.set_ylim(zoomed_axis[2:])
  ax1.set_xlim(zoomed_axis1[:2])
  ax1.set_ylim(zoomed_axis1[2:])
  
  ax.set_title(dat.Times.name
      + ' = {0:5.3f}'.format(dat.Times.data
      * dat.Times.unit_conversion), fontsize = fs)

  plt.show()


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
    print("Waring: Linear Interpolation!")
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
