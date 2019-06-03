# Duncan Barlow, Odin project, Warwick University, 01/19
import sdf_helper as sh
import matplotlib.pyplot as plt
import numpy as np
import glob
import import_sdf as isdf
import sys, os
from matplotlib.widgets import Slider, RadioButtons



def time_history(dat, fig, ax1, cs, *args, **kwargs):
  
  var_name = kwargs.get('var_name', "Fluid_Rho")
  
  ax1.clear() # This is nessasary for speed
  
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
  y_data, c_data = two_dim_grid(dat, c_data, cs)
  x_data, grid = np.meshgrid(dat.Times.all_time_data, y_data[0,:], indexing='ij')

  cmesh = ax1.pcolormesh(x_data, y_data, c_data, linewidth=0.1)
  cbar = getattr(ax1, 'cbar')
  if cbar == 'None':
    cbar = fig.colorbar(cmesh)
    setattr(ax1, 'cbar', cbar)
  cbar.set_clim(np.min(c_data), np.max(c_data))
  cbar.draw_all()
  
  plt.show()


def two_dim_grid(dat, data, cs):
  
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
    ax.set_ylabel(name + ' (' + units + ')')
  
    var = getattr(dat, var_name)
    unit_conv = getattr(var, "unit_conversion")
    units = getattr(var, "units_new")
    name = getattr(var, "name")
    y_data1 = getattr(var, "all_time_data")
  
    ax.set_xlabel(dat.Times.name + " (" + dat.Times.units_new + ")")
    
    x_data = dat.Times.all_time_data
    l1.set_xdata(x_data)
    l2.set_xdata(x_data)
    l3.set_xdata(x_data)
  
    ax1.set_ylabel(name + " (" + units + ")", color='tab:red')
    l3.set_ydata(y_data1)

    ax.set_xlim(np.min(x_data[:-1]), np.max(x_data[:-1]))
    ax.set_ylim(np.min(y_data[:-1]), 1.3 * np.max(y_data[:-1]))
    ax1.set_xlim(np.min(x_data[:-1]), np.max(x_data[:-1]))
    ax1.set_ylim(np.min(y_data1[:-1]), 1.3 * np.max(y_data1[:-1]))
    

    plt.show()

def adiabat(*args, **kwargs):
        """ The adiabat is plot on a time vs mass coordinate graph. The adiabat
        was derived as shown in Atzeni "The Physics of Inertial Fusion" under the
        name of Isentrope. The equation used for degenerate pressure is only
        accurate for DT but is applied to CH as well.
        """
        pathname = os.path.abspath(os.getcwd())
        print('From directory:')
        print(os.path.abspath(pathname))
        runs = glob.glob1(pathname,"*.sdf")
        RunCounter = len(runs)
        print('We have imported ', RunCounter, ' data sets.')
        
        SDFName=pathname+'/'+str(0).zfill(4)+'.sdf'
        dat = sh.getdata(SDFName,verbose=False)
        
        # Volume must be times by 2*pi in RZ
        fac = 1.0
        if dat.Logical_flags.use_rz:
                fac = 2*np.pi
        
        # If a cross section is not chosen by the user then it is taken through
        # the middle
        half = round(np.shape(dat.Fluid_Volume.data)[1] / 2)
        cross_section = kwargs.get('cross_section', half)
        
        # initialise variables
        N = len(xc)
        print('The cross-section has width ', N, ' cells')
        all_time = np.zeros((RunCounter,N))
        all_mass = np.zeros((RunCounter,N))
        all_isentrope = np.zeros((RunCounter,N))
        all_pressure = np.zeros((RunCounter,N))
        all_pressure_ls = np.zeros((RunCounter,N-1))
        
        # loop over time to create an array of the slices
        for n in range(0,RunCounter):
                SDFName=pathname+'/'+str(n).zfill(4)+'.sdf'
                dat = sh.getdata(SDFName,verbose=False)
                
                volume = dat.Fluid_Volume.data[:,cross_section] * fac
                density = dat.Fluid_Rho.data[:,cross_section]
                pressure = dat.Fluid_Pressure_ion.data[:,cross_section]

                mass = density * volume
                # The conversion to electron degeneracy pressure is only true for DT
                deg_pressure = 2.17e12 * (density / 1000)**(5.0 / 3.0) / 10
                
                # As used by Craxton et al 2015 the inverse pressure scale length
                # makes the discontinous shock clear. Requires similar spatial and
                # temporal resolution
                dx = xc[1:,cross_section] - xc[:-1,cross_section]
                dlnp = np.abs(np.log(pressure[1:]) - np.log(pressure[:-1]))
                pressure_ls = dlnp / dx
                
                # Save the cross sections, isentrope's maximum is set as 3
                all_time[n,:] =  t
                all_mass[n,:] =  np.cumsum(mass)
                all_isentrope[n,:] = pressure / deg_pressure#np.minimum(pressure/deg_pressure,3)
                all_pressure[n,:] = pressure
                all_pressure_ls[n,:] = pressure_ls

        x=all_time[:,1]
        y=all_mass[1,:]
        
        X,Y = np.meshgrid(x,y)
        
        plt.figure()
        plt.pcolormesh(X,Y,(np.transpose(all_isentrope)), vmin=0, vmax=4)
        plt.xlabel('Time (s)')
        plt.ylabel('Mass Coordinates (kg)')
        cbar = plt.colorbar()
        cbar.set_label('Adiabat/Isentrope')
        
        plt.figure()
        plt.pcolormesh(X,Y,np.transpose(all_pressure_ls), vmin=10000, vmax=100000)
        plt.xlabel('Time (s)')
        plt.ylabel('Mass Coordinates (kg)')
        cbar = plt.colorbar()
        cbar.set_label('|d(ln(p))/dr|')
        
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
  fs = 10
  
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
        
  x_data = getattr(var_grid, 'data')[0]
  y_data = getattr(var_grid, 'data')[1]
  x_label = 'R (' + getattr(var_grid, 'units')[0] + ')'
  y_label = 'Z (' + getattr(var_grid, 'units')[1] + ')'
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
  ax1.set_xlabel(x_label, fontsize = fs)
  ax1.set_ylabel(y_label, fontsize = fs)
  cbar.set_label(c_label, fontsize = fs)
  
  time = getattr(dat, "Times")
  t_data = getattr(time, "data") * getattr(time, 'unit_conversion')
  t_label = getattr(time, "name") + ' = {0:5.3f}'.format(t_data) + getattr(time, "units_new")
  ax1.tick_params(axis='x', labelsize = fs)
  ax1.tick_params(axis='y', labelsize = fs)
  ax1.set_title(t_label)
  
  ax1.set_xlim(zoomed_axis1[:2])
  ax1.set_ylim(zoomed_axis1[2:])
  
  cbar.set_clim(np.min(c_data), np.max(c_data))
  cbar.ax.tick_params(labelsize=fs)
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



def total_energy(output_number):
  """ Internal, kinetic and pv energy
  """
  dat=sh.getdata(output_number, verbose=False)

  fac = 1.0
  if dat.Logical_flags.use_rz:
          fac = 2*np.pi

  vol = dat.Fluid_Volume.data *   fac
  mass = rho[:,:] * vol[:,:]
  
  print("Total mass is: ", np.sum(np.sum(mass)))

  ei = dat.Fluid_Energy_ion.data
  ee = dat.Fluid_Energy_electron.data
  
  Ei = ei * mass
  Ee = ee * mass
  
  tot_Ei = np.sum(np.sum(Ei))
  tot_Ee = np.sum(np.sum(Ee))
  
  tot_IE = tot_Ei + tot_Ee
  print("Total internal energy in the system is: ", tot_IE)
  
  # pV energy = internal
  pressure = dat.Fluid_Pressure.data
  
  gamma = dat.material_real_flags.gamma
  tot_pV = np.sum(np.sum((pressure * vol) / (gamma - 1)))
  print("Total pV energy in the system is:       ", tot_pV)
  
  # Kinetic energy
  corner_mass = dat.Test_Corner_Mass.data *       fac
  all_mass = corner_mass[0:-1,0:-1] + corner_mass[1:,0:-1] + corner_mass[0:-1,1:] + corner_mass[1:,1:]
  node_mass = all_mass[::2,::2]
  
  if dat.Logical_flags.use_rz:
          v1 = dat.Velocity_VTheta.data
          v2 = dat.Velocity_Vr.data
          v3 = dat.Velocity_Vz.data
  else:
          v1 = dat.Velocity_Vx.data
          v2 = dat.Velocity_Vy.data
          v3 = dat.Velocity_Vz.data
  
  vel_sqr = v1**2 + v2**2 + v3**2
  
  #KE = 0.5 * node_mass * vel_sqr
  KE = 0.5 * corner_mass[::2,::2] * vel_sqr[:-1,:-1]\
     + 0.5 * corner_mass[1::2,1::2] * vel_sqr[1:,1:]\
     + 0.5 * corner_mass[1::2,::2] * vel_sqr[1:,:-1]\
     + 0.5 * corner_mass[::2,1::2] * vel_sqr[:-1,1:]
  tot_KE = np.sum(np.sum(KE))
  print("Total kinetic energy in the system is:  ", tot_KE)
  print("Total energy in the system is:          ", tot_KE + tot_IE)
  
  try:
          laser_dep = dat.Fluid_Energy_deposited_laser.data * mass
          tot_laser = np.sum(np.sum(laser_dep))
          print("Total laser energy in the system is:    ", tot_laser)
  except:
          print("Laser energy requires dump_mask(38)")



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
  ax.set_ylabel(name + ' (' + units + ')')
  
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
  
  ax.set_xlabel(grid_name + " (" + grid_units + ")")
  l1.set_xdata(x_data)
  l2.set_xdata(x_data)
  l3.set_xdata(x_data)
  
  ax1.set_ylabel(name + " (" + units + ")", color='tab:red')
  l3.set_ydata(y_data1)
  l3.set_label("Not Fixed")
  l1.set_marker(grid_style)
  
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
      * dat.Times.unit_conversion))

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
    

class axis_data:
  def __init__(self):
    self.new_axis = np.zeros(4)
    self.old_axis = np.zeros(4)
    self.zoomed_axis1 = np.zeros(4)
    self.zoomed_axis2 = np.zeros(4)



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
