# Duncan Barlow, Odin project, Warwick University, 01/19
import sdf_helper as sh
import matplotlib.pyplot as plt
import numpy as np
import glob
import import_sdf as isdf
import sys, os
from matplotlib.widgets import Slider, RadioButtons



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
                deg_pressure = 2.17e12 * (density / 1000)**(5.0/3.0) / 10
                
                # As used by Craxton et al 2015 the inverse pressure scale length
                # makes the discontinous shock clear. Requires similar spatial and
                # temporal resolution
                dx = xc[1:,cross_section] - xc[:-1,cross_section]
                dlnp = np.abs(np.log(pressure[1:]) - np.log(pressure[:-1]))
                pressure_ls = dlnp / dx
                
                # Save the cross sections, isentrope's maximum is set as 3
                all_time[n,:] =  t
                all_mass[n,:] =  np.cumsum(mass)
                all_isentrope[n,:] = pressure / deg_pressure# np.minimum(pressure / deg_pressure, 3)
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



def check_analysis(analysis):
        if analysis == True:
                print("starting analysis")
        elif analysis == False:
                print("set: <analysis = True>, for analysis")
        else:
                print("set: <analysis = True> or <False>")
                print("it requires certain dump_masks")
                sys.exit()



def snapshot(istart, *args, **kwargs):
  """
  """
  var_name = kwargs.get('var_name', "Fluid_Rho")
  analysis = kwargs.get('analysis', False)
  check_analysis(analysis)
    
  pathname = os.path.abspath(os.getcwd())
  runs = glob.glob1(pathname,"*.sdf")
  RunCounter = len(runs)
  
  dat1 = isdf.use_sdf(istart, pathname, use_analysis = analysis, istart = istart)
  
  plt.ion()
  plt.close('all')
       
  fig=plt.figure(num=1,figsize=(12,8),facecolor='white')
  ax1=plt.axes([0.1, 0.1, 0.5, 0.6])
  
  var = getattr(dat1, 'Fluid_Rho')
  var_grid = getattr(var, 'grid')
        
  x_data = getattr(var_grid, 'data')[0]
  y_data = getattr(var_grid, 'data')[1]
  c_data = getattr(var, 'data') * getattr(var, 'unit_conversion')
  
  cmesh = ax1.pcolormesh(x_data, y_data, c_data, linewidth=0.1)
  cbar = plt.colorbar(cmesh)
  ax1.set_xlim([np.min(x_data[:-1,:]),np.max(x_data[:-1,:])])
  ax1.set_ylim([np.min(y_data[:-1,:]),np.max(y_data[:-1,:])])
  
  axcolor='lightgoldenrodyellow' # slider background colour
  axtime=plt.axes([0.1, 0.85, 0.4, 0.03], facecolor=axcolor) # slider size
  stime=Slider(axtime, ' ', istart, RunCounter-1, valinit = istart, valfmt = '%1.0f')
  
  axcolor = 'lightgoldenrodyellow'
  rax = plt.axes([0.65, 0.05, 0.35, 0.9], facecolor=axcolor)
  radio = RadioButtons(rax, (dat1.variables))
  
  grid_colour = "none"
  setattr(fig, "grid_colour", grid_colour)
  
  def change_variable(label):
        max_x=stime.val
        stime.set_val(max_x)
  
  def press(event): # Defining a function for keyboard inputs for o and p
        max_x=stime.val
        if event.key == 'o' and max_x > istart + 0.5:
                stime.set_val(max_x-1) # define the increment
        elif event.key == 'p' and max_x < RunCounter - 1.5:
                stime.set_val(max_x+1)
        elif event.key == 'i':
                grid_colour = getattr(fig, "grid_colour")
                if grid_colour == "none":
                        grid_colour = "black"
                elif grid_colour == "black":
                        grid_colour = "none"
                setattr(fig, "grid_colour", grid_colour)
                stime.set_val(max_x)
        # publish the current view of the top figure
        elif event.key == 't':
                plt.draw()
                filename1 = 'test.pdf'
                extent = ax1.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
                fig.savefig(filename1, bbox_inches=extent.expanded(1.7, 1.25))
  
  def update(val):
    """
    """
    
    x=ax1.get_xlim()
    y=ax1.get_ylim()
    zoomed_axis1=np.array([x[0],x[1],y[0],y[1]])
  
    sdf_num = int(round(stime.val))
    
    var_name = radio.value_selected
    dat = isdf.use_sdf(sdf_num, pathname, use_analysis = analysis, istart = istart)

    var = getattr(dat, var_name)
    var_grid = getattr(var, 'grid')
    
    X_data = getattr(var_grid, 'data')[0]
    Y_data = getattr(var_grid, 'data')[1]
    C_data = getattr(var, 'data') * getattr(var, 'unit_conversion')
    
    X_label = 'X (' + getattr(var_grid, 'units')[0] + ')'
    Y_label = 'Y (' + getattr(var_grid, 'units')[1] + ')'
    C_label = getattr(var, "name") + " (" + getattr(var, "units_new") + ")"
    
    fs = 12
    ax1.clear()
    cmesh = ax1.pcolormesh(X_data, Y_data, C_data, linewidth=0.01)
    ax1.tick_params(axis='x', labelsize = fs)
    ax1.tick_params(axis='y', labelsize = fs)
    ax1.set_xlabel(X_label, fontsize = fs)
    ax1.set_ylabel(Y_label, fontsize = fs)
    grid_colour = getattr(fig, "grid_colour")
    cmesh.set_edgecolor(grid_colour)
    ax1.set_title('Time = {0:5.3f}'.format(dat.Header['time']*1e9))
    cbar.set_clim(np.min(C_data), np.max(C_data))
    cbar.set_label(C_label, fontsize = fs)
    cbar.ax.tick_params(labelsize=fs) 
    cbar.draw_all()
    
    ax1.set_xlim(zoomed_axis1[:2])
    ax1.set_ylim(zoomed_axis1[2:])

    fig.canvas.draw_idle()

  fig.canvas.mpl_connect('key_press_event', press) # this makes the keyboard function work (key press)
  stime.on_changed(update) # this makes the slider work
  radio.on_clicked(change_variable)
  update(istart)
  plt.show()



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



def lineout(istart, *args, **kwargs):
  """
  """
  use_analysis = kwargs.get('analysis', False)
  check_analysis(use_analysis)
  
  pathname = os.path.abspath(os.getcwd())
  dat = isdf.use_sdf(istart, pathname, use_analysis = use_analysis)
  
  nmat = dat.Integer_flags.nmat # assume nmat doesn't change
  runs = glob.glob1(pathname,"*.sdf")
  RunCounter = len(runs)-1
  RunCounter = kwargs.get('end', RunCounter) + 1
  
  len_x = np.shape(dat.Fluid_Rho.data)[0]
  len_y = np.shape(dat.Fluid_Rho.data)[1]
  half = round(len_y / 2)
  cs = kwargs.get('cross_section', half)
  
  dat = isdf.get_data_all(dat, istart, RunCounter, pathname, use_analysis, cs)
  
  plt.ion()
  plt.close('all')
  
  fig, (ax,ax2)=plt.subplots(2, 1, num=1, figsize=(12,10), facecolor='white')
  plt.subplots_adjust(left  = 0.3,  # the left side of the subplots of the figure
                      right = 0.9,    # the right side of the subplots of the figure
                      bottom = 0.15,   # the bottom of the subplots of the figure
                      top = 0.9,      # the top of the subplots of the figure
                      wspace = 0.25,   # the amount of width reserved for space between subplots
                      hspace = 0.3)
  fs = 12
  
  ax1 = ax.twinx()
  
  y_var = dat.Fluid_Rho
  label = y_var.name + " (" + y_var.units + ")"
  ax.set_ylabel(label, fontsize = fs)  
  l1, = ax.plot(1, lw = 2.5, color='black')
  ax.tick_params(axis='y', labelcolor='black', labelsize = fs)
  
  x_var = dat.Grid_Grid_mid
  x_data = x_var.data[0] * x_var.unit_conversion
  label = x_var.name + " (" + x_var.units_new + ")"
  ax.set_xlabel(label, fontsize = fs)
  ax.tick_params(axis='x', labelcolor = 'black', labelsize = fs)
  
  ax1.set_ylabel('Place holder', color='tab:red', fontsize = fs)
  l2, = ax1.plot(0, lw = 2, color = 'tab:blue')
  l3, = ax1.plot(0, lw = 2, color = 'tab:red')
  ax1.tick_params(axis='y', labelcolor = 'tab:red', labelsize = fs)

  xmin = 0 #np.min(all_time.radius[0,:])
  xmax = np.max(x_data)
  ax.set_xlim([xmin, xmax])
  
  ax3 = ax2.twinx()
  
  label = dat.Times.name + ' (' + dat.Times.units_new + ')' 
  ax2.set_xlabel(label, fontsize = fs)
  ax2.tick_params(axis='x', labelcolor = 'black', labelsize = fs)
  x_data = dat.Times.all_time_data * dat.Times.unit_conversion
  
  # Plot laser power deposited
  ax2.set_ylabel('Laser power deposited (W)', fontsize = fs)
  dt = dat.Times.all_time_data[1:] - dat.Times.all_time_data[:-1]
  var = dat.Total_Energy_Laser_deposited.all_time_data
  y_data = (var[1:] - var[:-1]) / dt
  y_data = np.insert(y_data, 0, 0)
  la = ax2.plot(x_data, y_data, lw = 2, color='black')
  ax2.tick_params(axis='y', labelcolor='black', labelsize = fs)
  
  # Plot maximum density
  var = dat.Fluid_Rho
  y_data = np.max(var.all_time_data, axis=1) * var.unit_conversion
  ax3.set_ylabel('Maximum Density' + var.units_new, color='tab:red', fontsize = fs)
  lb = ax3.plot(x_data, y_data, lw = 2.5, color='tab:red')
  ax3.tick_params(axis='y', labelcolor='tab:red', labelsize = fs)

  xmin = 0 #np.min(x_data)
  xmax = np.max(x_data)
  #ymin = 0
  #ymax = 1
  ax2.set_xlim([xmin, xmax])
  #ax2.set_ylim([ymin, ymax])
  
  # Setting up slider
  axcolor='lightgoldenrodyellow'
  axtime=plt.axes([0.3, 0.05, 0.6, 0.03], facecolor=axcolor)
  stime=Slider(axtime, 'SDF selector', istart, RunCounter-1, valinit = istart, valfmt = '%1.0f')
  stime.label.set_size(fs)
  
  axcolor = 'lightgoldenrodyellow'
  rax = plt.axes([0.0, 0.05, 0.18, 0.90], facecolor=axcolor)
  radio = RadioButtons(rax, (dat.variables))
  
  ax_data = axis_data()
  ax1_data = axis_data()
  
  ii = 0
  iv = 0
  # to avoid globals by setting figure attributes
  setattr(fig,'reset_axis',ii)
  setattr(fig,'reset_axis2',iv)
  setattr(fig,'log_on',0)
  

  def change_variable(label):
    max_x=stime.val
    iv=0
    setattr(fig,'iv',iv)
    stime.set_val(max_x)

  def press(event): #keyboard movement of slider, currently set to 'o' and 'p'
    t0 = stime.val
    # return one time step
    if event.key == 'o' and t0 > istart + 0.5:
      stime.set_val(t0-1)
    # advance one time step
    elif event.key == 'p' and t0 < RunCounter - 1.5:
      stime.set_val(t0+1)
    # turn on cell markers
    elif event.key == 'i':
      marker_onoff = l1.get_marker()
      if marker_onoff == 'None':
        l1.set_marker('x')
        #l2.set_marker('x')
        #l3.set_marker('x')
      else:
        l1.set_marker('None')
        #l2.set_marker('None')
        #l3.set_marker('None')
      stime.set_val(t0)
    # reset the axis (for this time step)
    elif event.key == 'u':
      ii = 0
      iv = 0
      setattr(fig,'reset_axis',ii)
      setattr(fig,'reset_axis2',iv)
      stime.set_val(t0)
    # turn on or off logarithmic density
    elif event.key == 'y':
      log_on = getattr(fig,'log_on')
      if log_on == 0:
        setattr(fig,'log_on',1)     
      else:
        setattr(fig,'log_on',0)  
      setattr(fig,'reset_axis',0)
      setattr(fig,'reset_axis2',0)
      stime.set_val(t0)
    # publish the current view of the top figure
    elif event.key == 't':
      plt.draw()
      filename1 = 'test.pdf'
      extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
      fig.savefig(filename1, bbox_inches=extent.expanded(1.3, 1.3))
  
  def update(val):
    '''update(val)
    This updates all the lines when the slider is moved
    '''
    x=ax.get_xlim()
    y=ax.get_ylim()
    ax_data.zoomed_axis1=np.array([x[0],x[1],y[0],y[1]])
    x=ax1.get_xlim()
    y=ax1.get_ylim()
    ax1_data.zoomed_axis1=np.array([x[0],x[1],y[0],y[1]])
    
    ii = getattr(fig,'reset_axis')
    iv = getattr(fig,'reset_axis2')
    t0 = int(np.around(stime.val))
    
    log_on = getattr(fig,'log_on')
    var = dat.Fluid_Rho
    y_data = var.all_time_data[t0,:] * var.unit_conversion
    name = var.name
    units = var.units_new
    
    if log_on == 0:
      l1.set_ydata(y_data)
    
      ylimL = 0.0
      ylimH = 1.5 * np.max(y_data)
      ax.set_ylabel(name + ' (' + units + ')')  
    else:
      l1.set_ydata(np.log10(y_data))
      ax.set_ylabel('Log[' + name + ' (' + units + ')]')   
    
      ylimL = np.log10(0.0001)#0.0
      ylimH = np.log10(10000)#1.5 * np.max(y_data)
    
    var_name = radio.value_selected
    var = getattr(dat, var_name)
    unit_conv = getattr(var, "unit_conversion")
    units = getattr(var, "units_new")
    name = getattr(var, "name")
    grid_name = getattr(var, "grid")
    grid = getattr(grid_name, "all_time_data")
    grid_conv = getattr(grid_name, "unit_conversion")
    
    x_data = dat.Grid_Grid_mid.all_time_data[0,t0,:] * dat.Grid_Grid_mid.unit_conversion
    y_data = getattr(var, "all_time_data")[t0,:] * unit_conv
    
    y_data = one_dim_grid(grid[:,t0,:], grid_conv, x_data, y_data)
    
    l1.set_xdata(x_data)
    l2.set_xdata(x_data)
    l3.set_xdata(x_data)
    
    xlimH = np.max(x_data)
    xlimL = 0.0
    ax.set_xlim([xlimL,xlimH])
    ax.set_ylim([ylimL,ylimH])
    ax_data.new_axis=np.array([xlimL,xlimH,ylimL,ylimH])
    
    ax1.set_ylabel(name + " (" + units + ")")
    #y_data = getattr(var, "all_time_data")[t0,:] * unit_conv
    #l2.set_ydata(y_data)
    #l2.set_label("Not Fixed")
    
    l3.set_ydata(y_data)
    l3.set_label("Not Fixed")
    
    ylimL = 0.0 # np.min(y_data)
    if np.max(y_data) != 0.0:
      ylimH = 1.5 * np.max(y_data)
    else:
      ylimH = 1.0
    ax1.set_ylim([ylimL,ylimH])
    xlimL = 0.0
    xlimH = np.max(x_data)
    ax1.set_xlim([xlimL,xlimH])
    ax1_data.new_axis=np.array([xlimL,xlimH,ylimL,ylimH])
    
    ax.set_title(dat.Times.name
        + ' = {0:5.3f}'.format(dat.Times.all_time_data[t0]
        * dat.Times.unit_conversion))
    
    if ii!=0:
      ax_data.zoomed_axis2 = ax_data.zoomed_axis1[0:2]
          #* ((ax_data.new_axis[1] - ax_data.new_axis[0]) 
          #/ (ax_data.old_axis[1] - ax_data.old_axis[0]))
      ax_data.zoomed_axis2 = np.append(ax_data.zoomed_axis2, ((ax_data.new_axis[3]
          - ax_data.new_axis[2]) / (ax_data.old_axis[3]- ax_data.old_axis[2]))
          * ax_data.zoomed_axis1[2:])
      
      ax.set_xlim(ax_data.zoomed_axis2[:2])
      ax.set_ylim(ax_data.zoomed_axis2[2:])
      
    if iv!=0 and ii!=0:
      ax1_data.zoomed_axis2 = ((ax1_data.new_axis[1] - ax1_data.new_axis[0])
          / (ax1_data.old_axis[1] - ax1_data.old_axis[0]))* ax1_data.zoomed_axis1[0:2]
      ax1_data.zoomed_axis2 = np.append(ax1_data.zoomed_axis2, ((ax1_data.new_axis[3]
          - ax1_data.new_axis[2]) / (ax1_data.old_axis[3]- ax1_data.old_axis[2]))
          * ax1_data.zoomed_axis1[2:])
      
      ax1.set_ylim(ax1_data.zoomed_axis2[2:])
    
    ax1.legend(loc='upper right')
    ax_data.old_axis = ax_data.new_axis
    ax1_data.old_axis = ax1_data.new_axis
    
    ii+=1
    iv+=1
    setattr(fig,'reset_axis',ii)
    setattr(fig,'reset_axis2',iv)
    
    fig.canvas.draw_idle()
  
  fig.canvas.mpl_connect('key_press_event', press)
  stime.on_changed(update)
  radio.on_clicked(change_variable)
  update(istart)

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
