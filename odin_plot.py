# Duncan Barlow, Odin project, Warwick University, 01/19
import sdf_helper as sh
import matplotlib.pyplot as plt
import numpy as np
import glob
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
	plt.pcolormesh(X,Y,np.log10(np.transpose(all_isentrope)), vmin=-1, vmax=0.5)
	plt.xlabel('Time (s)')
	plt.ylabel('Cumalitive Mass (kg)')
	cbar = plt.colorbar()
	cbar.set_label('log(Isentrope/Adiabat)')
	
	plt.figure()
	plt.pcolormesh(X,Y,np.transpose(all_pressure_ls), vmin=10000, vmax=100000)
	plt.xlabel('Time (s)')
	plt.ylabel('Cumalitive Mass (kg)')
	cbar = plt.colorbar()
	cbar.set_label('|d(ln(p))/dr|')
	
	plt.show()



def snapshot(start, *args, **kwargs):
	"""
	"""
	var_name = kwargs.get('var_name', "Fluid_Rho")

	pathname = os.path.abspath(os.getcwd())
	SDFName=pathname+'/'+str(start).zfill(4)+'.sdf'
	dat = sh.getdata(SDFName,verbose=False)
	runs = glob.glob1(pathname,"*.sdf")
	RunCounter = len(runs)
	minrun = start
	
	var=getattr(dat, var_name)
	X=var.grid.data[0][:,:]
	Y=var.grid.data[1][:,:]

	plt.ion()
	plt.close('all')

	fig=plt.figure(num=1,figsize=(6,6),facecolor='white')
	ax1=plt.axes([0.1, 0.1, 0.8, 0.8])
	cmesh = ax1.pcolormesh(X,Y,var.data[:,:], linewidth=0.1)
	#plt.gca().set_aspect('equal', adjustable='box')
	ax1.set_xlim([np.min(X[:-1,:]),np.max(X[:-1,:])])
	ax1.set_ylim([np.min(Y[:-1,:]),np.max(Y[:-1,:])])
	plt.xlabel(var.grid.labels[0] + ' $(' + var.grid.units[1] + ')$')
	plt.ylabel(var.grid.labels[1] + ' $(' + var.grid.units[1] + ')$')
	plt.title('Time '+str(t*10**9)+'ns')
	cbar = plt.colorbar(cmesh)
	cbar.set_label(var.name + ' $(' + var.units + ')$')

	axcolor='lightgoldenrodyellow' # slider background colour
	axtime=plt.axes([0.02, 0.02, 0.90, 0.03], facecolor=axcolor) # slider size
	stime=Slider(axtime, ' ', minrun, RunCounter-1, valinit=minrun, valfmt='%1.0f')
	
	axcolor = 'lightgoldenrodyellow'
	rax = plt.axes([0.15, 0.735, 0.15, 0.15], facecolor=axcolor)
	radio = RadioButtons(rax, (var_name, "Fluid_Pressure_ion",
		 "Fluid_Pressure_electron", "Fluid_Energy_ion",
		 "Fluid_Energy_electron", "Fluid_Rho", "Fluid_Temperature_ion",
		 "Fluid_Temperature_electron"))

	grid_colour = "none"
	setattr(fig, "grid_colour", grid_colour)

	def change_variable(label):
		max_x=stime.val
		stime.set_val(max_x)

	def press(event): # Defining a function for keyboard inputs for o and p
		max_x=stime.val
		if event.key == 'o' and max_x > minrun+0.5:
			stime.set_val(max_x-1) # define the increment
		elif event.key == 'p' and max_x < RunCounter-1.5:
			stime.set_val(max_x+1)
		elif event.key == 'i':
			grid_colour = getattr(fig, "grid_colour")
			if grid_colour == "none":
				grid_colour = "black"
			elif grid_colour == "black":
				grid_colour = "none"
			setattr(fig, "grid_colour", grid_colour)
			stime.set_val(max_x)

	def update(val):
		"""
		"""
		x=ax1.get_xlim()
		y=ax1.get_ylim()
		zoomed_axis1=np.array([x[0],x[1],y[0],y[1]])

		sdf_num=int(round(stime.val))
		SDFName=pathname+'/'+str(sdf_num).zfill(4)+'.sdf'
		dat = sh.getdata(SDFName,verbose=False)

		var_name = radio.value_selected
		var = getattr(dat, var_name)
		X=var.grid.data[0][:,:]
		Y=var.grid.data[1][:,:]
		
		ax1.clear()
		cmesh = ax1.pcolormesh(X,Y,var.data[:,:], linewidth=0.01)
		grid_colour = getattr(fig, "grid_colour")
		cmesh.set_edgecolor(grid_colour)
		ax1.set_title('Time {0:5.3f}'.format(t*10**9)+'ns')
		cbar.set_clim(np.min(var.data[:,:]), np.max(var.data[:,:]))
		cbar.set_label(var.name + ' $(' + var.units + ')$')
		cbar.draw_all()
		#cmesh.set_array(var.data[:,:].ravel())
		
		ax1.set_xlim(zoomed_axis1[:2])
		ax1.set_ylim(zoomed_axis1[2:])

		fig.canvas.draw_idle()

	fig.canvas.mpl_connect('key_press_event', press) # this makes the keyboard function work (key press)
	stime.on_changed(update) # this makes the slider work
	radio.on_clicked(change_variable)
	plt.show()



def total_energy(output_number):
	""" Internal, kinetic and pv energy
	"""
	dat=sh.getdata(output_number, verbose=False)

	fac = 1.0
	if dat.Logical_flags.use_rz:
		fac = 2*np.pi

	vol = dat.Fluid_Volume.data *	fac
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
	corner_mass = dat.Test_Corner_Mass.data *	fac
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



def lineout(start, *args, **kwargs):
	"""
	"""
	var_name = kwargs.get('var_name', "Fluid_Temperature")
	
	pathname = os.path.abspath(os.getcwd())
	SDFName=pathname+'/'+str(start).zfill(4)+'.sdf'
	dat = sh.getdata(SDFName,verbose=False)
	nmat = dat.Integer_flags.nmat # assume nmat doesn't change
	runs = glob.glob1(pathname,"*.sdf")
	RunCounter = len(runs)-1
	RunCounter = kwargs.get('end', RunCounter) + 1
	minrun = start
	
	len_x = np.shape(dat.Fluid_Rho.data)[0]
	len_y = np.shape(dat.Fluid_Rho.data)[1]
	half = round(len_y / 2)
	cs = kwargs.get('cross_section', half)
	
	all_time = all_sdf(RunCounter, nmat, len_x)
	all_time = get_data_all(minrun, RunCounter, nmat, pathname, cs, all_time)

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
	
	ax.set_xlabel('Radius (m)', fontsize = fs)
	ax.tick_params(axis='x', labelcolor = 'black', labelsize = fs)
	
	ax.set_ylabel('Density (g/cm$^3$)', fontsize = fs)	
	l1, = ax.plot(1, lw = 2.5, color='black')
	ax.tick_params(axis='y', labelcolor='black', labelsize = fs)
	
	ax1.set_ylabel('Temperature (keV)', color='tab:red', fontsize = fs)
	l2, = ax1.plot(0, lw = 2, color = 'tab:blue')
	l3, = ax1.plot(0, lw = 2, color = 'tab:red')
	ax1.tick_params(axis='y', labelcolor = 'tab:red', labelsize = fs)

	xmin = 0 #np.min(all_time.radius[0,:])
	xmax = np.max(all_time.radius[0,:])
	ax.set_xlim([xmin, xmax])
	
	ax3 = ax2.twinx()
	
	ax2.set_xlabel('Time (ns)', fontsize = fs)
	ax2.tick_params(axis='x', labelcolor = 'black', labelsize = fs)
	x_data = all_time.time * 1.0e9
	
	# Plot laser power deposited
	ax2.set_ylabel('Laser power deposited (W)', fontsize = fs)
	dt = all_time.time[1:] - all_time.time[:-1]
	y_data = (all_time.tot_laser_dep[1:] - all_time.tot_laser_dep[:-1]) / dt
	y_data = np.insert(y_data, 0, 0)
	la = ax2.plot(x_data, y_data, lw = 2, color='black')
	ax2.tick_params(axis='y', labelcolor='black', labelsize = fs)
	
	# Plot maximum density
	ax3.set_ylabel('Maximum density', color='tab:red', fontsize = fs)
	y_data = all_time.max_rho / 1000
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
	stime=Slider(axtime, 'SDF selector', minrun, RunCounter-1, valinit = minrun, valfmt = '%1.0f')
	stime.label.set_size(fs)
	
	axcolor = 'lightgoldenrodyellow'
	rax = plt.axes([0.02, 0.65, 0.18, 0.18], facecolor=axcolor)
	radio = RadioButtons(rax, ("Fluid_Temperature", "Fluid_Pressure",
	    "Fluid_Energy", var_name))
	
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
		if event.key == 'o' and t0 > minrun + 0.5:
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
				ax.set_ylabel('Log(Density (g/cm$^3$))')	
			else:
				setattr(fig,'log_on',0)
				ax.set_ylabel('Density (g/cm$^3$)')	
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
		
		l1.set_xdata(all_time.radius[t0,:])
		l2.set_xdata(all_time.radius[t0,:])
		l3.set_xdata(all_time.radius[t0,:])
		
		log_on = getattr(fig,'log_on')
		if log_on == 0:
			y_data = all_time.Fluid_Rho[t0,:] / 1000 # convert to g/cm^-3
			l1.set_ydata(y_data)
		
			ylimL = 0.0
			ylimH = 1.5 * np.max(y_data)
		else:
			y_data = all_time.Fluid_Rho[t0,:] / 1000 # convert to g/cm^-3
			l1.set_ydata(np.log10(y_data))
		
			ylimL = np.log10(0.0001)#0.0
			ylimH = np.log10(10000)#1.5 * np.max(y_data)
		
		ax.set_ylim([ylimL,ylimH])
		xlimH = np.max(all_time.radius[t0,:])
		xlimL = 0.0
		ax.set_xlim([xlimL,xlimH])
		ax_data.new_axis=np.array([xlimL,xlimH,ylimL,ylimH])
		
		var_name = radio.value_selected
		var_conv = var_name + '_conversion'
		unit_conv = getattr(all_time, var_conv)
		var_units = var_name + '_units'
		units = getattr(all_time, var_units)
		try:
			var_ion = var_name + '_ion'
			var_ele = var_name + '_electron'
			ti = getattr(all_time, var_ion)
			te = getattr(all_time, var_ele)
		except:
			var = getattr(all_time, var_name)
		
		ax1.set_ylabel(units)
		y_data = ti[t0,:] * unit_conv
		l2.set_ydata(y_data)
		l2.set_label('Ion')
		
		y_data = te[t0,:] * unit_conv
		l3.set_ydata(y_data)
		l3.set_label('Electron')
		
		ylimL = 0.0 # np.min(y_data)
		ylimH = 1.5 * np.max(y_data)
		ax1.set_ylim([ylimL,ylimH])
		xlimH = np.max(all_time.radius[t0,:])
		xlimL = 0.0
		ax1.set_xlim([xlimL,xlimH])
		ax1_data.new_axis=np.array([xlimL,xlimH,ylimL,ylimH])
		
		ax.set_title('Time {0:5.3f}'.format(all_time.time[t0]*1e9) + ' $ns$')
		
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
	update(minrun)

	plt.show()



class axis_data:
	""" This class will collect all the data from the sdf files
	"""
	def __init__(self):
		self.new_axis = np.zeros(4)
		self.old_axis = np.zeros(4)
		self.zoomed_axis1 = np.zeros(4)
		self.zoomed_axis2 = np.zeros(4)



class all_sdf:
	""" This class will collect all the data from the sdf files
	"""
	def __init__(self, RunCounter, nmat, len_x):
		self.com = np.zeros(RunCounter)
		self.time = np.zeros(RunCounter)
		self.max_rho = np.zeros(RunCounter)
		self.tot_laser_dep = np.zeros(RunCounter)
		
		self.nmat = nmat
		self.material_names = [None] * nmat
		self.material_Volume_fraction = np.zeros((RunCounter, len_x, nmat))
		
		self.radius = np.zeros((RunCounter,len_x))
		self.Fluid_Rho = np.zeros((RunCounter,len_x))
		
		self.Fluid_Temperature_ion = np.zeros((RunCounter,len_x))
		self.Fluid_Temperature_electron = np.zeros((RunCounter,len_x))
		self.Fluid_Temperature_units = 'Temperature (keV)'
		self.Fluid_Temperature_conversion = 1.0 / 11604.5 / 1000.0
		
		self.Fluid_Pressure_ion = np.zeros((RunCounter,len_x))
		self.Fluid_Pressure_electron = np.zeros((RunCounter,len_x))
		self.Fluid_Pressure_units = 'Pressure (Mbar)'
		self.Fluid_Pressure_conversion = 1.0e-11
		
		self.Fluid_Energy_ion = np.zeros((RunCounter,len_x))
		self.Fluid_Energy_electron = np.zeros((RunCounter,len_x))
		self.Fluid_Energy_units = 'Energy (J/kg)'
		self.Fluid_Energy_conversion = 1.0
		
		self.var = np.zeros((RunCounter,len_x))



def get_data_all(minrun, RunCounter, nmat, pathname, cs, all_time):
	"""
	"""
	for n in range(minrun,RunCounter):
		SDFName=pathname+'/'+str(n).zfill(4)+'.sdf'
		dat = sh.getdata(SDFName,verbose=False)
		fac = 1.0
		if dat.Logical_flags.use_rz:
			fac = 2*np.pi

		vol=dat.Fluid_Volume.data * fac
		rho = dat.Fluid_Rho.data
		mass = rho[:,:] * vol[:,:]
		rad = np.sqrt(xc**2 + yc**2)
		laser_dep = dat.Fluid_Energy_deposited_laser.data
	
		all_time.com[n] = np.sum(np.sum(mass * rad)) / np.sum(np.sum(mass))
		all_time.time[n] = t
		all_time.max_rho[n] = np.max(np.max(rho))
		all_time.tot_laser_dep[n] = np.sum(np.sum(mass * laser_dep))
		
		all_time.nmat = dat.Integer_flags.nmat
		for nm in range(1, nmat+1):
			all_time.material_names[nm-1] = getattr(getattr(dat, 'material_string_flags_'
			    + str(nm).zfill(3)),'name_')
			all_time.material_Volume_fraction[n,:,nm-1] = getattr(getattr(dat, 'Fluid_Volume_fraction_'
			    + all_time.material_names[nm-1]),'data')[:,cs]
		all_time.radius[n,:] = rad[:,cs]
		all_time.Fluid_Rho[n,:] = rho[:,cs]
		all_time.Fluid_Temperature_ion[n,:] = dat.Fluid_Temperature_ion.data[:,cs]
		all_time.Fluid_Temperature_electron[n,:] = dat.Fluid_Temperature_electron.data[:,cs]
		all_time.Fluid_Pressure_ion[n,:] = dat.Fluid_Pressure_ion.data[:,cs]
		all_time.Fluid_Pressure_electron[n,:] = dat.Fluid_Pressure_electron.data[:,cs]
		all_time.Fluid_Energy_ion[n,:] = dat.Fluid_Energy_ion.data[:,cs]
		all_time.Fluid_Energy_electron[n,:] = dat.Fluid_Energy_electron.data[:,cs]
	print(all_time.material_Volume_fraction)
	return all_time



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
