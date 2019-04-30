# Duncan Barlow, Odin project, Warwick University, 01/19
import sdf_helper as sh
import matplotlib.pyplot as plt
import numpy as np
import glob
import sys, os
from matplotlib.widgets import Slider, RadioButtons
import analysis_functions as afunc


def preallocate_dat(dat, iend):
  setattr(dat, "Times", np.zeros(RunCounter))
  
  for var_name in dat.variables:
    var = getattr(dat, var_name)
    data = getattr(var, "data")
    len_x = np.shape(data)[0]
    setattr(var, "all_time_data" , np.zeros((RunCounter, len_x)))
  
  setattr(dat.Radius_mid, "all_time", np.zeros((RunCounter, len_x)))
  
  for var_name in dat.variables_time:
    var = getattr(dat, var_name)
    setattr(var, "all_time_data" , np.zeros(RunCounter))


def add_label(dat):
  
  # default labels
  for var_name in dat.variables:
    var = getattr(dat, var_name)
    name = getattr(var, "name")
    units = getattr(var, "units")
    
    label = name + " (" + units + ")"
    setattr(var, "label", label)
    setattr(var, "unit_conversion", 1)
    
  # User defined
  var = getattr(dat, "Fluid_Temperature_electron")
  label = 'Electron Temperature (keV)'
  setattr(var, "label", label)
  unit_conversion = 1.0 / 11604.5 / 1000.0 # from Kelvin
  setattr(var, "unit_conversion", unit_conversion)
  
  return dat



class read_sdf:
  """ This class will collect all the data from an sdf file. The units and
  their conversions from SI are set here
  """
  def __init__(self, *args, **kwargs):
    RunCounter = kwargs.get('RunCounter', 0)
    nmat =  kwargs.get('nmat', 0)
    len_x =  kwargs.get('len_x', 0)
    
    self.com = np.zeros(RunCounter)
    self.time = np.zeros(RunCounter)
    self.time_units = 'Time (ns)'
    self.time_conversion = 1e9
    self.max_rho = np.zeros(RunCounter)
    self.tot_laser_dep = np.zeros(RunCounter)
    
    self.nmat = nmat
    self.material_names = [None] * nmat
    self.material_Volume_fraction = np.zeros((RunCounter, nmat, len_x))
                    
    # 1D
    self.radius = np.zeros((RunCounter,len_x))
    
    # 2D
    self.X = []
    self.X_units = 'x ($\mu$m)'
    self.X_conversion = 1e6
    self.Y = []
    self.Y_units = 'y ($\mu$m)'
    self.Y_conversion = 1e6
    
    self.Fluid_Rho = np.zeros((RunCounter,len_x))
    self.Fluid_Rho_units = 'Density (g/cm$^3$)'
    self.Fluid_Rho_conversion = 1.0 / 1000.0
    
    self.Fluid_Temperature_ion = np.zeros((RunCounter,len_x))
    self.Fluid_Temperature_electron = np.zeros((RunCounter,len_x))
    self.Fluid_Temperature_units = 'Temperature (keV)'
    self.Fluid_Temperature_conversion = 1.0 / 11604.5 / 1000.0 # from Kelvin
    
    self.Fluid_Pressure_ion = np.zeros((RunCounter,len_x))
    self.Fluid_Pressure_electron = np.zeros((RunCounter,len_x))
    self.Fluid_Pressure_units = 'Pressure (Mbar)'
    self.Fluid_Pressure_conversion = 1.0e-11 # from Pascal
    
    self.Fluid_Energy_ion = np.zeros((RunCounter,len_x))
    self.Fluid_Energy_electron = np.zeros((RunCounter,len_x))
    self.Fluid_Energy_units = 'Energy (J/kg)'
    self.Fluid_Energy_conversion = 1.0
    
    # These values are overwritten (ie do not change here)
    self.User_defined = np.zeros((RunCounter,len_x))
    self.User_defined_units = 'Place holder'
    self.User_defined_conversion =  1.0




def use_sdf(n, pathname, use_analysis):
  SDFName=pathname+'/'+str(n).zfill(4)+'.sdf'
  dat = sh.getdata(SDFName,verbose=False)
  
  dat_names = list(dat.__dict__.keys())
  variable_type = type(dat.Fluid_Rho)
  dat_variables_names = []
  
  for n in range(0, len(dat_names)):
    var = getattr(dat, dat_names[n])
    if type(var) == variable_type:
      dat_variables_names.append(dat_names[n])
  setattr(dat, "variables", dat_variables_names)
  
  if use_analysis:
    dat = afunc.basic(dat)
    dat = afunc.laser(dat)
    # energy
    # adiabat
  
  dat = add_label(dat)
  
  return dat



def get_data_all(dat1, istart, iend, pathname, use_analysis, cs):
  """
  """
  dat1 = use_sdf(dat1)
  dat1 = preallocate_dat(dat1, iend)
  for n in range(istart, iend):     
    dat = use_sdf(n, pathname, use_analysis)
    # grid for this data is either radius or X depending on rz t/f
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
  return all_time



def main():
        """
        """


if __name__ == "__main__":
        main()
