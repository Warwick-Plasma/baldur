# Duncan Barlow, Odin project, Warwick University, 01/19
import sdf_helper as sh
import matplotlib.pyplot as plt
import numpy as np
import glob
import sys, os
from matplotlib.widgets import Slider, RadioButtons
import analysis_functions as afunc



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
    print(dat_variables_names)
    setattr(dat, 'variables', dat_variables_names)
    return dat



def get_data_one(one_sdf, n, pathname, var_name, use_analysis):
        """
        """

        SDFName=pathname+'/'+str(n).zfill(4)+'.sdf'
        dat = sh.getdata(SDFName,verbose=False)

        len_x = np.shape(dat.Fluid_Rho.data)[0]
        len_y = np.shape(dat.Fluid_Rho.data)[1]

        rho = dat.Fluid_Rho.data
        one_sdf.X = x
        one_sdf.Y = y
        one_sdf.radius = np.sqrt(xc**2 + yc**2)
        
        if use_analysis:
                one_sdf = afunc.basic(dat, one_sdf)
                #one_sdf = afunc.energy(dat, one_sdf)
                one_sdf = afunc.laser(dat, one_sdf)
                #afunc.adiabat
        else:
                one_sdf.com = 0
                one_sdf.tot_laser_dep = 0

        one_sdf.time = t
        one_sdf.max_rho = np.max(np.max(rho))

        one_sdf.nmat = dat.Integer_flags.nmat
        one_sdf.material_names = [None] * one_sdf.nmat
        one_sdf.material_Volume_fraction = np.zeros((one_sdf.nmat, len_x, len_y))
        if one_sdf.nmat != 1:
                for nm in range(1, one_sdf.nmat+1):
                        one_sdf.material_names[nm-1] = getattr(getattr(dat, 'material_string_flags_'
                            + str(nm).zfill(3)),'name_')
                        one_sdf.material_Volume_fraction[nm-1,:,:] = getattr(getattr(dat, 'Fluid_Volume_fraction_'
                            + one_sdf.material_names[nm-1]),'data')[:,:]

        one_sdf.Fluid_Rho = rho[:,:]
        one_sdf.Fluid_Temperature_ion = dat.Fluid_Temperature_ion.data
        one_sdf.Fluid_Temperature_electron = dat.Fluid_Temperature_electron.data
        one_sdf.Fluid_Pressure_ion = dat.Fluid_Pressure_ion.data
        one_sdf.Fluid_Pressure_electron = dat.Fluid_Pressure_electron.data
        one_sdf.Fluid_Energy_ion = dat.Fluid_Energy_ion.data
        one_sdf.Fluid_Energy_electron = dat.Fluid_Energy_electron.data
        
        # If nessacery change units and conversion here
        var = getattr(dat, var_name)
        one_sdf.User_defined = var.data
        one_sdf.User_defined_units = var.name + ' $(' + var.units + ')$'
        one_sdf.User_defined_conversion = 1.0
        return one_sdf




def get_data_all(all_time, minrun, pathname, var_name, use_analysis, RunCounter, cs, nmat):
        """
        """
        for n in range(minrun,RunCounter):
                
                one_sdf = read_sdf()
                one_sdf = get_data_one(one_sdf, n, pathname, var_name, use_analysis)
                
                all_time.nmat = nmat
                
                all_time.com[n] = one_sdf.com
                all_time.time[n] = one_sdf.time
                all_time.max_rho[n] = one_sdf.max_rho
                all_time.tot_laser_dep[n] = one_sdf.tot_laser_dep
                all_time.material_names = one_sdf.material_names
                all_time.material_Volume_fraction[n,:,:] = one_sdf.material_Volume_fraction[:,:,cs]
                all_time.radius[n,:] = one_sdf.radius[:,cs]
                all_time.Fluid_Rho[n,:] = one_sdf.Fluid_Rho[:,cs]
                all_time.Fluid_Temperature_ion[n,:] = one_sdf.Fluid_Temperature_ion[:,cs]
                all_time.Fluid_Temperature_electron[n,:] = one_sdf.Fluid_Temperature_electron[:,cs]
                all_time.Fluid_Pressure_ion[n,:] = one_sdf.Fluid_Pressure_ion[:,cs]
                all_time.Fluid_Pressure_electron[n,:] = one_sdf.Fluid_Pressure_electron[:,cs]
                all_time.Fluid_Energy_ion[n,:] = one_sdf.Fluid_Energy_ion[:,cs]
                all_time.Fluid_Energy_electron[n,:] = one_sdf.Fluid_Energy_electron[:,cs]
                
                all_time.User_defined = one_sdf.User_defined[:,cs]
                all_time.User_defined_units = one_sdf.User_defined_units
                all_time.User_defined_conversion = one_sdf.User_defined_conversion
        return all_time



def main():
        """
        """


if __name__ == "__main__":
        main()
