# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

import netCDF4
import numpy as np
from scipy.interpolate import griddata

def do_it():

    dataset_ini = netCDF4.Dataset('data/ini.nc',mode='w',format='NETCDF4_CLASSIC')
    dataset_bry = netCDF4.Dataset('data/bry.nc', mode='w', format='NETCDF4_CLASSIC')

    dataset_grid = netCDF4.Dataset('data/Campania_max200m_withC3andC4_angle0_hmin5.nc')
    var_lat_rho_dst = dataset_grid.variables['lat_rho'][:]
    var_lon_rho_dst = dataset_grid.variables['lon_rho'][:]

    print("lat_rho.shape",var_lat_rho_dst.shape, "lon_rho.shape",var_lon_rho_dst.shape)

    dataset_ssh = netCDF4.Dataset('data/myoc_d00_20240215_ssh.nc')
    var_lat = dataset_ssh.variables['latitude'][:]
    var_lon = dataset_ssh.variables['longitude'][:]
    var_zos = dataset_ssh.variables['zos'][0]

    print(var_zos.shape)
    print(var_lat.shape,var_lon.shape)

    var_lon_src, var_lat_src = np.meshgrid(var_lon, var_lat)

    #print(var_lat_src)

    var_zos_dst = griddata(
        (var_lat_src.flatten(), var_lon_src.flatten()),
        var_zos.flatten(),
        (var_lat_rho_dst, var_lon_rho_dst),
        fill_value=1e20, method="linear")

    print (var_zos_dst)

    dataset_bry.close()
    dataset_ini.close()

    import matplotlib.pyplot as plt

    cs = plt.contour(var_zos_dst)
    plt.show()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    do_it()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
