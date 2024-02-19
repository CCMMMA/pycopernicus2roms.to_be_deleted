import netCDF4
import numpy as np
from scipy.interpolate import griddata


def main():
    dataset_ini = netCDF4.Dataset('data/ini.nc', mode='w', format='NETCDF4_CLASSIC')
    dataset_bry = netCDF4.Dataset('data/bry.nc', mode='w', format='NETCDF4_CLASSIC')

    dataset_grid = netCDF4.Dataset('data/Campania_max200m_withC3andC4_angle0_hmin5.nc')
    var_lat_rho_dst = dataset_grid.variables['lat_rho'][:]
    var_lon_rho_dst = dataset_grid.variables['lon_rho'][:]

    dataset_ssh = netCDF4.Dataset('data/myoc_d00_20240215_ssh.nc')
    var_lat = dataset_ssh.variables['latitude'][:]
    var_lon = dataset_ssh.variables['longitude'][:]
    var_zos = dataset_ssh.variables['zos'][0]
    var_lon_src, var_lat_src = np.meshgrid(var_lon, var_lat)

    var_zos_dst = griddata(
        (var_lat_src.flatten(), var_lon_src.flatten()),
        var_zos.filled(fill_value=np.nan).flatten(),
        (var_lat_rho_dst, var_lon_rho_dst),
        fill_value=1.e+20,
        method="linear"
    )

    dataset_bry.close()
    dataset_ini.close()

    import matplotlib.pyplot as plt

    var_zos_dst[var_zos_dst == 1e20] = np.nan
    plt.figure(figsize=(10, 6))
    plt.imshow(var_zos_dst, cmap='jet', origin='lower')
    plt.colorbar()
    plt.show()


if __name__ == '__main__':
    main()
