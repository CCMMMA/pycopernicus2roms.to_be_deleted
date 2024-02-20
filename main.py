import netCDF4
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import time


def interpolate(lon_dst, lat_dst, mask, lon_src, lat_src, var):
    lon_src_mesh, lat_src_mesh = np.meshgrid(lon_src, lat_src)

    var_dst = griddata(
        (lat_src_mesh.flatten(), lon_src_mesh.flatten()),
        var.filled(fill_value=np.nan).flatten(),
        (lat_dst, lon_dst),
        method="linear"
    )

    indices = np.where(mask == 0)
    var_dst[indices] = 1e37

    rows, cols = np.where(~np.isnan(var_dst) & (var_dst != 1e37))
    values = var_dst[rows, cols]
    interp_rows, interp_cols = np.where(np.isnan(var_dst))
    interpolated_values = griddata((rows, cols), values, (interp_rows, interp_cols), method='nearest')
    var_dst[np.isnan(var_dst)] = interpolated_values

    indices = np.where(var_dst == 1e37)
    var_dst[indices] = np.nan

    return var_dst


def process_variable(var, lon_dst, lat_dst, mask, var_lon_src, var_lat_src):
    var_interpolated = interpolate(lon_dst, lat_dst, mask, var_lon_src, var_lat_src, var)

    # fig = plt.figure()
    # ax = fig.add_axes([0, 0, 1, 1])
    # ax.pcolormesh(lon_dst, lat_dst, var_interpolated, cmap='jet')
    # plt.show()


def main():
    dataset_grid = netCDF4.Dataset('data/Campania_max200m_withC3andC4_angle0_hmin5.nc')
    lat_dst = dataset_grid.variables['lat_rho'][:]
    lon_dst = dataset_grid.variables['lon_rho'][:]
    mask = dataset_grid.variables['mask_rho'][:]

    start_time = time.time()

    directory = "data"
    for filename in tqdm(os.listdir(directory)):
        if filename.startswith('myoc_d00'):
            file_path = os.path.join(directory, filename)
            dataset = netCDF4.Dataset(file_path)

            for var_name in dataset.variables:
                var = dataset.variables[var_name]
                dimensions = var.dimensions

                var_lat_src = dataset.variables['latitude'][:]
                var_lon_src = dataset.variables['longitude'][:]

                print("var: ", var_name)

                if len(dimensions) == 3:
                    k, _, _ = dataset.variables[var_name].shape

                    for i in tqdm(range(k)):
                        process_variable(var[i], lon_dst, lat_dst, mask, var_lon_src, var_lat_src)

                elif len(dimensions) == 4:
                    t, k, _, _ = dataset.variables[var_name].shape

                    for i in tqdm(range(t)):
                        for j in tqdm(range(k)):
                            process_variable(var[i][j], lon_dst, lat_dst, mask, var_lon_src, var_lat_src)

    # dataset_ini = netCDF4.Dataset('data/ini.nc', mode='w', format='NETCDF4_CLASSIC')
    # dataset_bry = netCDF4.Dataset('data/bry.nc', mode='w', format='NETCDF4_CLASSIC')
    # dataset_bry.close()
    # dataset_ini.close()

    end_time = time.time()
    execution_time = end_time - start_time
    print("Execution time:", execution_time, "seconds")


if __name__ == '__main__':
    main()
