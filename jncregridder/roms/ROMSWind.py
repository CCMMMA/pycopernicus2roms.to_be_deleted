import netCDF4 as nc
import os
import multiprocessing
import numpy as np
from concurrent.futures import ProcessPoolExecutor

from jncregridder.util.Interpolator import BilinearInterpolator
from jncregridder.wrf.WRFData import WRFData


def calcStress(u10mMatrix, v10mMatrix, ust, angleMatrix, slpMatrix, t2mMatrix):
    etaU, xiU = u10mMatrix.shape

    rotU10m = u10mMatrix * np.cos(angleMatrix) + v10mMatrix * np.sin(angleMatrix)
    rotV10m = v10mMatrix * np.cos(angleMatrix) - u10mMatrix * np.sin(angleMatrix)

    aRotUV10m = np.arctan2(-rotV10m, -rotU10m + 1e-8)
    rotUV10m = np.sqrt(rotU10m ** 2 + rotV10m ** 2)
    rhoAir = slpMatrix * 100 / (287.058 * (t2mMatrix + 273.15))

    condition_ust_999 = (ust == 999.9)
    condition_ust_888 = (ust == 888.8)

    dcU = np.where(condition_ust_999 | condition_ust_888, ust, np.nan)

    if condition_ust_999:
        dcU_ust_999 = np.where(condition_ust_999 & (rotUV10m > 7.5), (0.8 + 0.065 * rotUV10m) / 1000, 0)
        dcU = np.where(condition_ust_999, dcU_ust_999, dcU)

    elif condition_ust_888:
        dcU_lt_3 = np.where(rotUV10m < 3, 2.17 / 1000, 0)
        dcU_between_3_and_6 = np.where((3 <= rotUV10m) & (rotUV10m <= 6),
                                       (0.29 + 3.1 / rotUV10m + 7.7 / (rotUV10m * rotUV10m)) / 1000, 0)
        dcU_between_6_and_26 = np.where((6 < rotUV10m) & (rotUV10m <= 26), (0.60 + 0.070 * rotUV10m) / 1000, 0)
        dcU_else = np.where(
            ~((rotUV10m < 3) | ((3 <= rotUV10m) & (rotUV10m <= 6)) | ((6 < rotUV10m) & (rotUV10m <= 26))), 2.42 / 1000,
            0)

        dcU = dcU_lt_3 + dcU_between_3_and_6 + dcU_between_6_and_26 + dcU_else

    result = np.zeros((etaU, xiU, 2), dtype=np.float32)
    result[..., 0] = -(rhoAir * (dcU ** 2) * np.sin(aRotUV10m))*1000
    result[..., 1] = -(rhoAir * (dcU ** 2) * np.cos(aRotUV10m))*1000

    return result[..., :]


def process_wrf_file(wrfData, LATRHO, LONRHO, MASKRHO):
    print(f"Computing {wrfData['path']} ...")

    bilinearInterpolatorRho = BilinearInterpolator(wrfData["XLAT"], wrfData["XLONG"], LATRHO, LONRHO, MASKRHO)

    # print("Interpolating T2")
    T2 = bilinearInterpolatorRho.simpleInterp(wrfData["T2"][0], 1e37)

    # print("Interpolating SLP")
    SLP = bilinearInterpolatorRho.simpleInterp(wrfData["SLP"], 1e37)

    # print("Interpolating U10M")
    U10M = bilinearInterpolatorRho.simpleInterp(wrfData["U10M"], 1e37)

    # print("Interpolating V10M")
    V10M = bilinearInterpolatorRho.simpleInterp(wrfData["V10M"], 1e37)

    # print("Interpolating RH2")
    RH2 = bilinearInterpolatorRho.simpleInterp(wrfData["RH2"], 1e37)

    # print("Interpolating SWDOWN")
    SWDOWN = bilinearInterpolatorRho.simpleInterp(wrfData["SWDOWN"], 1e37)

    # print("Interpolating GLW")
    GLW = bilinearInterpolatorRho.simpleInterp(wrfData["GLW"], 1e37)

    result = {
        "T2": T2,
        "SLP": SLP,
        "U10M": U10M,
        "V10M": V10M,
        "RH2": RH2,
        "SWDOWN": SWDOWN,
        "GLW": GLW,
        "date": wrfData['dModDate']
    }

    return result


class ROMSWind:
    def __init__(self, url, romsGrid, wrfFilenameOrFilesPath):
        self.url = url
        self.romsGrid = romsGrid
        self.totalTime = 1
        self.wrfFiles = []

        folder = os.path.abspath(wrfFilenameOrFilesPath)
        if os.path.isdir(folder):
            listOfFiles = sorted(os.listdir(folder))
            for filename in listOfFiles:
                filepath = os.path.join(folder, filename)
                if os.path.isfile(filepath) and filename.startswith("wrf"):
                    wrfData = WRFData(filepath)
                    wrfFile = wrfData.getData()
                    self.wrfFiles.append(wrfFile)

            self.totalTime = len(self.wrfFiles)
        else:
            wrfData = WRFData(wrfFilenameOrFilesPath)
            wrfFile = wrfData.getData()
            self.wrfFiles.append(wrfFile)

        self.ncfWritable = nc.Dataset(url, 'w', format='NETCDF4_CLASSIC')
        self.ncfWritable.createDimension(self.romsGrid.dimEtaRho.name, len(self.romsGrid.dimEtaRho))
        self.ncfWritable.createDimension(self.romsGrid.dimXiRho.name, len(self.romsGrid.dimXiRho))
        self.ncfWritable.createDimension(self.romsGrid.dimEtaU.name, len(self.romsGrid.dimEtaU))
        self.ncfWritable.createDimension(self.romsGrid.dimXiU.name, len(self.romsGrid.dimXiU))
        self.ncfWritable.createDimension(self.romsGrid.dimEtaV.name, len(self.romsGrid.dimEtaV))
        self.ncfWritable.createDimension(self.romsGrid.dimXiV.name, len(self.romsGrid.dimXiV))

        self.ncfWritable.createDimension('ocean_time', self.totalTime)

        self.lat = self.ncfWritable.createVariable('lat', 'f8', ('eta_rho', 'xi_rho'))
        self.lat.long_name = "latitude of RHO-points"
        self.lat.units = "degree_north"
        self.lat.field = "lat_rho, scalar"
        self.lat.standard_name = "latitude"
        self.lat._CoordinateAxisType = "Lat"

        self.lon = self.ncfWritable.createVariable('lon', 'f8', ('eta_rho', 'xi_rho'))
        self.lon.long_name = "longitude of RHO-points"
        self.lon.units = "degree_east"
        self.lon.field = "lon_rho, scalar"
        self.lon.standard_name = "longitude"
        self.lon._CoordinateAxisType = "Lon"

        self.time = self.ncfWritable.createVariable('time', 'f8', ('ocean_time',))
        self.time.long_name = "atmospheric forcing time"
        self.time.units = "days since 1968-05-23 00:00:00 GMT"
        self.time.calendar = "gregorian"

        self.Pair = self.ncfWritable.createVariable('Pair', 'f4', ('ocean_time', 'eta_rho', 'xi_rho'))
        self.Pair.long_name = "Mean Sea Level Pressure"
        self.Pair.units = "millibar"
        self.Pair.time = "time"
        self.Pair.missing_value = 1e37

        self.Tair = self.ncfWritable.createVariable('Tair', 'f4', ('ocean_time', 'eta_rho', 'xi_rho'))
        self.Tair.long_name = "Air Temperature (2m)"
        self.Tair.units = "Celsius"
        self.Tair.time = "time"
        self.Tair.missing_value = 1e37

        self.Qair = self.ncfWritable.createVariable('Qair', 'f4', ('ocean_time', 'eta_rho', 'xi_rho'))
        self.Qair.long_name = "Relative Humidity (2m)"
        self.Qair.units = "percentage"
        self.Qair.time = "time"
        self.Qair.missing_value = 1e37

        self.rain = self.ncfWritable.createVariable('rain', 'f4', ('ocean_time', 'eta_rho', 'xi_rho'))
        self.rain.long_name = "Rain fall rate"
        self.rain.units = "kilogram meter-2 second-1"
        self.rain.time = "time"
        self.rain.missing_value = 1e37

        self.swrad = self.ncfWritable.createVariable('swrad', 'f4', ('ocean_time', 'eta_rho', 'xi_rho'))
        self.swrad.long_name = "Solar shortwave radiation"
        self.swrad.units = "watt meter-2"
        self.swrad.time = "time"
        self.swrad.positive_value = "downward flux, heating"
        self.swrad.negative_value = "upward flux, cooling"
        self.swrad.missing_value = 1e37

        self.lwrad_down = self.ncfWritable.createVariable('lwrad_down', 'f4', ('ocean_time', 'eta_rho', 'xi_rho'))
        self.lwrad_down.long_name = "Net longwave radiation flux"
        self.lwrad_down.units = "watt meter-2"
        self.lwrad_down.time = "time"
        self.lwrad_down.positive_value = "downward flux, heating"
        self.lwrad_down.negative_value = "upward flux, cooling"
        self.lwrad_down.missing_value = 1e37

        self.Uwind = self.ncfWritable.createVariable('Uwind', 'f4', ('ocean_time', 'eta_rho', 'xi_rho'))
        self.Uwind.long_name = "Wind velocity, u-component (m s-1)"
        self.Uwind.units = "m s-1"
        self.Uwind.time = "time"
        self.Uwind.missing_value = 1e37

        self.Vwind = self.ncfWritable.createVariable('Vwind', 'f4', ('ocean_time', 'eta_rho', 'xi_rho'))
        self.Vwind.long_name = "Wind velocity, v-component (m s-1)"
        self.Vwind.units = "m s-1"
        self.Vwind.time = "time"
        self.Vwind.missing_value = 1e37

        self.lat_u = self.ncfWritable.createVariable('lat_u', 'f8', ('eta_u', 'xi_u'))
        self.lat_u.long_name = "latitude of RHO-points"
        self.lat_u.units = "degree_north"
        self.lat_u.standard_name = "latitude"
        self.lat_u._CoordinateAxisType = "Lat"

        self.lon_u = self.ncfWritable.createVariable('lon_u', 'f8', ('eta_u', 'xi_u'))
        self.lon_u.long_name = "longitude of U-points"
        self.lon_u.units = "degree_east"
        self.lon_u.standard_name = "longitude"
        self.lon_u._CoordinateAxisType = "Lon"

        self.lat_v = self.ncfWritable.createVariable('lat_v', 'f8', ('eta_v', 'xi_v'))
        self.lat_v.long_name = "latitude of V-points"
        self.lat_v.units = "degree_north"
        self.lat_v.standard_name = "latitude"
        self.lat_v._CoordinateAxisType = "Lat"

        self.lon_v = self.ncfWritable.createVariable('lon_v', 'f8', ('eta_v', 'xi_v'))
        self.lon_v.long_name = "longitude of V-points"
        self.lon_v.units = "degree_east"
        self.lon_v.standard_name = "longitude"
        self.lon_v._CoordinateAxisType = "Lon"

        self.ocean_time = self.ncfWritable.createVariable('ocean_time', 'f8', ('ocean_time',))
        self.ocean_time.long_name = "surface ocean time"
        self.ocean_time.units = "days since 1968-05-23 00:00:00 GMT"
        self.ocean_time.calendar = "gregorian"

        self.sustr = self.ncfWritable.createVariable('sustr', 'f4', ('ocean_time', 'eta_u', 'xi_u'))
        self.sustr.long_name = "Kinematic wind stress, u-component (m2 s-2)"
        self.sustr.units = "Newton meter-2"
        self.sustr.scale_factor = 1000.0
        self.sustr.time = "ocean_time"
        self.sustr.missing_value = 1e37

        self.svstr = self.ncfWritable.createVariable('svstr', 'f4', ('ocean_time', 'eta_v', 'xi_v'))
        self.svstr.long_name = "Kinematic wind stress, v-component (m2 s-2)"
        self.svstr.units = "Newton meter-2"
        self.svstr.scale_factor = 1000.0
        self.svstr.time = "ocean_time"
        self.svstr.missing_value = 1e37

        self.lon[:] = self.romsGrid.LONRHO
        self.lat[:] = self.romsGrid.LATRHO
        self.lon_u[:] = self.romsGrid.LONU
        self.lat_u[:] = self.romsGrid.LATU
        self.lon_v[:] = self.romsGrid.LONV
        self.lat_v[:] = self.romsGrid.LATV

    def make(self):
        num_processors = multiprocessing.cpu_count()
        with ProcessPoolExecutor(max_workers=num_processors) as executor:
            futures = [executor.submit(process_wrf_file, wrfFile, self.romsGrid.LATRHO, self.romsGrid.LONRHO, self.romsGrid.MASKRHO)
                       for wrfFile in self.wrfFiles]

            for k, future in enumerate(futures):
                res = future.result()

                self.time[k] = res["date"]
                self.ocean_time[k] = res["date"]

                self.Uwind[k] = res["U10M"]
                self.Vwind[k] = res["V10M"]
                self.Tair[k] = res["T2"]
                self.Pair[k] = res["SLP"]
                self.Qair[k] = res["RH2"]
                self.swrad[k] = res["SWDOWN"]
                self.lwrad_down[k] = res["GLW"]

                self.romsGrid.ANGLE[:] = 0
                stress = calcStress(res["U10M"], res["V10M"], 888.8, self.romsGrid.ANGLE, res["SLP"], res["T2"])
                interpolator2DU = BilinearInterpolator(self.lat[:], self.lon[:], self.lat_u[:], self.lon_u[:], self.romsGrid.MASKU, method="nearest")
                # print("Interpolating stressU")
                stressU = interpolator2DU.simpleInterp(stress[..., 0], 1e37)
                interpolator2DV = BilinearInterpolator(self.lat[:], self.lon[:], self.lat_v[:], self.lon_v[:], self.romsGrid.MASKV, method="nearest")
                # print("Interpolating stressV")
                stressV = interpolator2DV.simpleInterp(stress[..., 1], 1e37)
                self.sustr[k] = stressU
                self.svstr[k] = stressV

    def close(self):
        if self.ncfWritable:
            self.ncfWritable.close()
