import netCDF4 as nc
import numpy as np
from datetime import datetime


class ROMSInit:
    def __init__(self, url, romsGrid, forcingTimeSteps):
        self.url = url
        self.romsGrid = romsGrid
        self.forcingTimeSteps = forcingTimeSteps

        self.OCEAN_TIME = None
        self.SCRUM_TIME = None
        self.UBAR = None
        self.VBAR = None
        self.ZETA = None
        self.TEMP = None
        self.SALT = None
        self.U = None
        self.V = None

        self.ncfWritable = nc.Dataset(url, 'w', format='NETCDF4_CLASSIC')
        self.ncfWritable.setncattr('type', 'Initial file')
        self.ncfWritable.setncattr('title', 'Initialization file (INI) used for forcing of the ROMS model')
        self.ncfWritable.setncattr('grd_file', romsGrid.url)
        self.ncfWritable.setncattr('source', 'University of Napoli Parthenope Weather Centre http://meteo.uniparthenope.it')

        self.ncfWritable.setncattr('date', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

        self.dimOceanTime = self.ncfWritable.createDimension('ocean_time', forcingTimeSteps)
        self.dimScrumTime = self.ncfWritable.createDimension('scrum_time', forcingTimeSteps)
        self.dimEtaRho = self.ncfWritable.createDimension(romsGrid.dimEtaRho.name, len(romsGrid.dimEtaRho))
        self.dimXiRho = self.ncfWritable.createDimension(romsGrid.dimXiRho.name, len(romsGrid.dimXiRho))
        self.dimEtaU = self.ncfWritable.createDimension(romsGrid.dimEtaU.name, len(romsGrid.dimEtaU))
        self.dimXiU = self.ncfWritable.createDimension(romsGrid.dimXiU.name, len(romsGrid.dimXiU))
        self.dimEtaV = self.ncfWritable.createDimension(romsGrid.dimEtaV.name, len(romsGrid.dimEtaV))
        self.dimXiV = self.ncfWritable.createDimension(romsGrid.dimXiV.name, len(romsGrid.dimXiV))
        self.dimSRho = self.ncfWritable.createDimension('s_rho', len(romsGrid.s_rho))

        self.dimOne = self.ncfWritable.createDimension('one', 1)

        self.h = self.ncfWritable.createVariable('h', np.float64, ('eta_rho', 'xi_rho'))
        self.h.long_name = 'Final bathymetry at RHO-points'
        self.h.units = 'meters'
        self.h.field = 'bath, scalar'

        self.lat_rho = self.ncfWritable.createVariable('lat_rho', np.float64, ('eta_rho', 'xi_rho'))
        self.lat_rho.long_name = 'latitude of RHO-points'
        self.lat_rho.units = 'degree_north'
        self.lat_rho.field = 'lat_rho, scalar'
        self.lat_rho.standard_name = 'latitude'
        self.lat_rho._CoordinateAxisType = 'Lat'

        self.lon_rho = self.ncfWritable.createVariable('lon_rho', np.float64, ('eta_rho', 'xi_rho'))
        self.lon_rho.long_name = 'longitude of RHO-points'
        self.lon_rho.units = 'degree_east'
        self.lon_rho.field = 'lon_rho, scalar'
        self.lon_rho.standard_name = 'longitude'
        self.lon_rho._CoordinateAxisType = 'Lon'

        self.lat_u = self.ncfWritable.createVariable('lat_u', np.float64, ('eta_u', 'xi_u'))
        self.lat_u.long_name = 'latitude of U-points'
        self.lat_u.units = 'degree_north'
        self.lat_u.standard_name = 'latitude'
        self.lat_u._CoordinateAxisType = 'Lat'

        self.lon_u = self.ncfWritable.createVariable('lon_u', np.float64, ('eta_u', 'xi_u'))
        self.lon_u.long_name = 'longitude of U-points'
        self.lon_u.units = 'degree_east'
        self.lon_u.standard_name = 'longitude'
        self.lon_u._CoordinateAxisType = 'Lon'

        self.lat_v = self.ncfWritable.createVariable('lat_v', np.float64, ('eta_v', 'xi_v'))
        self.lat_v.long_name = 'latitude of V-points'
        self.lat_v.units = 'degree_north'
        self.lat_v.standard_name = 'latitude'
        self.lat_v._CoordinateAxisType = 'Lat'

        self.lon_v = self.ncfWritable.createVariable('lon_v', np.float64, ('eta_v', 'xi_v'))
        self.lon_v.long_name = 'longitude of V-points'
        self.lon_v.units = 'degree_east'
        self.lon_v.standard_name = 'longitude'
        self.lon_v._CoordinateAxisType = 'Lon'

        self.temp = self.ncfWritable.createVariable('temp', np.float64, ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'))
        self.temp.long_name = 'potential temperature'
        self.temp.units = 'Celsius'
        self.temp.coordinates = 'lon_rho lat_rho sc_r ocean_time'
        self.temp.missing_value = 1e37
        self.temp.time = 'ocean_time'

        self.salt = self.ncfWritable.createVariable('salt', np.float64, ('ocean_time', 's_rho', 'eta_rho', 'xi_rho'))
        self.salt.long_name = 'salinity'
        self.salt.units = 'PSU'
        self.salt.coordinates = 'lon_rho lat_rho sc_r ocean_time'
        self.salt.missing_value = 1e37
        self.salt.time = 'ocean_time'

        self.ubar = self.ncfWritable.createVariable('ubar', np.float64, ('ocean_time', 'eta_u', 'xi_u'))
        self.ubar.long_name = 'vertically integrated u-momentum component'
        self.ubar.units = 'meter second-1'
        self.ubar.coordinates = 'lon_u lat_u ocean_time'
        self.ubar.missing_value = 1e37
        self.ubar.time = 'ocean_time'

        self.vbar = self.ncfWritable.createVariable('vbar', np.float64, ('ocean_time', 'eta_v', 'xi_v'))
        self.vbar.long_name = 'vertically integrated v-momentum component'
        self.vbar.units = 'meter second-1'
        self.vbar.coordinates = 'lon_v lat_v ocean_time'
        self.vbar.missing_value = 1e37
        self.vbar.time = 'ocean_time'

        self.u = self.ncfWritable.createVariable('u', np.float64, ('ocean_time', 's_rho', 'eta_u', 'xi_u'))
        self.u.long_name = 'u-momentum component'
        self.u.units = 'meter second-1'
        self.u.coordinates = 'lon_u lat_u sc_r ocean_time'
        self.u.missing_value = 1e37
        self.u.time = 'ocean_time'

        self.v = self.ncfWritable.createVariable('v', np.float64, ('ocean_time', 's_rho', 'eta_v', 'xi_v'))
        self.v.long_name = 'v-momentum component'
        self.v.units = 'meter second-1'
        self.v.coordinates = 'lon_v lat_v sc_r ocean_time'
        self.v.missing_value = 1e37
        self.v.time = 'ocean_time'

        self.zeta = self.ncfWritable.createVariable('zeta', np.float64, ('ocean_time', 'eta_rho', 'xi_rho'))
        self.zeta.long_name = 'free-surface'
        self.zeta.units = 'meter'
        self.zeta.coordinates = 'lon_rho lat_rho ocean_time'
        self.zeta.missing_value = 1e37
        self.zeta.time = 'ocean_time'

        self.theta_b = self.ncfWritable.createVariable('theta_b', np.float64, ('one',))
        self.theta_b.long_name = 'S-coordinate surface control parameter'
        self.theta_b.units = 'nondimensional'

        self.theta_s = self.ncfWritable.createVariable('theta_s', np.float64, ('one',))
        self.theta_s.long_name = 'S-coordinate bottom control parameter'
        self.theta_s.units = 'nondimensional'

        self.Tcline = self.ncfWritable.createVariable('Tcline', np.float64, ('one',))

        self.ocean_time = self.ncfWritable.createVariable('ocean_time', np.float64, ('ocean_time',))
        self.ocean_time.long_name = 'ocean forcing time'
        self.ocean_time.units = 'days since 1968-05-23 00:00:00 GMT'
        self.ocean_time.calendar = 'gregorian'

        self.hc = self.ncfWritable.createVariable('hc', np.float64, ('one',))
        self.hc.long_name = 'S-coordinate parameter, critical depth'
        self.hc.units = 'meter'

        self.scrum_time = self.ncfWritable.createVariable('scrum_time', np.float64, ('scrum_time',))
        self.scrum_time.long_name = 'time since initialization'
        self.scrum_time.units = 'second'

        self.tend = self.ncfWritable.createVariable('tend', np.float64, ('one',))
        self.tend.long_name = 'end processing day'
        self.tend.units = 'day'

        self.Cs_r = self.ncfWritable.createVariable('Cs_r', np.float64, ('s_rho',))
        self.Cs_r.long_name = 'S-coordinate stretching curves at RHO-points'
        self.Cs_r.units = 'nondimensional'

        self.sc_r = self.ncfWritable.createVariable('sc_r', np.float64, ('s_rho',))
        self.sc_r.long_name = 'S-coordinate at RHO-points'
        self.sc_r.units = 'nondimensional'

        self.s_rho = self.ncfWritable.createVariable('s_rho', np.float64, ('s_rho',))
        self.s_rho.long_name = 'oS-coordinate at RHO-points'
        self.s_rho.valid_min = -1.0
        self.s_rho.valid_max = 0.0
        self.s_rho.positive = 'up'
        self.s_rho.standard_name = 'ocean_s_coordinate_g1'
        self.s_rho.formula_terms = ''
        self.s_rho.field = 's_rho, scalar'
        self.s_rho._CoordinateTransformType = 'Vertical'
        self.s_rho._CoordinateAxes = 's_rho'
        self.s_rho._CoordinateAxisType = 'GeoZ'
        self.s_rho._CoordinateZisPositive = 'up'

    def make(self):
        one = len(self.dimOne)

        self.h[:] = self.romsGrid.H
        self.lon_rho[:] = self.romsGrid.LONRHO
        self.lat_rho[:] = self.romsGrid.LATRHO
        self.lon_u[:] = self.romsGrid.LONU
        self.lat_u[:] = self.romsGrid.LATU
        self.lon_v[:] = self.romsGrid.LONV
        self.lat_v[:] = self.romsGrid.LATV

        self.theta_s[:] = self.romsGrid.theta_s
        self.theta_b[:] = self.romsGrid.theta_b
        self.tend[:] = np.zeros(one)
        self.Tcline[:] = self.romsGrid.TCLINE

        self.ocean_time[:] = self.OCEAN_TIME
        self.scrum_time[:] = self.SCRUM_TIME

        self.sc_r[:] = self.romsGrid.s_rho
        self.s_rho[:] = self.romsGrid.s_rho
        self.hc[:] = self.romsGrid.HC
        self.Cs_r[:] = self.romsGrid.cs_r

    def write(self, time):
        self.u[time] = self.U
        self.v[time] = self.V
        self.ubar[time] = self.UBAR
        self.vbar[time] = self.VBAR
        self.zeta[time] = self.ZETA
        self.temp[time] = self.TEMP
        self.salt[time] = self.SALT

    def close(self):
        if self.ncfWritable:
            self.ncfWritable.close()
