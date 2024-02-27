import netCDF4 as nc
import numpy as np
from datetime import datetime


class ROMSBoundary:
    def __init__(self, url, romsGrid, forcingTimeSteps):
        self.url = url
        self.romsGrid = romsGrid
        self.forcingTimeSteps = forcingTimeSteps

        self.OCEAN_TIME = None
        self.UBAR = None
        self.VBAR = None
        self.ZETA = None
        self.TEMP = None
        self.SALT = None
        self.U = None
        self.V = None

        self.ncfWritable = nc.Dataset(url, 'w', format='NETCDF4_CLASSIC')
        self.ncfWritable.setncattr('type', 'Boundary forcing file')
        self.ncfWritable.setncattr('title', 'Boundary forcing file (BRY) used for forcing of the ROMS model')
        self.ncfWritable.setncattr('grd_file', romsGrid.url)
        self.ncfWritable.setncattr('source', 'University of Napoli Parthenope Weather Centre http://meteo.uniparthenope.it')

        self.ncfWritable.setncattr('date', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

        self.dimEtaRho = self.ncfWritable.createDimension(romsGrid.dimEtaRho.name, len(romsGrid.dimEtaRho))
        self.dimXiRho = self.ncfWritable.createDimension(romsGrid.dimXiRho.name, len(romsGrid.dimXiRho))
        self.dimEtaU = self.ncfWritable.createDimension(romsGrid.dimEtaU.name, len(romsGrid.dimEtaU))
        self.dimXiU = self.ncfWritable.createDimension(romsGrid.dimXiU.name, len(romsGrid.dimXiU))
        self.dimEtaV = self.ncfWritable.createDimension(romsGrid.dimEtaV.name, len(romsGrid.dimEtaV))
        self.dimXiV = self.ncfWritable.createDimension(romsGrid.dimXiV.name, len(romsGrid.dimXiV))
        self.dimSRho = self.ncfWritable.createDimension('s_rho', len(romsGrid.s_rho))
        self.dimSw = self.ncfWritable.createDimension('s_w', len(romsGrid.s_w))

        self.dimOceanTime = self.ncfWritable.createDimension('ocean_time', forcingTimeSteps)

        self.dimOne = self.ncfWritable.createDimension('one', 1)
        self.dimTwo = self.ncfWritable.createDimension('two', 2)
        self.dimFour = self.ncfWritable.createDimension('four', 4)
        self.dimBath = self.ncfWritable.createDimension('bath', 1)

        self.ocean_time = self.ncfWritable.createVariable('ocean_time', np.float64, ('ocean_time',))
        self.ocean_time.long_name = 'ocean forcing time'
        self.ocean_time.units = 'days since 1968-05-23 00:00:00 GMT'
        self.ocean_time.calendar = 'gregorian'

        self.angle = self.ncfWritable.createVariable("angle", np.double, ("eta_rho", "xi_rho"))
        self.angle.long_name = "angle between xu axis and east"
        self.angle.units = "radiant"

        self.theta_b = self.ncfWritable.createVariable("theta_b", np.double, ("one",))
        self.theta_b.long_name = "S-coordinate surface control parameter"
        self.theta_b.units = "nondimensional"

        self.theta_s = self.ncfWritable.createVariable("theta_s", np.double, ("one",))
        self.theta_s.long_name = "S-coordinate bottom control parameter"
        self.theta_s.units = "nondimensional"

        self.Tcline = self.ncfWritable.createVariable("Tcline", np.double, ())

        self.z_r = self.ncfWritable.createVariable("z_r", np.double, ("s_rho", "eta_rho", "xi_rho"))
        self.z_r.long_name = "Sigma layer to depth matrix"
        self.z_r.units = "meter"

        self.hc = self.ncfWritable.createVariable("hc", np.double, ("one",))
        self.hc.long_name = "S-coordinate parameter, critical depth"
        self.hc.units = "meter"

        self.Cs_w = self.ncfWritable.createVariable("Cs_w", np.double, ("s_w",))
        self.Cs_w.long_name = "S-coordinate stretching curves at W-points"
        self.Cs_w.valid_min = -1
        self.Cs_w.valid_max = 0
        self.Cs_w.field = "s_w, scalar"

        self.Cs_r = self.ncfWritable.createVariable("Cs_r", np.double, ("s_rho",))
        self.Cs_r.long_name = "S-coordinate stretching curves at RHO-points"
        self.Cs_r.units = "nondimensional"

        self.s_w = self.ncfWritable.createVariable("s_w", np.double, ("s_w",))
        self.s_w.long_name = "S-coordinate at W-points"
        self.s_w.valid_min = -1
        self.s_w.valid_max = 0
        self.s_w.standard_name = "ocean_s_coordinate_g1"
        self.s_w.formula_terms = "s: s_w C: Cs_w eta: zeta depth: h depth_c: hc"
        self.s_w.field = "s_w, scalar"

        self.sc_r = self.ncfWritable.createVariable("sc_r", np.double, ("s_rho",))
        self.sc_r.long_name = "S-coordinate at RHO-points"
        self.sc_r.units = "nondimensional"

        self.s_rho = self.ncfWritable.createVariable("s_rho", np.double, ("s_rho",))
        self.s_rho.long_name = "oS-coordinate at RHO-points"
        self.s_rho.valid_min = -1.0
        self.s_rho.valid_max = 0.0
        self.s_rho.positive = "up"
        self.s_rho.standard_name = "ocean_s_coordinate_g1"
        self.s_rho.formula_terms = ""
        self.s_rho.field = "s_rho, scalar"
        self.s_rho._CoordinateTransformType = "Vertical"
        self.s_rho._CoordinateAxes = "s_rho"
        self.s_rho._CoordinateAxisType = "GeoZ"
        self.s_rho._CoordinateZisPositive = "up"

        self.h = self.ncfWritable.createVariable("h", np.double, ("eta_rho", "xi_rho"))
        self.h.long_name = "Final bathymetry at RHO-points"
        self.h.units = "meters"
        self.h.field = "bath, scalar"

        self.lat_rho = self.ncfWritable.createVariable("lat_rho", np.double, ("eta_rho", "xi_rho"))
        self.lat_rho.long_name = "latitude of RHO-points"
        self.lat_rho.units = "degree_north"
        self.lat_rho.field = "lat_rho, scalar"
        self.lat_rho.standard_name = "latitude"
        self.lat_rho._CoordinateAxisType = "Lat"

        self.lon_rho = self.ncfWritable.createVariable("lon_rho", np.double, ("eta_rho", "xi_rho"))
        self.lon_rho.long_name = "longitude of RHO-points"
        self.lon_rho.units = "degree_east"
        self.lon_rho.field = "lon_rho, scalar"
        self.lon_rho.standard_name = "longitude"
        self.lon_rho._CoordinateAxisType = "Lon"

        self.lat_u = self.ncfWritable.createVariable("lat_u", np.double, ("eta_u", "xi_u"))
        self.lat_u.long_name = "latitude of RHO-points"
        self.lat_u.units = "degree_north"
        self.lat_u.standard_name = "latitude"
        self. lat_u._CoordinateAxisType = "Lat"

        self.lon_u = self.ncfWritable.createVariable("lon_u", np.double, ("eta_u", "xi_u"))
        self.lon_u.long_name = "longitude of U-points"
        self.lon_u.units = "degree_east"
        self.lon_u.standard_name = "longitude"
        self.lon_u._CoordinateAxisType = "Lon"

        self.lat_v = self.ncfWritable.createVariable("lat_v", np.double, ("eta_v", "xi_v"))
        self.lat_v.long_name = "latitude of V-points"
        self.lat_v.units = "degree_north"
        self.lat_v.standard_name = "latitude"
        self.lat_v._CoordinateAxisType = "Lat"

        self.lon_v = self.ncfWritable.createVariable("lon_v", np.double, ("eta_v", "xi_v"))
        self.lon_v.long_name = "longitude of V-points"
        self.lon_v.units = "degree_east"
        self.lon_v.standard_name = "longitude"
        self.lon_v._CoordinateAxisType = "Lon"

        self.temp_west = self.ncfWritable.createVariable("temp_west", np.double, ("ocean_time", "s_rho", "eta_rho"))
        self.temp_west.long_name = "potential temperature western boundary conditions"
        self.temp_west.units = "Celsius"
        self.temp_west.field = "temp_west, scalar, series"
        self.temp_west.missing_value = romsGrid.no_data
        self.temp_west.time = "ocean_time"

        self.temp_east = self.ncfWritable.createVariable("temp_east", np.double, ("ocean_time", "s_rho", "eta_rho"))
        self.temp_east.long_name = "potential temperature eastern boundary conditions"
        self.temp_east.units = "Celsius"
        self.temp_east.field = "temp_east, scalar, series"
        self.temp_east.missing_value = romsGrid.no_data
        self.temp_east.time = "ocean_time"

        self.temp_south = self.ncfWritable.createVariable("temp_south", np.double, ("ocean_time", "s_rho", "xi_rho"))
        self.temp_south.long_name = "potential temperature southern boundary conditions"
        self.temp_south.units = "Celsius"
        self.temp_south.field = "temp_south, scalar, series"
        self.temp_south.missing_value = romsGrid.no_data
        self.temp_south.time = "ocean_time"

        self.temp_north = self.ncfWritable.createVariable("temp_north", np.double, ("ocean_time", "s_rho", "xi_rho"))
        self.temp_north.long_name = "potential temperature northern boundary conditions"
        self.temp_north.units = "Celsius"
        self.temp_north.field = "temp_north, scalar, series"
        self.temp_north.missing_value = romsGrid.no_data
        self.temp_north.time = "ocean_time"

        self.salt_west = self.ncfWritable.createVariable("salt_west", np.double, ("ocean_time", "s_rho", "eta_rho"))
        self.salt_west.long_name = "salinity western boundary conditions"
        self.salt_west.units = "PSU"
        self.salt_west.field = "salt_west, scalar, series"
        self.salt_west.missing_value = romsGrid.no_data
        self.salt_west.time = "ocean_time"

        self.salt_east = self.ncfWritable.createVariable("salt_east", np.double, ("ocean_time", "s_rho", "eta_rho"))
        self.salt_east.long_name = "salinity eastern boundary conditions"
        self.salt_east.units = "PSU"
        self.salt_east.field = "salt_east, scalar, series"
        self.salt_east.missing_value = romsGrid.no_data
        self.salt_east.time = "ocean_time"

        self.salt_south = self.ncfWritable.createVariable("salt_south", np.double, ("ocean_time", "s_rho", "xi_rho"))
        self.salt_south.long_name = "salinity southern boundary conditions"
        self.salt_south.units = "PSU"
        self.salt_south.field = "salt_south, scalar, series"
        self.salt_south.missing_value = romsGrid.no_data
        self.salt_south.time = "ocean_time"

        self.salt_north = self.ncfWritable.createVariable("salt_north", np.double, ("ocean_time", "s_rho", "xi_rho"))
        self.salt_north.long_name = "salinity northern boundary conditions"
        self.salt_north.units = "PSU"
        self.salt_north.field = "salt_north, scalar, series"
        self.salt_north.missing_value = romsGrid.no_data
        self.salt_north.time = "ocean_time"

        self.zeta_west = self.ncfWritable.createVariable("zeta_west", np.double, ("ocean_time", "eta_rho"))
        self.zeta_west.long_name = "free-surface western boundary conditions"
        self.zeta_west.units = "meter"
        self.zeta_west.field = "zeta_west, scalar, series"
        self.zeta_west.missing_value = romsGrid.no_data
        self.zeta_west.time = "ocean_time"

        self.zeta_east = self.ncfWritable.createVariable("zeta_east", np.double, ("ocean_time", "eta_rho"))
        self.zeta_east.long_name = "free-surface eastern boundary conditions"
        self.zeta_east.units = "meter"
        self.zeta_east.field = "zeta_east, scalar, series"
        self.zeta_east.missing_value = romsGrid.no_data
        self.zeta_east.time = "ocean_time"

        self.zeta_south = self.ncfWritable.createVariable("zeta_south", np.double, ("ocean_time", "xi_rho"))
        self.zeta_south.long_name = "free-surface southern boundary conditions"
        self.zeta_south.units = "meter"
        self.zeta_south.field = "zeta_south, scalar, series"
        self.zeta_south.missing_value = romsGrid.no_data
        self.zeta_south.time = "ocean_time"

        self.zeta_north = self.ncfWritable.createVariable("zeta_north", np.double, ("ocean_time", "xi_rho"))
        self.zeta_north.long_name = "free-surface northern boundary conditions"
        self.zeta_north.units = "meter"
        self.zeta_north.field = "zeta_north, scalar, series"
        self.zeta_north.missing_value = romsGrid.no_data
        self.zeta_north.time = "ocean_time"

        self.u_west = self.ncfWritable.createVariable("u_west", np.double, ("ocean_time", "s_rho", "eta_u"))
        self.u_west.long_name = "3D U-momentum western boundary conditions"
        self.u_west.units = "meter second-1"
        self.u_west.field = "u_west, scalar, series"
        self.u_west.missing_value = romsGrid.no_data
        self.u_west.time = "ocean_time"

        self.u_east = self.ncfWritable.createVariable("u_east", np.double, ("ocean_time", "s_rho", "eta_u"))
        self.u_east.long_name = "3D U-momentum eastern boundary conditions"
        self.u_east.units = "meter second-1"
        self.u_east.field = "u_east, scalar, series"
        self.u_east.missing_value = romsGrid.no_data
        self.u_east.time = "ocean_time"

        self.u_south = self.ncfWritable.createVariable("u_south", np.double, ("ocean_time", "s_rho", "xi_u"))
        self.u_south.long_name = "3D U-momentum southern boundary conditions"
        self.u_south.units = "meter second-1"
        self.u_south.field = "u_south, scalar, series"
        self.u_south.missing_value = romsGrid.no_data
        self.u_south.time = "ocean_time"

        self.u_north = self.ncfWritable.createVariable("u_north", np.double, ("ocean_time", "s_rho", "xi_u"))
        self.u_north.long_name = "3D U-momentum northern boundary conditions"
        self.u_north.units = "meter second-1"
        self.u_north.field = "u_north, scalar, series"
        self.u_north.missing_value = romsGrid.no_data
        self.u_north.time = "ocean_time"

        self.v_west = self.ncfWritable.createVariable("v_west", np.double, ("ocean_time", "s_rho", "eta_v"))
        self.v_west.long_name = "3D V-momentum western boundary conditions"
        self.v_west.units = "meter second-1"
        self.v_west.field = "v_west, scalar, series"
        self.v_west.missing_value = romsGrid.no_data
        self.v_west.time = "ocean_time"

        self.v_east = self.ncfWritable.createVariable("v_east", np.double, ("ocean_time", "s_rho", "eta_v"))
        self.v_east.long_name = "3D V-momentum eastern boundary conditions"
        self.v_east.units = "meter second-1"
        self.v_east.field = "v_east, scalar, series"
        self.v_east.missing_value = romsGrid.no_data
        self.v_east.time = "ocean_time"

        self.v_south = self.ncfWritable.createVariable("v_south", np.double, ("ocean_time", "s_rho", "xi_v"))
        self.v_south.long_name = "3D V-momentum southern boundary conditions"
        self.v_south.units = "meter second-1"
        self.v_south.field = "v_south, scalar, series"
        self.v_south.missing_value = romsGrid.no_data
        self.v_south.time = "ocean_time"

        self.v_north = self.ncfWritable.createVariable("v_north", np.double, ("ocean_time", "s_rho", "xi_v"))
        self.v_north.long_name = "3D V-momentum northern boundary conditions"
        self.v_north.units = "meter second-1"
        self.v_north.field = "v_north, scalar, series"
        self.v_north.missing_value = romsGrid.no_data
        self.v_north.time = "ocean_time"

        self.vbar_west = self.ncfWritable.createVariable("vbar_west", np.double, ("ocean_time", "eta_v"))
        self.vbar_west.long_name = "2D V-momentum western boundary conditions"
        self.vbar_west.units = "meter second-1"
        self.vbar_west.field = "vbar_west, scalar, series"
        self.vbar_west.missing_value = romsGrid.no_data
        self.vbar_west.time = "ocean_time"

        self.vbar_east = self.ncfWritable.createVariable("vbar_east", np.double, ("ocean_time", "eta_v"))
        self.vbar_east.long_name = "2D V-momentum eastern boundary conditions"
        self.vbar_east.units = "meter second-1"
        self.vbar_east.field = "vbar_east, scalar, series"
        self.vbar_east.missing_value = romsGrid.no_data
        self.vbar_east.time = "ocean_time"

        self.vbar_south = self.ncfWritable.createVariable("vbar_south", np.double, ("ocean_time", "xi_v"))
        self.vbar_south.long_name = "2D V-momentum southern boundary conditions"
        self.vbar_south.units = "meter second-1"
        self.vbar_south.field = "vbar_south, scalar, series"
        self.vbar_south.missing_value = romsGrid.no_data
        self.vbar_south.time = "ocean_time"

        self.vbar_north = self.ncfWritable.createVariable("vbar_north", np.double, ("ocean_time", "xi_v"))
        self.vbar_north.long_name = "2D V-momentum northern boundary conditions"
        self.vbar_north.units = "meter second-1"
        self.vbar_north.field = "vbar_north, scalar, series"
        self.vbar_north.missing_value = romsGrid.no_data
        self.vbar_north.time = "ocean_time"

        self.ubar_west = self.ncfWritable.createVariable("ubar_west", np.double, ("ocean_time", "eta_u"))
        self.ubar_west.long_name = "2D U-momentum western boundary conditions"
        self.ubar_west.units = "meter second-1"
        self.ubar_west.field = "ubar_west, scalar, series"
        self.ubar_west.missing_value = romsGrid.no_data
        self.ubar_west.time = "ocean_time"

        self.ubar_east = self.ncfWritable.createVariable("ubar_east", np.double, ("ocean_time", "eta_u"))
        self.ubar_east.long_name = "2D U-momentum eastern boundary conditions"
        self.ubar_east.units = "meter second-1"
        self.ubar_east.field = "ubar_east, scalar, series"
        self.ubar_east.missing_value = romsGrid.no_data
        self.ubar_east.time = "ocean_time"

        self.ubar_south = self.ncfWritable.createVariable("ubar_south", np.double, ("ocean_time", "xi_u"))
        self.ubar_south.long_name = "2D U-momentum southern boundary conditions"
        self.ubar_south.units = "meter second-1"
        self.ubar_south.field = "ubar_south, scalar, series"
        self.ubar_south.missing_value = romsGrid.no_data
        self.ubar_south.time = "ocean_time"

        self.ubar_north = self.ncfWritable.createVariable("ubar_north", np.double, ("ocean_time", "xi_u"))
        self.ubar_north.long_name = "2D U-momentum northern boundary conditions"
        self.ubar_north.units = "meter second-1"
        self.ubar_north.field = "ubar_north, scalar, series"
        self.ubar_north.missing_value = romsGrid.no_data
        self.ubar_north.time = "ocean_time"

    def make(self):
        self.h[:] = self.romsGrid.H
        self.lon_rho[:] = self.romsGrid.LONRHO
        self.lat_rho[:] = self.romsGrid.LATRHO
        self.lon_u[:] = self.romsGrid.LONU
        self.lat_u[:] = self.romsGrid.LATU
        self.lon_v[:] = self.romsGrid.LONV
        self.lat_v[:] = self.romsGrid.LATV

        self.theta_s[:] = self.romsGrid.theta_s
        self.theta_b[:] = self.romsGrid.theta_b
        self.Tcline[:] = self.romsGrid.TCLINE

        self.ocean_time[:] = self.OCEAN_TIME

        self.sc_r[:] = self.romsGrid.s_rho
        self.s_rho[:] = self.romsGrid.s_rho
        self.hc[:] = self.romsGrid.HC
        self.Cs_r[:] = self.romsGrid.cs_r

    def write(self, time):
        self.u_west[time] = self.U[:, :, 0]
        self.u_east[time] = self.U[:, :, -1]
        self.u_south[time] = self.U[:, 0, :]
        self.u_north[time] = self.U[:, -1, :]

        self.v_west[time] = self.V[:, :, 0]
        self.v_east[time] = self.V[:, :, -1]
        self.v_south[time] = self.V[:, 0, :]
        self.v_north[time] = self.V[:, -1, :]

        self.salt_west[time] = self.SALT[:, :, 0]
        self.salt_east[time] = self.SALT[:, :, -1]
        self.salt_south[time] = self.SALT[:, 0, :]
        self.salt_north[time] = self.SALT[:, -1, :]

        self.temp_west[time] = self.TEMP[:, :, 0]
        self.temp_east[time] = self.TEMP[:, :, -1]
        self.temp_south[time] = self.TEMP[:, 0, :]
        self.temp_north[time] = self.TEMP[:, -1, :]

        self.zeta_west[time] = self.ZETA[:, 0]
        self.zeta_east[time] = self.ZETA[:, -1]
        self.zeta_south[time] = self.ZETA[0, :]
        self.zeta_north[time] = self.ZETA[-1, :]

        self.ubar_west[time] = self.UBAR[:, 0]
        self.ubar_east[time] = self.UBAR[:, -1]
        self.ubar_south[time] = self.UBAR[0, :]
        self.ubar_north[time] = self.UBAR[-1, :]

        self.vbar_west[time] = self.VBAR[:, 0]
        self.vbar_east[time] = self.VBAR[:, -1]
        self.vbar_south[time] = self.VBAR[0, :]
        self.vbar_north[time] = self.VBAR[-1, :]

    def close(self):
        if self.ncfWritable:
            self.ncfWritable.close()