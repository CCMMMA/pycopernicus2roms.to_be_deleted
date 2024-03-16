from netCDF4 import Dataset
import numpy as np


class CopernicusBase:
    def __init__(self, url):
        self.ncDataset = Dataset(url)

        self.dimTime = self.ncDataset.dimensions.get("time")
        self.dimLat = self.ncDataset.dimensions.get("latitude")
        self.dimLon = self.ncDataset.dimensions.get("longitude")
        self.dimDepth = self.ncDataset.dimensions.get("depth")
        self.b3D = self.dimDepth is not None

        self.time = self.dimTime.size
        self.lat = self.dimLat.size
        self.lon = self.dimLon.size

        self.LAT = self.__load__("latitude")
        self.LON = self.__load__("longitude")
        self.TIME = self.__load__("time")

        self.DEPTH = None
        self.Z = None
        if self.b3D:
            self.depth = self.dimDepth.size
            self.DEPTH = self.__load__("depth")
            self.Z = self.__getZ__()

    def __getZ__(self):
        Z = np.ones((self.depth, self.lat, self.lon)) * self.DEPTH[:, np.newaxis, np.newaxis]
        return Z

    def __load__(self, variable_name):
        try:
            return self.ncDataset.variables[variable_name][:]
        except KeyError:
            raise Exception(f"Variable {variable_name} not found in the dataset.")
        except Exception as e:
            raise Exception(f"Error loading variable {variable_name}: {str(e)}")


class CopernicusCur(CopernicusBase):
    def __init__(self, url):
        super().__init__(url)
        self.VO = self.__load__("vo")
        self.UO = self.__load__("uo")
        self.FillValue = 1e20


class CopernicusSSH(CopernicusBase):
    def __init__(self, url):
        super().__init__(url)
        self.ZOS = self.__load__("zos")
        self.FillValue = 1e20


class CopernicusSal(CopernicusBase):
    def __init__(self, url):
        super().__init__(url)
        self.SO = self.__load__("so")
        self.FillValue = 1e20


class CopernicusTem(CopernicusBase):
    def __init__(self, url):
        super().__init__(url)
        self.THETAO = self.__load__("thetao")
        self.BOTTOMT = self.__load__("bottomT")
        self.FillValue = 1e20
