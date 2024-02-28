from netCDF4 import Dataset
from wrf import getvar
import numpy as np
import jdcal


class WRFData:
    def __init__(self, url, interpLevels):
        self.url = url
        self.interpLevels = interpLevels

        self.simStartDate = None
        self.dModDate = None

        self.ncDataset = Dataset(url)
        self.dimTime = self.ncDataset.dimensions.get("Time")
        self.timeDim = len(self.dimTime)

        attrs = self.ncDataset.__dict__
        for attr_name, attr_value in attrs.items():
            if attr_name == "SIMULATION_START_DATE":
                self.simStartDate = attr_value
                year = int(self.simStartDate[0:4])
                month = int(self.simStartDate[5:7])
                day = int(self.simStartDate[8:10])
                hour = int(self.simStartDate[11:13])
                minute = int(self.simStartDate[14:16])
                second = int(self.simStartDate[17:19])

                dDate = sum(jdcal.gcal2jd(year, month, day)) + hour / 24.0 + minute / 1440.0 + second / 86400.0
                dModOffset = sum(jdcal.gcal2jd(1968, 5, 23))
                self.dModDate = dDate - dModOffset

        self.XLAT = np.array(self.__loadWrf("XLAT"))
        self.XLONG = np.array(self.__loadWrf("XLONG"))
        self.T2 = self.__load("T2")#-273.15
        self.SLP = self.__loadWrf("slp")
        self.UVMET10 = self.__loadWrf("uvmet10")
        self.U10M = self.UVMET10[0]
        self.V10M = self.UVMET10[1]
        self.RH2 = self.__loadWrf("rh2")
        self.CLFR = self.__loadWrf("cloudfrac")
        self.CLF = np.maximum(self.CLFR[0], self.CLFR[1])
        self.SWDOWN = self.__load("SWDOWN")
        self.GLW = self.__load("GLW")

    def __loadWrf(self, variable_name):
        try:
            return getvar(self.ncDataset, variable_name, meta=False)
        except KeyError:
            raise Exception(f"Variable {variable_name} not found in the dataset.")
        except Exception as e:
            raise Exception(f"Error loading variable {variable_name}: {str(e)}")

    def __load(self, variable_name):
        try:
            return self.ncDataset.variables[variable_name][:]
        except KeyError:
            raise Exception(f"Variable {variable_name} not found in the dataset.")
        except Exception as e:
            raise Exception(f"Error loading variable {variable_name}: {str(e)}")
