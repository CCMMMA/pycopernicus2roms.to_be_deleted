import datetime
import numpy as np
from netCDF4 import Dataset


class WRFData():
    def __init__(self, url, interpMethod, interpLevels):
        self.url = url
        self.ncDataset = Dataset(url)
        self.interpMethod = interpMethod
        self.interpLevels = interpLevels
        self.tacc = 0
        self.numberOfZLevs = 0
        self.extrapolate = True

        self.VERTICALTYPE_NONE = 0
        self.VERTICALTYPE_P = 1
        self.VERTICALTYPE_Z = 2

        self.simStartDate = None
        self.mapProj = None
        self.trueLat1 = None
        self.trueLat2 = None
        self.standLon = None
        self.title = None
        self.iProgram = None
        self.westEastDim = None
        self.southNorthDim = None
        self.bottomTopDim = None
        self.dX = None
        self.dY = None
        self.cenLat = None
        self.cenLon = None
        self.moadCenLat = None
        self.bucketMM = None
        self.bucketJ = None
        self.poleLat = None
        self.poleLon = None
        self.verticalType = None
        
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

                gcSimStartDate = datetime.datetime(year, month, day, hour, minute, second)
                # print("Simulation start date:", gcSimStartDate.year, "-", gcSimStartDate.month, "-", gcSimStartDate.day, "Z", gcSimStartDate.hour)
            elif attr_name == "MAP_PROJ":
                self.mapProj = int(attr_value)
            elif attr_name == "TRUELAT1":
                self.trueLat1 = float(attr_value)
            elif attr_name == "TRUELAT2":
                self.trueLat2 = float(attr_value)
            elif attr_name == "STAND_LON":
                self.standLon = float(attr_value)
            elif attr_name == "TITLE":
                self.title = attr_value
                if "OUTPUT FROM GEOGRID" in self.title:
                    self.iProgram = 1
                elif "OUTPUT FROM GRIDGEN" in self.title:
                    self.iProgram = 1
                elif "OUTPUT FROM METGRID" in self.title:
                    self.iProgram = 3
                elif "OUTPUT FROM OBSGRID" in self.title:
                    self.iProgram = 3
                elif "OUTPUT FROM REAL_EM" in self.title:
                    self.iProgram = 6
                elif "OUTPUT FROM WRF" in self.title:
                    self.iProgram = 8
                else:
                    raise Exception("Unknown file format: " + title)
            elif attr_name == "WEST-EAST_GRID_DIMENSION":
                self.westEastDim = int(attr_value) - 1
            elif attr_name == "SOUTH-NORTH_GRID_DIMENSION":
                self.southNorthDim = int(attr_value) - 1
            elif attr_name == "BOTTOM-TOP_GRID_DIMENSION":
                self.bottomTopDim = int(attr_value) - 1
                if self.iProgram <= 1:
                    self.bottomTopDim = 24
            elif attr_name == "DX":
                self.dX = float(attr_value)
            elif attr_name == "DY":
                self.dY = float(attr_value)
            elif attr_name == "CEN_LAT":
                self.cenLat = float(attr_value)
            elif attr_name == "CEN_LON":
                self.cenLon = float(attr_value)
            elif attr_name == "MOAD_CEN_LAT":
                self.moadCenLat = float(attr_value)
            elif attr_name == "BUCKET_MM":
                self.bucketMM = float(attr_value)
                if self.bucketMM < 0:
                    self.bucketMM = 0
            elif attr_name == "BUCKET_J":
                self.bucketJ = float(attr_value)
                if self.bucketJ < 0:
                    self.bucketJ = 0
            elif attr_name == "POLE_LAT":
                self.poleLat = float(attr_value)
                if self.poleLat < 0:
                    self.poleLat = 0
            elif attr_name == "POLE_LON":
                self.poleLon = float(attr_value)
                if self.poleLon < 0:
                    self.poleLon = 0

        self.dimTime = self.ncDataset.dimensions.get("Time")
        self.dimDateStrLen = self.ncDataset.dimensions.get("DateStrLen")
        self.dimWestEastStag = self.ncDataset.dimensions.get("west_east_stag")
        self.dimSouthNordStag = self.ncDataset.dimensions.get("south_north_stag")
        self.dimBottomTopStag = self.ncDataset.dimensions.get("bottom_top_stag")

        self.timeDim = len(self.dimTime)
        self.dateStrLenDim = len(self.dimDateStrLen)
        self.southNorthStagDim = len(self.dimSouthNordStag)
        self.westEastStagDim = len(self.dimWestEastStag)
        self.bottomTopStagDim = len(self.dimBottomTopStag)

        if self.interpMethod == -1 or self.interpMethod == 0:
            self.interpLevels = []
            for i in range(100):
                self.interpLevels = -99999.

        if abs(self.interpMethod) == 1:
            self.verticalType = self.VERTICALTYPE_Z

        if self.interpMethod == 1 and self.interpLevels[0] > 100:
            self.verticalType = self.VERTICALTYPE_P

        if self.interpMethod == 1:
            self.numberOfZLevs = len(self.interpLevels)

        if self.extrapolate:
            if self.verticalType != self.VERTICALTYPE_P and self.verticalType != self.VERTICALTYPE_Z:
                self.extrapolate = False
                # print("WARNING: Can only extrapolate when interpolating to pressure/height fields")

        if self.iProgram > 6 and self.interpMethod != 0:
            self.__getInterpInfo()
        else:
            self.extrapolate = False
            self.verticalType = self.VERTICALTYPE_NONE
            self.numberOfZLevs = self.bottomTopDim

        self.XLAT = self.__load("XLONG")
        self.XLONG = self.__load("XLAT")

    def __getInterpInfo(self):
        found = 0
        locOfMinZ = []

        # print("getInterpInfo...")

        if self.verticalType == self.VERTICALTYPE_P:
            # print("getInterpInfo: verticalTyep=p")

            if "P" not in self.ncDataset.variables:
                found += 1
            if "PB" not in self.ncDataset.variables:
                found += 1

            if found != 0 and self.iProgram == 6:
                found = 0
                # print("INFO: probably old wrfinput data - Try getting MU and MUB")

                if "QVAPOR" not in self.ncDataset.variables:
                    found += 1
                if "MU" not in self.ncDataset.variables:
                    found += 1
                if "MUB" not in self.ncDataset.variables:
                    found += 1
                if "P_TOP" not in self.ncDataset.variables:
                    found += 1
                if "ZNU" not in self.ncDataset.variables:
                    found += 1
                if "ZNW" not in self.ncDataset.variables:
                    found += 1
            
            if found != 0:
                # print("WARNING: Asked to interpolate to PRESSURE, but we don't have enough information. Will output data on MODEL LEVELS")
                self.verticalType = self.VERTICALTYPE_NONE
                self.extrapolate = False
            # else:
            #    print("Interpolating to PRESSURE levels")

        if self.verticalType == self.VERTICALTYPE_Z and self.interpMethod == 1:
            # print("getInterpInfo: verticalTyep=z and interpMethod==1")

            if "P" not in self.ncDataset.variables:
                found += 1
            if "PB" not in self.ncDataset.variables:
                found += 1

            if found != 0:
                # print("WARNING: Asked to interpolate to USER SPECIFIED HEIGHT, but we don't have enough information. Will output data on MODEL LEVELS")
                self.verticalTye = self.VERTICALTYPE_NONE
                self.extrapolate = False
            # else:
            #    print("Interpolating to USER SPECIFIED HEIGHT levels")

        if self.verticalType == self.VERTICALTYPE_Z and self.interpMethod == -1:
            # print("getInterpInfo: verticalTyep=z and interpMethod=-1")

            found = 0

            varPH = self.ncDataset.variables.get("PH")
            if varPH is not None:
                PHtmp = np.zeros(varPH.shape, dtype=float)
                aTmp3D = varPH[:]
                for k in range(1, varPH.shape[0]):
                    for i in range(varPH.shape[1]):
                        for j in range(varPH.shape[2]):
                            PHtmp[k][i][j] = 0.5 * (aTmp3D[k - 1, i, j] - 1) + aTmp3D[k, i, j]
            else:
                found += 1

            varPHB = ncDataset.variables.get("PHB")
            if varPHB is not None:
                PHBtmp = np.zeros(varPHB.shape, dtype=float)
                aTmp3D = varPHB[:]
                for k in range(1, varPHB.shape[0]):
                    for i in range(varPHB.shape[1]):
                        for j in range(varPHB.shape[2]):
                            PHBtmp[k][i][j] = 0.5 * (aTmp3D[k - 1, i, j] - 1) + aTmp3D[k, i, j]
            else:
                found += 1

            if found == 0:
                for k in range(varPH.shape[0]):
                    for i in range(varPH.shape[1]):
                        for j in range(varPH.shape[2]):
                            PHtmp[k][i][j] = (PHtmp[k][i][j] + PHBtmp[k][i][j]) / (9.81 * 1000)

                self.numberOfZLevs = self.bottomTopDim
                locOfMinZ = np.unravel_index(np.argmin(PHtmp[0]), PHtmp[0].shape)

                self.interpLevels = np.zeros(self.numberOfZLevs)
                for k in range(self.numberOfZLevs):
                    self.interpLevels[k] = PHtmp[k][locOfMinZ[0]][locOfMinZ[1]]

                self.interpLevels[0] += 0.002
                self.interpLevels[0] = max(self.interpLevels[0], self.interpLevels[1] / 2)  # No negative value
                self.interpLevels[self.numberOfZLevs - 1] -= 0.002

                # print("Interpolating to GENERATED HEIGHT LEVELS.")
            else:
                # print("WARNING: Asked to interpolate to generated height, but we do not have information.\nWill output data on MODEL LEVELS")
                self.verticalTye = self.VERTICALTYPE_NONE
                self.extrapolate = False
            
        # print(f"getInterpInfo: verticalType={self.verticalType} extrapolate={self.extrapolate}")

    def __load(self, variable_name):
        try:
            variable_data = self.ncDataset.variables[variable_name][:]
            return variable_data
        except KeyError:
            raise Exception(f"Variable {variable_name} not found in the dataset.")
        except Exception as e:
            raise Exception(f"Error loading variable {variable_name}: {str(e)}")
