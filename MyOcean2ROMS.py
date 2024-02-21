import os
from datetime import datetime, timedelta

from jncregridder.data.copernicus.Copernicus import CopernicusTem, CopernicusSal, CopernicusSSH, CopernicusCur
from jncregridder.roms.ROMSGrid import ROMSGrid
from jncregridder.util.Interpolator import BilinearInterpolator, BilinearInterpolator3D


def getModSimStartDate(ncepDate):
    year = int(ncepDate[:4])
    month = int(ncepDate[4:6])
    day = int(ncepDate[6:8])

    # Create a datetime object for the simulation start date
    sim_start_date = datetime(year, month, day)

    # Define the simulation start date for the model
    model_start_date = datetime(1968, 5, 23)

    # Calculate the timedelta between the simulation start date and the model start date
    offset = sim_start_date - model_start_date

    return model_start_date + offset


class MyOcean2ROMS:
    def __init__(self, gridPath, dataPath, ncepDate, initPath, boundaryPath):
        # Path to the ROMS grid
        self.romsGridPath = gridPath

        # Path to the MyOcean current data
        self.myOceanPathCur = os.path.join(dataPath, f"myoc_d00_{ncepDate}_cur.nc")

        # Path to the MyOcean temperature data
        self.myOceanPathTem = os.path.join(dataPath, f"myoc_d00_{ncepDate}_tem.nc")

        # Path to the MyOcean salinity data
        self.myOceanPathSal = os.path.join(dataPath, f"myoc_d00_{ncepDate}_sal.nc")

        # Path to the MyOcean sea surface height data
        self.myOceanPathSSH = os.path.join(dataPath, f"myoc_d00_{ncepDate}_ssh.nc")

        # Path to the output init file
        self.romsInitPath = initPath

        # Path to the output boundary file
        self.romsBoundaryPath = boundaryPath

        # Open ROMS grid data
        romsGrid = ROMSGrid(gridPath)

        # Get dimension size
        etaRho = romsGrid.etaRho
        xiRho = romsGrid.xiRho
        etaU = romsGrid.etaU
        xiU = romsGrid.xiU
        etaV = romsGrid.etaV
        xiV = romsGrid.xiV

        print("Rho:\n\teta:", etaRho, "\txi:", xiRho)
        print("U:\n\teta:", etaU, "\txi:", xiU)
        print("V:\n\teta:", etaV, "\txi:", xiV)

        # MASK at rho points
        MASKRHO = romsGrid.MASKRHO

        # MASK at u points
        MASKU = romsGrid.MASKU

        # MASK at v points
        MASKV = romsGrid.MASKV

        # LAT,LON at rho points
        LATRHO = romsGrid.LATRHO
        LONRHO = romsGrid.LONRHO

        # Z at rho/sigma pints
        romsZ = romsGrid.Z

        # LAT,LON at u points
        LATU = romsGrid.LATU
        LONU = romsGrid.LONU

        # LAT,LON at v points
        LATV = romsGrid.LATV
        LONV = romsGrid.LONV

        # Get the angle between the xi axis and the real east
        ANGLE = romsGrid.ANGLE

        print("MASKRHO:", MASKRHO.shape)
        print("MASKU:", MASKU.shape)
        print("MASKV:", MASKV.shape)
        print("LATRHO:", LATRHO.shape)
        print("LONRHO:", LONRHO.shape)
        print("LATU:", LATU.shape)
        print("LONU:", LONU.shape)
        print("LATV:", LATV.shape)
        print("LONV:", LONV.shape)
        print("ANGLE:", ANGLE.shape)
        print("ZETA:", romsZ.shape)

        dataTem = CopernicusTem(self.myOceanPathTem)
        dataSal = CopernicusSal(self.myOceanPathSal)
        dataSSH = CopernicusSSH(self.myOceanPathSSH)
        dataCur = CopernicusCur(self.myOceanPathCur)

        # LON at XY points
        LATXY = dataCur.LAT
        LONXY = dataCur.LON
        myOceanZ = dataCur.Z

        forcingTimeSteps = len(dataCur.TIME)

        # TODO: add ROMS INIT file

        oceanTime = []
        for t in range(forcingTimeSteps):
            oceanTime.append((getModSimStartDate(ncepDate) + timedelta(days=t)).strftime("%Y-%m-%d"))

        for t in range(forcingTimeSteps):
            print(f"Time: {t} {oceanTime[t]}")

            # 2D sea surface height
            valuesSSH = dataSSH.ZOS[t]

            # 3D current U component
            valuesU = dataCur.UO[t]

            # 3D current V component
            valuesV = dataCur.VO[t]

            # 3D temperature
            valuesTem = dataTem.ZOS[t]

            # 3D salinity
            valuesSal = dataSal.SO[t]

            # Create a 2D bilinear interpolator on Rho points
            bilinearInterpolatorRho = BilinearInterpolator(LATXY, LONXY, LATRHO, LONRHO, MASKRHO)

            # Interpolate the SSH
            print("Interpolating SSH")
            SSH_ROMS = bilinearInterpolatorRho.interp(valuesSSH, dataSSH.FillValue)

            # Create a 3D bilinear interpolator on Rho points
            interpolator3DRho = BilinearInterpolator3D(LATXY, LONXY, myOceanZ, LATRHO, LONRHO, romsZ, MASKRHO)
            # Create a 3D bilinear interpolator on U points
            interpolator3DU = BilinearInterpolator3D(LATXY, LONXY, myOceanZ, LATRHO, LONRHO, romsZ, MASKU)
            # Create a 3D bilinear interpolator on V points
            interpolator3DV = BilinearInterpolator3D(LATXY, LONXY, myOceanZ, LATRHO, LONRHO, romsZ, MASKV)

            print("Interpolating U")
            U_ROMS = interpolator3DU.interp(valuesU, dataCur.FillValue)

            print("Interpolating V")
            V_ROMS = interpolator3DV.interp(valuesV, dataCur.FillValue)

            print("Interpolating TEMP")
            TEM_ROMS = interpolator3DRho.interp(valuesTem, dataTem.FillValue)

            print("Interpolating SAL")
            SAL_ROMS = interpolator3DRho.interp(valuesSal, dataSal.FillValue)

            # TODO: Add saving to netcdf file


def main():
    gridParh = "data/Campania_new.nc"
    dataPath = "data"
    ncepDate = "20240215"
    initPath = ""
    boundaryPath = ""

    MyOcean2ROMS(gridParh, dataPath, ncepDate, initPath, boundaryPath)


if __name__ == '__main__':
    main()
