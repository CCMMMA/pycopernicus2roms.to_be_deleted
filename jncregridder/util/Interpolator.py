from scipy.interpolate import griddata
import numpy as np


class Interpolator:
    def __init__(self, srcLAT, srcLON, dstLAT, dstLON, srcMASK):
        self.srcLAT = srcLAT
        self.srcLON = srcLON
        self.dstLAT = dstLAT
        self.dstLON = dstLON
        self.srcMASK = srcMASK

        self.dstSNDim = len(dstLAT)
        self.dstWEDim = len(dstLAT[0])

    def interp(self, values, fillValue):
        lon_src_mesh, lat_src_mesh = np.meshgrid(self.srcLON, self.srcLAT)

        valuesInterp = griddata(
            (lat_src_mesh.flatten(), lon_src_mesh.flatten()),
            values.filled(fill_value=np.nan).flatten(),
            (self.dstLAT, self.dstLON),
            fill_value=fillValue,
            method="linear"
        )

        indices = np.where(self.srcMASK == 0)
        valuesInterp[indices] = 1e37

        rows, cols = np.where(~np.isnan(valuesInterp) & (valuesInterp != 1e37))
        values = valuesInterp[rows, cols]
        interp_rows, interp_cols = np.where(np.isnan(valuesInterp))
        interpolated_values = griddata((rows, cols), values, (interp_rows, interp_cols), method='nearest')
        valuesInterp[np.isnan(valuesInterp)] = interpolated_values

        indices = np.where(valuesInterp == 1e37)
        valuesInterp[indices] = np.nan

        return valuesInterp


class BilinearInterpolator(Interpolator):
    def __init__(self, srcLAT, srcLON, dstLAT, dstLON, srcMASK):
        super().__init__(srcLAT, srcLON, dstLAT, dstLON, srcMASK)


class BilinearInterpolator3D(Interpolator):
    def __init__(self, srcLAT, srcLON, srcZ, dstLAT, dstLON, dstZ, srcMASK):
        super().__init__(srcLAT, srcLON, dstLAT, dstLON, srcMASK)
        self.srcZ = srcZ
        self.dstZ = dstZ
        self.srcLevs = 0
        self.dstMinZ = 0
        self.weight3Ds = None

        self.prepare()

    def prepare(self):
        dstLevs = len(self.dstZ)
        self.dstMinZ = np.min(self.dstZ)
        # print("dstMinZ:" + str(self.dstMinZ))

        k = 0
        while k < len(self.srcZ) and abs(self.srcZ[k][0][0] <= abs(self.dstMinZ)):
            # print("srcDEPTH:" + str(self.srcZ[k][0][0]))
            k += 1
        self.srcLevs = k + 1
        print("srcLevs:", str(self.srcLevs))

    def interp(self, values, fillValue):
        tSrc = []
        for k in range(self.srcLevs):
            print(f"<k={k} depth:{self.srcZ[k][0][0]:.2f}")
            tSrc.append(super().interp(values[k], fillValue))
            print(" >")

        # TODO: add vertical interpolation

        return tSrc
