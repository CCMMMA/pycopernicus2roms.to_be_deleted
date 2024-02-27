from scipy.interpolate import griddata, interp1d
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
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
    def __init__(self, srcLAT, srcLON, srcZ, dstLAT, dstLON, dstZ, srcMASK, romsGrid):
        super().__init__(srcLAT, srcLON, dstLAT, dstLON, srcMASK)
        self.srcZ = srcZ
        self.dstZ = dstZ
        self.romsGrid = romsGrid

        self.srcLevs = 0
        self.dstLevs = 0
        self.sigma = None

        self.prepare()

    def prepare(self):
        self.srcLevs = len(self.srcZ)
        self.dstLevs = len(self.romsGrid.s_rho)
        self.sigma = abs(self.dstZ * self.romsGrid.s_rho[:, np.newaxis, np.newaxis])

    def interp(self, values, fillValue):
        tSrc = np.empty((self.srcLevs, self.dstSNDim, self.dstWEDim))
        tDst = np.empty((self.dstLevs, self.dstSNDim, self.dstWEDim))

        def parallel_interp(k):
            print(f"<k={k} depth:{self.srcZ[k][0][0]:.2f}>")
            result = super(type(self), self).interp(values[k], fillValue)
            # result = np.random.rand(self.dstSNDim, self.dstWEDim)  # JUST FOR TESTING
            return result

        num_processors = multiprocessing.cpu_count()
        print(f"Number of threads/processes used: {num_processors}")
        with ThreadPoolExecutor(max_workers=num_processors) as executor:
            futures = [executor.submit(parallel_interp, k) for k in range(self.srcLevs)]

            for k, future in enumerate(futures):
                tSrc[k] = future.result()

        print("Interpolating vertically...")
        for j in range(self.dstSNDim):
            for i in range(self.dstWEDim):
                depthCopernicus = self.srcZ[:, 0, 0].filled(np.nan)
                srcCopernicus = tSrc[:, j, i]
                interp = interp1d(depthCopernicus, srcCopernicus, kind='linear', fill_value="extrapolate")
                depthSigma = self.sigma[:, j, i].filled(np.nan)
                dstSigma = interp(depthSigma)
                tDst[:, j, i] = dstSigma

        return tDst
