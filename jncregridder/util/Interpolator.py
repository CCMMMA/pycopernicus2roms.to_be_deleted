from scipy.interpolate import griddata, interp1d
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
from numba import njit
import numpy as np


def interp(srcLAT, srcLON, srcMASK, values, dstLAT, dstLON, fillValue):
    lon_src_mesh, lat_src_mesh = np.meshgrid(srcLON, srcLAT)

    valuesInterp = griddata(
        (lat_src_mesh.flatten(), lon_src_mesh.flatten()),
        values.filled(fill_value=np.nan).flatten(),
        (dstLAT, dstLON),
        fill_value=fillValue,
        method="linear"
    )

    indices = np.where(srcMASK == 0)
    valuesInterp[indices] = 1e37

    rows, cols = np.where(~np.isnan(valuesInterp) & (valuesInterp != 1e37))
    values = valuesInterp[rows, cols]
    interp_rows, interp_cols = np.where(np.isnan(valuesInterp))
    interpolated_values = griddata((rows, cols), values, (interp_rows, interp_cols), method='nearest')
    valuesInterp[np.isnan(valuesInterp)] = interpolated_values

    indices = np.where(valuesInterp == 1e37)
    valuesInterp[indices] = np.nan

    return valuesInterp


def interp_horizontal(k, srcLAT, srcLON, srcZ, srcMASK, values, dstLAT, dstLON, fillValue):
    print(f"<k={k} depth:{srcZ[k][0][0]:.2f}>")
    result = interp(srcLAT, srcLON, srcMASK, values[k], dstLAT, dstLON, fillValue)
    # result = np.random.rand(len(dstLAT), len(dstLAT[0]))
    return result


@njit
def interp_vertical(dstSNDim, dstWEDim, srcZ, tSrc, sigma, tDst):
    """
    depthCopernicus = srcZ[:, 0, 0].filled(np.nan)
    srcCopernicus = tSrc[:, j, i]
    interpV = interp1d(depthCopernicus, srcCopernicus, kind='linear', fill_value="extrapolate")
    depthSigma = sigma[:, j, i].filled(np.nan)
    dstSigma = interpV(depthSigma)
    tDst[:, j, i] = dstSigma
    """
    for j in range(dstSNDim):
        for i in range(dstWEDim):
            depthCopernicus = srcZ[:, 0, 0]
            srcCopernicus = tSrc[:, j, i]
            depthSigma = sigma[:, j, i]
            dstSigma = np.interp(depthSigma, depthCopernicus, srcCopernicus)
            tDst[:, j, i] = dstSigma


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
        interp_values = interp(self.srcLAT, self.srcLON, self.srcMASK, values, self.dstLAT, self.dstLON, fillValue)
        return interp_values


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

        num_processors = multiprocessing.cpu_count()
        print(f"Number of processes used: {num_processors}")
        with ProcessPoolExecutor(max_workers=num_processors) as executor:
            futures = [executor.submit(interp_horizontal, k, self.srcLAT, self.srcLON, self.srcZ, self.srcMASK,
                                       values, self.dstLAT, self.dstLON, fillValue) for k in range(self.srcLevs)]

            for k, future in enumerate(futures):
                tSrc[k] = future.result()

        print("Interpolating vertically...")
        interp_vertical(self.dstSNDim, self.dstWEDim, self.srcZ.filled(np.nan), tSrc, self.sigma.filled(np.nan), tDst)

        return tDst

