from scipy.interpolate import griddata
from concurrent.futures import ProcessPoolExecutor
from numpy.ma.core import MaskedArray
import multiprocessing
from numba import njit
import numpy as np


def gridInterp(srcLAT, srcLON, values, dstLAT, dstLON, fillValue, method):
    if isinstance(values, MaskedArray):
        values = values.filled(fill_value=np.nan)

    return griddata(
        (srcLAT.flatten(), srcLON.flatten()),
        values.flatten(),
        (dstLAT, dstLON),
        fill_value=fillValue,
        method=method
    )


def interp(srcLAT, srcLON, srcMASK, values, dstLAT, dstLON, fillValue):
    lon_src_mesh, lat_src_mesh = np.meshgrid(srcLON, srcLAT)

    valuesInterp = gridInterp(lat_src_mesh, lon_src_mesh, values, dstLAT, dstLON, fillValue)

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
    for j in range(dstSNDim):
        for i in range(dstWEDim):
            depthCopernicus = srcZ[:, 0, 0]
            srcCopernicus = tSrc[:, j, i]
            depthSigma = sigma[:, j, i]
            dstSigma = np.interp(depthSigma, depthCopernicus, srcCopernicus)
            tDst[:, j, i] = dstSigma


class Interpolator:
    def __init__(self, srcLAT, srcLON, dstLAT, dstLON, srcMASK, method="linear"):
        self.srcLAT = srcLAT
        self.srcLON = srcLON
        self.dstLAT = dstLAT
        self.dstLON = dstLON
        self.srcMASK = srcMASK
        self.method = method

        self.dstSNDim = len(dstLAT)
        self.dstWEDim = len(dstLAT[0])

    def interp(self, values, fillValue):
        interp_values = interp(self.srcLAT, self.srcLON, self.srcMASK, values, self.dstLAT, self.dstLON, fillValue)
        return interp_values

    def simpleInterp(self, values, fillValue):
        interp_values = gridInterp(self.srcLAT, self.srcLON, values, self.dstLAT, self.dstLON, fillValue, self.method)
        return interp_values


class BilinearInterpolator(Interpolator):
    def __init__(self, srcLAT, srcLON, dstLAT, dstLON, srcMASK, method="linear"):
        super().__init__(srcLAT, srcLON, dstLAT, dstLON, srcMASK, method)


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
