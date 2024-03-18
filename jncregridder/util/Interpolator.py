from scipy.interpolate import griddata, interp1d
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


def interp(srcLAT, srcLON, srcMASK, values, dstLAT, dstLON, fillValue, method="linear"):
    minIdxLat = np.argmin(np.abs(srcLAT - dstLAT.min())) - 1
    maxIdxLat = np.argmin(np.abs(srcLAT - dstLAT.max())) + 1
    minIdxLon = np.argmin(np.abs(srcLON - dstLON.min())) - 1
    maxIdxLon = np.argmin(np.abs(srcLON - dstLON.max())) + 1

    srcLAT = srcLAT[minIdxLat:maxIdxLat]
    srcLON = srcLON[minIdxLon:maxIdxLon]
    values = values[minIdxLat:maxIdxLat, minIdxLon:maxIdxLon]

    lon_src_mesh, lat_src_mesh = np.meshgrid(srcLON, srcLAT)

    valuesInterp = gridInterp(lat_src_mesh, lon_src_mesh, values, dstLAT, dstLON, np.nan, method)

    indices = np.where(srcMASK == 0)
    valuesInterp[indices] = 1e37

    rows, cols = np.where(~np.isnan(valuesInterp) & (valuesInterp != 1e37))
    values = valuesInterp[rows, cols]
    interp_rows, interp_cols = np.where(np.isnan(valuesInterp))
    interpolated_values = griddata((rows, cols), values, (interp_rows, interp_cols), fill_value=fillValue, method='nearest')
    valuesInterp[np.isnan(valuesInterp)] = interpolated_values

    # indices = np.where(valuesInterp == 1e37)
    # valuesInterp[indices] = np.nan

    return valuesInterp


def interp_horizontal(k, srcLAT, srcLON, srcZ, values, dstLAT, dstLON, dstZ, fillValue):
    print(f"<k={k} depth:{srcZ[k][0][0]:.2f}>")
    copernicusDepth = srcZ[k][0][0]

    lats = np.array([row[0] for row in dstLAT])
    lons = np.array(dstLON[0])

    @njit
    def createMask():
        mask = np.full((len(lats), len(lons)), 1)
        for j in range(len(lats)):
            for i in range(len(lons)):
                if dstZ[j][i] <= copernicusDepth:
                    mask[j][i] = 0

        return mask

    mask = createMask()
    result = interp(srcLAT, srcLON, mask, values[k], dstLAT, dstLON, fillValue)
    return result


# @njit
def interp_vertical(dstSNDim, dstWEDim, srcZ, tSrc, sigma, srcMask, indexes, tDst):
    depth = srcZ[:, 0, 0]
    for j in range(dstSNDim):
        for i in range(dstWEDim):
            if srcMask[j][i] == 1:
                depthCopernicus = depth[indexes[:, j, i]]
                srcCopernicus = tSrc[:, j, i][indexes[:, j, i]]
                depthSigma = sigma[:, j, i]

                f_interp = interp1d(depthCopernicus, srcCopernicus, kind='nearest', fill_value="extrapolate")
                dstSigma = f_interp(depthSigma)
                tDst[:, j, i] = dstSigma
            else:
                tDst[:, j, i] = 1e37


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
        self.maxK = 0
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

        for k in range(self.srcLevs):
            copernicusDepth = self.srcZ[k][0][0]

            if copernicusDepth > self.dstZ.max():
                self.maxK = k + 1
                break

    def interp(self, values, fillValue):
        tSrc = np.empty((self.srcLevs, self.dstSNDim, self.dstWEDim))
        tDst = np.empty((self.dstLevs, self.dstSNDim, self.dstWEDim))

        if isinstance(self.sigma, MaskedArray):
            self.sigma = self.sigma.filled(np.nan)
        if isinstance(self.srcZ, MaskedArray):
            self.srcZ = self.srcZ.filled(np.nan)
        if isinstance(self.dstZ, MaskedArray):
            self.dstZ = self.dstZ.filled(np.nan)
        if isinstance(self.srcMASK, MaskedArray):
            self.srcMASK = self.srcMASK.filled(np.nan)

        num_processors = multiprocessing.cpu_count()
        print(f"Number of processes used: {num_processors}")
        with ProcessPoolExecutor(max_workers=num_processors) as executor:
            futures = [executor.submit(interp_horizontal, k, self.srcLAT, self.srcLON, self.srcZ, values,
                                       self.dstLAT, self.dstLON, self.dstZ, fillValue) for k in range(self.maxK)]

            for k, future in enumerate(futures):
                tSrc[k] = future.result()

        print("Interpolating vertically...")
        @njit
        def createMatrix(values, dstLevs, dstSNDim, dstWEDim, sigma, srcZ):
            indexes = np.empty((dstLevs, dstSNDim, dstWEDim))
            for j in range(dstSNDim):
                for i in range(dstWEDim):
                    non_nan_indexes = np.where(np.isnan(values[:, j, i]))[0]
                    if len(non_nan_indexes) > 0 and non_nan_indexes[0] > 0:
                        distances = np.abs(sigma[:, j, i][:, np.newaxis] - srcZ[:non_nan_indexes[0], 0, 0])
                        indexes[:, j, i] = np.argmin(distances, axis=1)
                    else:
                        indexes[:, j, i] = 0

            return indexes
        indexes = np.floor(createMatrix(tSrc, self.dstLevs, self.dstSNDim, self.dstWEDim, self.sigma, self.srcZ)).astype(int)
        interp_vertical(self.dstSNDim, self.dstWEDim, self.srcZ, tSrc, self.sigma, self.srcMASK, indexes, tDst)

        return tDst
