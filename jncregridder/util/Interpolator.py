from scipy.interpolate import griddata, interp1d
from concurrent.futures import ProcessPoolExecutor
from numpy.ma.core import MaskedArray
from numba import jit, types
from numba.experimental import jitclass
import numpy as np
import multiprocessing


@jitclass([('w', types.float64[:]),
           ('KK', types.int32[:]),
           ('masked', types.boolean)])
class Weight3DStruct:
    def __init__(self):
        self.w = np.zeros(2, dtype=np.float64)
        self.KK = np.zeros(2, dtype=np.int32)
        self.masked = False


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


@jit(nopython=True)
def haversine(lat1, lon1, lat2, lon2):
    # convert decimal degrees to radians
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0) ** 2
    c = 2 * np.arcsin(np.sqrt(a))
    km = 6367 * c
    return km


@jit(nopython=True)
def interp(src, srcMissingValue, dstMissingValue, dstLAT, dstLON, srcLAT, srcLON, dstMASK):
    USE_IDW = True

    dstEta = len(dstLAT)
    dstXi = len(dstLAT[0])
    srcEta = len(srcLAT)
    srcXi = len(srcLAT[0])

    srcLonMin = srcLON[0][0]
    srcLatMin = srcLAT[0][0]
    srcLonMax = srcLON[srcEta - 1][srcXi - 1]
    srcLatMax = srcLAT[srcEta - 1][srcXi - 1]
    srcLatDelta = srcLatMax - srcLatMin
    srcLonDelta = srcLonMax - srcLonMin
    srcLatStep = srcLatDelta / srcEta
    srcLonStep = srcLonDelta / srcXi

    dst = np.full((dstEta, dstXi), dstMissingValue)

    for dstJ in range(dstEta):
        for dstI in range(dstXi):
            if dstMASK[dstJ][dstI] == 1:
                dstLon = dstLON[dstJ][dstI]
                dstLat = dstLAT[dstJ][dstI]

                srcII = (dstLon - srcLonMin) / srcLonStep
                srcJJ = (dstLat - srcLatMin) / srcLatStep

                iR = int(srcII)
                jR = int(srcJJ)

                pointsBilinear = []

                for j in [jR - 1, jR + 1]:
                    for i in [iR - 1, iR + 1]:
                        jj = min(max(j, 0), srcEta - 1)
                        ii = min(max(i, 0), srcXi - 1)
                        if not np.isnan(src[jj][ii]) and src[jj][ii] != srcMissingValue:
                            pointsBilinear.append((jj, ii, src[jj][ii]))

                if len(pointsBilinear) == 4:
                    lon1 = srcLON[pointsBilinear[0][0]][pointsBilinear[0][1]]
                    lat1 = srcLAT[pointsBilinear[0][0]][pointsBilinear[0][1]]
                    lon2 = srcLON[pointsBilinear[1][0]][pointsBilinear[1][1]]
                    lat2 = srcLAT[pointsBilinear[2][0]][pointsBilinear[2][1]]

                    dLon = lon2 - lon1
                    dLat = lat2 - lat1

                    FXY1 = ((lon2 - dstLon) / dLon) * pointsBilinear[0][2] + ((dstLon - lon1) / dLon) * \
                           pointsBilinear[1][2]
                    FXY2 = ((lon2 - dstLon) / dLon) * pointsBilinear[2][2] + ((dstLon - lon1) / dLon) * \
                           pointsBilinear[3][2]

                    dst[dstJ][dstI] = ((lat2 - dstLat) / dLat) * FXY1 + ((dstLat - lat1) / dLat) * FXY2

                elif USE_IDW:
                    pointsIDW = []
                    size = 0

                    while len(pointsIDW) == 0 and size < 4:
                        size += 1
                        for j in range(jR - size, jR + size + 1):
                            jj = min(max(j, 0), srcEta - 1)
                            for i in range(iR - size, iR + size + 1):
                                ii = min(max(i, 0), srcXi - 1)
                                if not np.isnan(src[jj][ii]) and src[jj][ii] != srcMissingValue:
                                    pointsIDW.append(
                                        (jj, ii, haversine(
                                            dstLat, dstLon,
                                            srcLAT[jj][ii], srcLON[jj][ii]
                                        ), src[jj][ii])
                                    )

                    if len(pointsIDW) > 0:
                        weighted_values_sum = 0.0
                        sum_of_weights = 0.0
                        for point in pointsIDW:
                            weight = 1 / point[2]
                            # weight = 1 / (point[2] ** 2)
                            sum_of_weights += weight
                            weighted_values_sum += weight * point[3]
                        dst[dstJ][dstI] = weighted_values_sum / sum_of_weights
    return dst


def interp_horizontal(k, srcLAT, srcLON, srcZ, values, dstLAT, dstLON, dstMask, fillValue):
    print(f"<k={k} depth:{srcZ[k][0][0]:.2f}>")

    result = interp(values[k].filled(np.nan), fillValue, 1e37, dstLAT, dstLON, srcLAT, srcLON, dstMask)
    return result


@jit(nopython=True)
def vertical_interp(dstLevs, srcEta, srcXi, dstMASK, weight3Ds, tSrc):
    tDst = np.empty((dstLevs, srcEta, srcXi))
    idx = 0
    for dstK in range(dstLevs):
        for dstJ in range(srcEta):
            for dstI in range(srcXi):
                if dstMASK[dstJ][dstI] == 1:
                    srcI = dstI
                    srcJ = dstJ

                    srcK0 = weight3Ds[idx].KK[0]
                    srcK1 = weight3Ds[idx].KK[1]
                    srcW0 = weight3Ds[idx].w[0]
                    srcW1 = weight3Ds[idx].w[1]
                    tSrcK0 = tSrc[srcK0][srcJ][srcI]
                    tSrcK1 = tSrc[srcK1][srcJ][srcI]

                    if tSrcK0 != 1e37 and tSrcK1 != 1e37:
                        tDst[dstK][dstJ][dstI] = tSrcK0 * srcW0 + tSrcK1 * srcW1
                    elif tSrcK0 != 1e37:
                        tDst[dstK][dstJ][dstI] = tSrcK0
                    elif tSrcK1 != 1e37:
                        tDst[dstK][dstJ][dstI] = tSrcK1
                    else:
                        tDst[dstK][dstJ][dstI] = 1e37
                else:
                    tDst[dstK][dstJ][dstI] = 1e37

                idx += 1

    return tDst


class Interpolator:
    def __init__(self, srcLAT, srcLON, dstLAT, dstLON, dstMASK, method="linear"):
        self.srcLAT = srcLAT
        self.srcLON = srcLON
        self.dstLAT = dstLAT.filled(np.nan)
        self.dstLON = dstLON.filled(np.nan)
        self.dstMASK = dstMASK.filled(np.nan)
        self.method = method

        self.dstSNDim = len(dstLAT)
        self.dstWEDim = len(dstLAT[0])

    def interp(self, values, fillValue):
        interp_values = interp(values.filled(np.nan), fillValue, 1e37, self.dstLAT, self.dstLON, self.srcLAT,
                               self.srcLON, self.dstMASK)
        return interp_values

    def simpleInterp(self, values, fillValue):
        interp_values = gridInterp(self.srcLAT, self.srcLON, values, self.dstLAT, self.dstLON, fillValue, self.method)
        return interp_values


class BilinearInterpolator(Interpolator):
    def __init__(self, srcLAT, srcLON, dstLAT, dstLON, dstMASK, method="linear"):
        super().__init__(srcLAT, srcLON, dstLAT, dstLON, dstMASK, method)


class BilinearInterpolator3D(Interpolator):
    def __init__(self, srcLAT, srcLON, srcZ, dstLAT, dstLON, dstZ, dstMASK, romsGrid):
        super().__init__(srcLAT, srcLON, dstLAT, dstLON, dstMASK)
        self.srcZ = srcZ.filled(np.nan)
        self.dstZ = dstZ.filled(np.nan)
        self.srcLevs = len(srcZ)
        self.dstLevs = len(romsGrid.s_rho)

        self.maxK = 0
        for k in range(self.srcLevs):
            copernicusDepth = self.srcZ[k][0][0]

            if copernicusDepth > romsGrid.H.max():
                self.maxK = k + 1
                break

        self.weight3Ds = self.prepare(self.maxK, self.dstLevs, self.srcZ, self.dstZ, self.dstSNDim, self.dstWEDim)

    @staticmethod
    @jit(nopython=True)
    def prepare(srcLevs, dstLevs, srcZ, dstZ, dstSNDim, dstWEDim):
        weight3Ds = []
        for dstK in range(dstLevs):
            for dstJ in range(dstSNDim):
                for dstI in range(dstWEDim):
                    dstZatKJI = abs(dstZ[dstK][dstJ][dstI])

                    srcK = 0
                    srcZat00 = 0.0
                    while srcK < srcLevs:
                        srcZat00 = abs(srcZ[srcK][0][0])
                        if srcZat00 > abs(dstZatKJI):
                            break
                        srcK += 1

                    srcKmin = srcK - 1
                    srcKmax = srcK

                    if srcKmax == srcLevs:
                        srcKmax = srcKmin

                    if srcKmin < 0:
                        srcKmin = 0

                    srcZmin = abs(srcZ[srcKmin][0][0])
                    srcZmax = abs(srcZ[srcKmax][0][0])
                    delta = srcZmax - srcZmin

                    weight3D = Weight3DStruct()

                    if delta != 0.0:
                        weight3D.w[0] = (dstZatKJI - srcZmin) / delta
                        weight3D.w[1] = (srcZmax - dstZatKJI) / delta
                    else:
                        weight3D.w[0] = 0.0
                        weight3D.w[1] = 1.0

                    weight3D.KK[0] = srcKmin
                    weight3D.KK[1] = srcKmax

                    weight3Ds.append(weight3D)

        return weight3Ds

    def interp(self, values, fillValue):
        tSrc = np.empty((self.srcLevs, self.dstSNDim, self.dstWEDim))

        if isinstance(self.srcZ, MaskedArray):
            self.srcZ = self.srcZ.filled(np.nan)
        if isinstance(self.dstZ, MaskedArray):
            self.dstZ = self.dstZ.filled(np.nan)
        if isinstance(self.dstMASK, MaskedArray):
            self.dstMASK = self.dstMASK.filled(np.nan)

        num_processors = multiprocessing.cpu_count()
        print(f"Number of processes used: {num_processors}")
        with ProcessPoolExecutor(max_workers=num_processors) as executor:
            futures = [executor.submit(interp_horizontal, k, self.srcLAT, self.srcLON, self.srcZ, values,
                                       self.dstLAT, self.dstLON, self.dstMASK, fillValue) for k in range(self.maxK)]

            for k, future in enumerate(futures):
                tSrc[k] = future.result()

        print("Interpolating vertically...")
        tDst = vertical_interp(self.dstLevs, self.dstSNDim, self.dstWEDim, self.dstMASK, self.weight3Ds, tSrc)

        return tDst
