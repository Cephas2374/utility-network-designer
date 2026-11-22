"""
Elevation from the SRTM 1-arc-second tile (big-endian int16, 3601 x 3601, rows north to south).
"""
import numpy as np
from pyproj import Transformer
from scipy.ndimage import uniform_filter

from .config import CRS


class SRTM:
    def __init__(self, path, lat_south=-5, lon_west=15, smooth=3):
        a = np.fromfile(path, dtype=">i2")
        n = int(round(np.sqrt(a.size)))
        z = a.reshape(n, n).astype(float)
        z[z < -1000] = np.nan
        if np.isnan(z).any():
            z[np.isnan(z)] = np.nanmean(z)
        # a 3x3 mean damps the surface (tree and roof) noise of the radar DEM
        self.z = uniform_filter(z, size=smooth) if smooth else z
        self.n = n
        self.lat_north = lat_south + 1
        self.lon_west = lon_west
        self.to_geo = Transformer.from_crs(CRS, 4326, always_xy=True)

    def sample(self, xy):
        """Bilinear elevation (m) at projected points xy (N x 2, EPSG:32733).

        The papers' surface. GridDEM below samples any raster the same way (3x3 mean, then bilinear
        between pixel centres), so that the QGIS plugin reproduces the papers on the same grid."""
        xy = np.atleast_2d(np.asarray(xy, dtype=float))
        lon, lat = self.to_geo.transform(xy[:, 0], xy[:, 1])
        c = (lon - self.lon_west) * (self.n - 1)
        r = (self.lat_north - lat) * (self.n - 1)
        c0, r0 = np.floor(c).astype(int), np.floor(r).astype(int)
        dc, dr = c - c0, r - r0
        z = self.z
        return ((1 - dr) * ((1 - dc) * z[r0, c0] + dc * z[r0, c0 + 1])
                + dr * ((1 - dc) * z[r0 + 1, c0] + dc * z[r0 + 1, c0 + 1]))


class GridDEM:
    """Elevation on any regular raster, sampled as SRTM.sample does: optional 3x3 mean, then bilinear
    interpolation between pixel centres.

    z       : 2-D array, row 0 at the top (north)
    x0, y0  : coordinates of the centre of pixel (0, 0) in the raster CRS
    dx, dy  : pixel size (dy < 0 for north-up rasters)
    to_grid : callable (x array, y array) -> (x array, y array) from the design CRS to the raster CRS,
              or None if they are the same
    nodata  : value treated as missing (filled with the mean, as for SRTM voids)

    Fix 2026-11-22: the plugin used to read the raw pixel under each vertex (QGIS provider.sample) of a
    DEM warped to UTM, while the papers sample the 3x3-smoothed SRTM bilinearly.  Ground levels differed
    by 0.65 m on average (7.2 m at most); on Gombe that turned the K=3 water design of the papers into
    K=4 and one fire check into a failure (77/78).  Sampling through this class on the native SRTM grid
    gives the papers' design (see papers/P5_software, Section 4).
    """

    def __init__(self, z, x0, y0, dx, dy, to_grid=None, nodata=None, smooth=3):
        z = np.asarray(z, dtype=float).copy()
        if nodata is not None:
            z[z == nodata] = np.nan
        z[z < -1000] = np.nan
        if np.isnan(z).all():
            raise ValueError("the elevation raster has no valid pixel")
        if np.isnan(z).any():
            z[np.isnan(z)] = np.nanmean(z)
        self.z = uniform_filter(z, size=smooth) if smooth and smooth > 1 else z
        self.x0, self.y0, self.dx, self.dy = float(x0), float(y0), float(dx), float(dy)
        self.to_grid = to_grid

    def sample(self, xy):
        xy = np.atleast_2d(np.asarray(xy, dtype=float))
        x, y = (self.to_grid(xy[:, 0], xy[:, 1]) if self.to_grid else (xy[:, 0], xy[:, 1]))
        c = (np.asarray(x) - self.x0) / self.dx
        r = (np.asarray(y) - self.y0) / self.dy
        nr, nc = self.z.shape
        c = np.clip(c, 0, nc - 1 - 1e-9)
        r = np.clip(r, 0, nr - 1 - 1e-9)
        c0, r0 = np.floor(c).astype(int), np.floor(r).astype(int)
        dc, dr = c - c0, r - r0
        c1, r1 = np.minimum(c0 + 1, nc - 1), np.minimum(r0 + 1, nr - 1)
        z = self.z
        return ((1 - dr) * ((1 - dc) * z[r0, c0] + dc * z[r0, c1])
                + dr * ((1 - dc) * z[r1, c0] + dc * z[r1, c1]))
