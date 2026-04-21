#!/usr/bin/env python3
"""Generate test NetCDF maps of the eddy diffusion coefficient Kzz at P0."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset

OUTPUT_DIR = "/glade/u/home/garimam/gitm_03_26/GITM/run"

NLAT, NLON = 73, 144

lat = np.linspace(-90.0, 90.0, NLAT, dtype="f4")
lon = np.linspace(0.0, 360.0, NLON, endpoint=False, dtype="f4")


def write_nc(path, kzz):
    with Dataset(path, "w", format="NETCDF4_CLASSIC") as nc:
        nc.createDimension("lat", NLAT)
        nc.createDimension("lon", NLON)
        v_lat = nc.createVariable("lat", "f4", ("lat",))
        v_lon = nc.createVariable("lon", "f4", ("lon",))
        v_kzz = nc.createVariable("Kzz_P0", "f4", ("lat", "lon"))
        v_lat[:] = lat
        v_lon[:] = lon
        v_kzz[:] = kzz
        v_lat.units = "degrees_north"
        v_lon.units = "degrees_east"
        v_kzz.units = "m^2 s^-1"
        v_kzz.long_name = "Eddy diffusion coefficient at homopause (P0)"
    print(f"wrote {path}: min={kzz.min():.1f} max={kzz.max():.1f}")


def plot_nc(nc_path, png_path=None, vmin=None, vmax=None):
    if png_path is None:
        png_path = nc_path.replace(".nc", ".png")

    with Dataset(nc_path, "r") as nc:
        latv = nc.variables["lat"][:]
        lonv = nc.variables["lon"][:]
        kzz  = nc.variables["Kzz_P0"][:]
        units = getattr(nc.variables["Kzz_P0"], "units", "")

    fig, ax = plt.subplots(figsize=(8, 4))
    im = ax.pcolormesh(lonv, latv, kzz, shading="auto",
                       cmap="viridis", vmin=vmin, vmax=vmax)
    ax.set_xlabel("Longitude (deg)")
    ax.set_ylabel("Latitude (deg)")
    ax.set_xlim(lonv.min(), lonv.max())
    ax.set_ylim(-90, 90)
    ax.set_title(f"{nc_path}: Kzz at P0")
    fig.colorbar(im, ax=ax, label=f"Kzz [{units}]")
    fig.tight_layout()
    fig.savefig(png_path, dpi=120)
    plt.close(fig)
    print(f"wrote {png_path}")


os.makedirs(OUTPUT_DIR, exist_ok=True)
const_path = os.path.join(OUTPUT_DIR, "eddy_map_const.nc")
var_path   = os.path.join(OUTPUT_DIR, "eddy_map_varying.nc")

kzz_const = np.full((NLAT, NLON), 50.0, dtype="f4")
write_nc(const_path, kzz_const)

LAT2D, _ = np.meshgrid(lat, lon, indexing="ij")
kzz_var = (25.0 + 75.0 * np.sin(np.deg2rad(LAT2D)) ** 2).astype("f4")
write_nc(var_path, kzz_var)

gmin = min(kzz_const.min(), kzz_var.min())
gmax = max(kzz_const.max(), kzz_var.max())
plot_nc(const_path, vmin=gmin, vmax=gmax)
plot_nc(var_path,   vmin=gmin, vmax=gmax)
