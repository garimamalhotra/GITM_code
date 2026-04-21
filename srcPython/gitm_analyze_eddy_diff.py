#!/usr/bin/env python3
"""Compare two GITM 3DALL runs as lat-lon maps (A | B | B-A) at fixed altitudes."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from gitm_routines import read_gitm_one_file

RUN_A = "/glade/u/home/garimam/gitm_03_26/GITM/run/UA/data.eddymap_const"
RUN_B = "/glade/u/home/garimam/gitm_03_26/GITM/run/UA/data.eddymap_varying"
OUT_DIR = "/glade/u/home/garimam/gitm_03_26/GITM/run/analysis_const_varying"
LABEL_A = "baseline (scalar Kzz=50)"
LABEL_B = "const Kzz(lat) map"
SNAPSHOTS = ["3DALL_t021221_000000.bin", "3DALL_t021221_000500.bin"]
ALT_INDICES = [22, 34]

VARS = [
    ("Tn",   "Temperature",   "inferno", False),
    ("Rho",  "Rho",           "viridis", True),
    ("O",    "[O(3P)]",       "viridis", True),
    ("O2",   "[O2]",          "viridis", True),
    ("N2",   "[N2]",          "viridis", True),
    ("Ne",   "[e-]",          "plasma",  True),
]


def load_run(path):
    head = read_gitm_one_file(path, vars_to_read=[0])
    nVars = head["nVars"]
    data = read_gitm_one_file(path, vars_to_read=list(range(nVars)))

    rad2deg = 180.0 / np.pi
    lon = data[0][:, 0, 0] * rad2deg
    lat = data[1][0, :, 0] * rad2deg
    alt = data[2][0, 0, :] / 1000.0

    result = {"lon": lon, "lat": lat, "alt": alt, "time": data["time"],
              "vars": data["vars"]}
    for short, name, _, _ in VARS:
        result[short] = data[data["vars"].index(name)]
    return result


def panel_map(ax, lon, lat, field, title, cmap, vmin=None, vmax=None,
              log=False, nlevels=30):
    field = field.T
    if log:
        field = np.where(field > 0, field, np.nan)
        if vmin is None:
            vmin = np.nanmin(field)
        if vmax is None:
            vmax = np.nanmax(field)
        levels = np.geomspace(vmin, vmax, nlevels)
        norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
        im = ax.contourf(lon, lat, field, levels=levels, cmap=cmap,
                         norm=norm, extend="both")
    else:
        if vmin is None:
            vmin = np.nanmin(field)
        if vmax is None:
            vmax = np.nanmax(field)
        levels = np.linspace(vmin, vmax, nlevels)
        im = ax.contourf(lon, lat, field, levels=levels, cmap=cmap,
                         vmin=vmin, vmax=vmax, extend="both")
    ax.set_xlabel("Longitude (deg)")
    ax.set_ylabel("Latitude (deg)")
    ax.set_xlim(lon.min(), lon.max())
    ax.set_ylim(-90, 90)
    ax.set_title(title)
    return im


def make_figure(short, cmap, log, runA, runB, k, out_path,
                snap_label, var_full_name):
    fA = runA[short]
    fB = runB[short]
    sA = fA[:, :, k]
    sB = fB[:, :, k]

    vmin = min(np.nanmin(sA), np.nanmin(sB))
    vmax = max(np.nanmax(sA), np.nanmax(sB))
    diff = sB - sA
    dmax = np.nanmax(np.abs(diff))
    if dmax == 0:
        dmax = 1.0

    fig, axs = plt.subplots(1, 3, figsize=(18, 4.5))
    imA = panel_map(axs[0], runA["lon"], runA["lat"], sA,
                    f"{LABEL_A}", cmap, vmin=vmin, vmax=vmax, log=log)
    imB = panel_map(axs[1], runB["lon"], runB["lat"], sB,
                    f"{LABEL_B}", cmap, vmin=vmin, vmax=vmax, log=log)
    imD = panel_map(axs[2], runB["lon"], runB["lat"], diff,
                    f"{LABEL_B} − {LABEL_A}", "RdBu_r",
                    vmin=-dmax, vmax=+dmax, log=False)

    fig.colorbar(imA, ax=axs[0], label=var_full_name)
    fig.colorbar(imB, ax=axs[1], label=var_full_name)
    fig.colorbar(imD, ax=axs[2], label="Δ " + var_full_name)

    fig.suptitle(f"{short} at z = {runA['alt'][k]:.0f} km   |   {snap_label}",
                 fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"wrote {out_path}")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    for snap in SNAPSHOTS:
        pathA = os.path.join(RUN_A, snap)
        pathB = os.path.join(RUN_B, snap)
        if not (os.path.exists(pathA) and os.path.exists(pathB)):
            print(f"skipping {snap}: missing in {RUN_A} or {RUN_B}")
            continue
        print(f"\n=== {snap} ===")
        runA = load_run(pathA)
        runB = load_run(pathB)
        snap_label = f"{snap}  ({runA['time']})"
        for short, name, cmap, log in VARS:
            for k in ALT_INDICES:
                alt_km = runA["alt"][k]
                tag = snap.replace("3DALL_", "").replace(".bin", "")
                fname = f"{tag}_{short}_z{int(alt_km)}km.png"
                make_figure(short, cmap, log, runA, runB, k,
                            os.path.join(OUT_DIR, fname),
                            snap_label, name)


if __name__ == "__main__":
    main()
