# Eddy Diffusion Map

Optional 2D (lat, lon) map of Kzz at the homopause pressure P0. When turned
on, GITM reads the map from NetCDF at init and uses it instead of the scalar
EddyDiffusionCoef.

The NetCDF file needs dims lat, lon and variables lat, lon,
Kzz_P0(lat, lon) in m^2/s.

Enable it in UAM.in:

```
#EDDYMAP
T                       UseEddyMap
eddy_map_varying.nc     EddyMapFile
```

## Helper scripts

- srcPython/gitm_make_eddy_map.py    -- writes the NetCDF test maps.
- srcPython/gitm_analyze_eddy_diff.py -- compares two runs as 3-panel
  lat/lon plots (run A | run B | difference).

Paths inside the helper scripts are absolute -- edit OUTPUT_DIR, RUN_A,
RUN_B, OUT_DIR, ALT_INDICES at the top to match your setup.

## Workflow

1. Build GITM as usual.
2. Make test maps: python srcPython/gitm_make_eddy_map.py.
3. Run GITM
4. Post-process: python post_process.py.
5. Compare two runs: python srcPython/gitm_analyze_eddy_diff.py.

## Sanity checks

Three back-to-back runs from the same run directory. GITM always writes to
UA/data/, so after each run rename the output and recreate an empty
UA/data/ before starting the next.

1. Baseline -- set UseEddyMap F in UAM.in, run GITM, post-process:
   ```
   mpiexec -n <np> ./GITM.exe
   python post_process.py
   mv UA/data UA/data.eddymap_off
   mkdir UA/data
   ```

2. Uniform map -- set UseEddyMap T and EddyMapFile eddy_map_const.nc,
   rerun:
   ```
   mpiexec -n <np> ./GITM.exe
   python post_process.py
   mv UA/data UA/data.eddymap_const
   mkdir UA/data
   ```
   Compare against baseline by setting RUN_A/RUN_B in
   gitm_analyze_eddy_diff.py to data.eddymap_off and data.eddymap_const
   (and OUT_DIR to something unique, e.g. analysis_off_vs_const, so
   plots don't overwrite a previous comparison), then run it -- all three
   panels should look identical and the difference panel should be uniform
   zero.

3. Varying map -- set EddyMapFile eddy_map_varying.nc, rerun:
   ```
   mpiexec -n <np> ./GITM.exe
   python post_process.py
   mv UA/data UA/data.eddymap_varying
   mkdir UA/data
   ```
   Then set RUN_A/RUN_B in gitm_analyze_eddy_diff.py to
   data.eddymap_off and data.eddymap_varying (and give OUT_DIR a
   unique name, e.g. analysis_off_vs_varying), then run it. Expect the
   largest differences near the poles (where Kzz = 100 vs. 25 at the
   equator).
