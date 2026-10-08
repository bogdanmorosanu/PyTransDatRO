#!/usr/bin/env python3
"""
Extract and report metadata from Romgeo SPG grid files (.spg).

This tool reads the compiled 3D binary grid file (.spg) located in the 
pytransdatro/grids/ directory, extracts all embedded metadata, parameters, 
and grid layer properties, and generates a structured text report.
"""

import os
import sys
import glob
import datetime
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

try:
    import numpy as np
except ImportError:
    print("Error: numpy is required to read and unpickle .spg files.")
    sys.exit(1)

from pytransdatro.utils import rad_to_sexa, DEG_TO_RAD_FACTOR


def decdeg_to_dms(dd: float, is_lat: bool = True) -> str:
    """Convert decimal degrees to human-readable DMS string using rad_to_sexa from utils."""
    rad_val = dd * DEG_TO_RAD_FACTOR
    sexa_str = rad_to_sexa(rad_val)
    parts = sexa_str.split(' ')
    deg = int(parts[0])
    minute = int(parts[1])
    second = float(parts[2])
    
    sign = -1 if dd < 0 else 1
    if is_lat:
        hemi = "N" if sign >= 0 else "S"
    else:
        hemi = "E" if sign >= 0 else "W"
    return f"{abs(deg)}° {minute:02d}' {second:05.2f}\" {hemi} ({dd:.8f}°)"


def extract_metadata(spg_path: Path, output_report_path: Path):
    """Extract metadata from an SPG file and write an annotated report."""
    if not spg_path.exists():
        raise FileNotFoundError(f"SPG grid file not found at: {spg_path}")

    # Load SPG data dictionary
    raw_data = np.load(str(spg_path), allow_pickle=True)
    if not isinstance(raw_data, dict):
        if hasattr(raw_data, "item"):
            data = raw_data.item()
        else:
            raise ValueError(f"Unexpected data format in {spg_path}: {type(raw_data)}")
    else:
        data = raw_data

    params = data.get("params", {})
    grids = data.get("grids", {})
    metadata = data.get("metadata", {})

    file_size_bytes = spg_path.stat().st_size
    file_size_kb = file_size_bytes / 1024.0

    lines = []
    def add_line(text=""):
        lines.append(text)

    def add_section(title):
        add_line("=" * 80)
        add_line(title.upper())
        add_line("=" * 80)

    # ---------------------------------------------------------
    # Header
    # ---------------------------------------------------------
    add_section("Romgeo SPG Grid Metadata Report")
    add_line(f"Generated On        : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    add_line(f"Source File         : {spg_path.name}")
    add_line(f"File Path           : {spg_path}")
    add_line(f"File Size           : {file_size_kb:.2f} KB ({file_size_bytes:,} bytes)")
    add_line()

    # ---------------------------------------------------------
    # 1. General Metadata & Provenance
    # ---------------------------------------------------------
    add_section("1. General Metadata & Provenance")
    add_line(f"Grid File Name      : {metadata.get('file', 'N/A')}")
    add_line(f"Author / Creator    : {metadata.get('created_by', 'N/A')}")
    add_line(f"License             : {metadata.get('license', 'N/A')}")
    add_line(f"Abstract / Copyright: {metadata.get('abstract', 'N/A')}")
    
    release = metadata.get("release", {})
    if release:
        add_line(f"Release Version     : {release.get('major', 'N/A')}.{release.get('minor', 'N/A'):02d} (Revision: {release.get('revision', 0)}, Legacy: {release.get('legacy', 'no')})")
    add_line(f"Release Date        : {metadata.get('release_date', 'N/A')}")
    add_line(f"Validity Period     : From {metadata.get('valid_from', 'N/A')} to {metadata.get('valid_to', 'Indefinite (None)')}")
    
    notes = metadata.get("notes", "")
    if notes:
        add_line(f"Release Notes       : {notes}")
    add_line()

    attribution = metadata.get("attribution", "")
    if attribution:
        add_line("Attribution HTML    :")
        add_line(f"  {attribution}")
    add_line()

    # ---------------------------------------------------------
    # 2. Build & Interpolation Parameters (params)
    # ---------------------------------------------------------
    add_section("2. Build Parameters & Interpolation Methods")
    add_line(f"Version String      : {params.get('version', 'N/A')}")
    add_line(f"Target Output File  : {params.get('output_file', 'N/A')}")
    add_line(f"Geodetic Shifts Src : {params.get('geodetic_shifts_file', 'N/A')}")
    add_line(f"Geoid Heights Src   : {params.get('geoid_heights_file', 'N/A')}")
    add_line()

    # Interpolation flags
    interp = params.get("interpolation", {})
    h_interp = interp.get("horizontal", "N/A")
    v_interp = interp.get("vertical", "N/A")

    method_names = {
        0: "INTERP_COLOCATE (Nearest-Neighbor)",
        1: "INTERP_BILINEAR (Bilinear)",
        2: "INTERP_BICUBIC (Bicubic Spline)"
    }
    h_desc = method_names.get(h_interp, f"Method {h_interp}")
    v_desc = method_names.get(v_interp, f"Method {v_interp}")

    add_line("Interpolation Strategy:")
    add_line(f"  - Horizontal (2D) : {h_interp} -> {h_desc}")
    add_line(f"  - Vertical (1D)   : {v_interp} -> {v_desc}")
    add_line()

    # Helmert Parameters
    helmert = params.get("helmert", {})
    add_line("2D Helmert Transformation Parameters:")
    add_line()

    for direction, h_params in helmert.items():
        add_line(f"  Direction: [{direction}]")
        tE = h_params.get("tE", 0.0)
        tN = h_params.get("tN", 0.0)
        dm = h_params.get("dm", 0.0)
        Rz = h_params.get("Rz", 0.0)
        add_line(f"    - Translation East (tE)  : {tE:+.7f} m")
        add_line(f"    - Translation North (tN) : {tN:+.7f} m")
        add_line(f"    - Scale Correction (dm)  : {dm:+.8f} ppm")
        add_line(f"    - Z-Axis Rotation (Rz)   : {Rz:+.8f} arcseconds")
    add_line()

    # ---------------------------------------------------------
    # 3. Layer: Geodetic Shifts (2D Planimetric Grid)
    # ---------------------------------------------------------
    add_section("3. Grid Layer: Geodetic Shifts (2D Horizontal Distortion)")
    geo_shifts = grids.get("geodetic_shifts", {})
    geo_meta = geo_shifts.get("metadata", {})
    geo_grid = geo_shifts.get("grid")

    add_line(f"Layer Name          : {geo_shifts.get('name', 'N/A')}")
    add_line(f"Source Coordinate   : {geo_shifts.get('source', 'N/A')}")
    add_line(f"Target Coordinate   : {geo_shifts.get('target', 'N/A')}")
    add_line(f"CRS Type            : {geo_meta.get('crs_type', 'N/A')}")
    add_line(f"Dimension Count     : {geo_meta.get('ndim', 'N/A')} dimensions")
    add_line()
    add_line("Planimetric Extents:")
    add_line(f"  - Easting  (E) Min : {geo_meta.get('mine', 0.0):12.3f} m")
    add_line(f"  - Easting  (E) Max : {geo_meta.get('maxe', 0.0):12.3f} m (Span: {geo_meta.get('maxe', 0.0) - geo_meta.get('mine', 0.0):,.1f} m)")
    add_line(f"  - Northing (N) Min : {geo_meta.get('minn', 0.0):12.3f} m")
    add_line(f"  - Northing (N) Max : {geo_meta.get('maxn', 0.0):12.3f} m (Span: {geo_meta.get('maxn', 0.0) - geo_meta.get('minn', 0.0):,.1f} m)")
    add_line()
    add_line("Grid Resolution & Dimensions:")
    add_line(f"  - Step Easting (E) : {geo_meta.get('stepe', 0.0):12.3f} m")
    add_line(f"  - Step Northing(N) : {geo_meta.get('stepn', 0.0):12.3f} m")
    ncols_geo = geo_meta.get('ncols', 0)
    nrows_geo = geo_meta.get('nrows', 0)
    add_line(f"  - Columns (East)   : {ncols_geo}")
    add_line(f"  - Rows (North)     : {nrows_geo}")
    add_line(f"  - Total Grid Nodes : {ncols_geo * nrows_geo:,} nodes")
    
    if isinstance(geo_grid, np.ndarray):
        add_line()
        add_line("Data Array Properties:")
        add_line(f"  - Array Shape      : {geo_grid.shape} (Bands, Rows, Columns)")
        add_line(f"  - Data Type        : {geo_grid.dtype}")
        add_line(f"  - Memory Usage     : {geo_grid.nbytes / 1024.0:.2f} KB")

        band_n = geo_grid[0]
        band_e = geo_grid[1]
        valid_n = band_n[~np.isnan(band_n)]
        valid_e = band_e[~np.isnan(band_e)]

        add_line()
        add_line("Statistical Distribution:")
        add_line(f"  - Valid Nodes      : {len(valid_n):,} / {ncols_geo * nrows_geo:,} ({len(valid_n) / (ncols_geo * nrows_geo) * 100:.1f}%)")
        add_line(f"  - NoData / NaN     : {np.isnan(band_n).sum():,} nodes")
        if len(valid_n) > 0:
            add_line(f"  - Northing Shift dN: Min = {valid_n.min():+.4f} m, Max = {valid_n.max():+.4f} m, Mean = {valid_n.mean():+.4f} m, Std = {valid_n.std():.4f} m")
        if len(valid_e) > 0:
            add_line(f"  - Easting Shift dE : Min = {valid_e.min():+.4f} m, Max = {valid_e.max():+.4f} m, Mean = {valid_e.mean():+.4f} m, Std = {valid_e.std():.4f} m")
    add_line()

    # ---------------------------------------------------------
    # 4. Layer: Geoid Heights (1D Vertical Quasigeoid Grid)
    # ---------------------------------------------------------
    add_section("4. Grid Layer: Geoid Heights (1D Vertical Quasigeoid Model)")
    geoid = grids.get("geoid_heights", {})
    geoid_meta = geoid.get("metadata", {})
    geoid_grid = geoid.get("grid")

    add_line(f"Layer Name          : {geoid.get('name', 'N/A')}")
    add_line(f"Source Coordinate   : {geoid.get('source', 'N/A')}")
    add_line(f"Target Coordinate   : {geoid.get('target', 'N/A')}")
    add_line(f"CRS Type            : {geoid_meta.get('crs_type', 'N/A')}")
    add_line(f"Dimension Count     : {geoid_meta.get('ndim', 'N/A')} dimension")
    add_line()

    minla = geoid_meta.get('minla', 0.0)
    maxla = geoid_meta.get('maxla', 0.0)
    minphi = geoid_meta.get('minphi', 0.0)
    maxphi = geoid_meta.get('maxphi', 0.0)
    stepla = geoid_meta.get('stepla', 0.0)
    stepphi = geoid_meta.get('stepphi', 0.0)
    ncols_geoid = geoid_meta.get('ncols', 0)
    nrows_geoid = geoid_meta.get('nrows', 0)

    add_line("Geodetic Extents:")
    add_line(f"  - Longitude (Lambda) Min : {decdeg_to_dms(minla, is_lat=False)}")
    add_line(f"  - Longitude (Lambda) Max : {decdeg_to_dms(maxla, is_lat=False)}")
    add_line(f"    * Longitude Span       : {maxla - minla:.6f}°")
    add_line(f"  - Latitude  (Phi)    Min : {decdeg_to_dms(minphi, is_lat=True)}")
    add_line(f"  - Latitude  (Phi)    Max : {decdeg_to_dms(maxphi, is_lat=True)}")
    add_line(f"    * Latitude Span        : {maxphi - minphi:.6f}°")
    add_line()
    add_line("Grid Resolution & Dimensions:")
    add_line(f"  - Step Longitude (dLambda): {stepla:.6f}° ({stepla * 60.0:.2f} arcminutes)")
    add_line(f"  - Step Latitude  (dPhi)   : {stepphi:.6f}° ({stepphi * 60.0:.2f} arcminutes)")
    add_line(f"  - Columns (Longitude)     : {ncols_geoid}")
    add_line(f"  - Rows (Latitude)         : {nrows_geoid}")
    add_line(f"  - Total Grid Nodes        : {ncols_geoid * nrows_geoid:,} nodes")

    if isinstance(geoid_grid, np.ndarray):
        add_line()
        add_line("Data Array Properties:")
        add_line(f"  - Array Shape      : {geoid_grid.shape} (Bands, Rows, Columns)")
        add_line(f"  - Data Type        : {geoid_grid.dtype}")
        add_line(f"  - Memory Usage     : {geoid_grid.nbytes / 1024.0:.2f} KB")

        h_vals = geoid_grid[0]
        valid_h = h_vals[~np.isnan(h_vals)]

        add_line()
        add_line("Statistical Distribution:")
        add_line(f"  - Valid Nodes      : {len(valid_h):,} / {ncols_geoid * nrows_geoid:,} ({len(valid_h) / (ncols_geoid * nrows_geoid) * 100:.1f}%)")
        add_line(f"  - NoData / NaN     : {np.isnan(h_vals).sum():,} nodes")
        if len(valid_h) > 0:
            add_line(f"  - Quasigeoid Height: Min = {valid_h.min():.4f} m, Max = {valid_h.max():.4f} m, Mean = {valid_h.mean():.4f} m, Std = {valid_h.std():.4f} m")
    add_line()

    # ---------------------------------------------------------
    # 5. Transformation Pipeline Integration Summary
    # ---------------------------------------------------------
    add_section("5. Transformation Pipeline Integration Summary")
    add_line("How PyTransDatRO applies these metadata parameters:")
    add_line()
    add_line(f"1. Horizontal (2D) Shifts: {h_desc}")
    add_line(f"2. Vertical (1D) Heights : {v_desc}")
    add_line(f"3. 2D Helmert Constants  : Dynamically loaded from params['helmert']")
    add_line()
    add_line("=" * 80)
    add_line("End of Report")
    add_line("=" * 80)

    # Write output report
    output_report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    print(f"Metadata report successfully created at: {output_report_path}")
    return output_report_path


def main():
    grids_dir = PROJECT_ROOT / "pytransdatro" / "grids"
    spg_candidates = list(grids_dir.glob("*.spg"))

    if len(sys.argv) > 1:
        spg_file = Path(sys.argv[1])
    elif len(spg_candidates) == 1:
        spg_file = spg_candidates[0]
    else:
        spg_file = grids_dir / "rom_grid3d_25.09.spg"

    default_report = PROJECT_ROOT / "reports" / "spg_grid_metadata.txt"
    report_file = Path(sys.argv[2]) if len(sys.argv) > 2 else default_report

    extract_metadata(spg_file, report_file)


if __name__ == "__main__":
    main()
