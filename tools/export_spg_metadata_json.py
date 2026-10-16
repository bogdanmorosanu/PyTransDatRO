#!/usr/bin/env python3
"""
Dynamic and resilient SPG grid metadata extractor.

Reads binary 3D SPG grid files (.spg), exhaustively extracts 100% of all
embedded parameters, metadata, and grid layer properties (dynamically adapting
to new, altered, or removed attributes in future releases), computes layer
diagnostics and statistics, and serializes the complete dataset into a
clean, web-ready JSON file.
"""

import datetime
import hashlib
import json
import sys
from pathlib import Path
from typing import Any, Dict

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

try:
    import numpy as np
except ImportError:
    print("Error: numpy is required to read and process .spg files.")
    sys.exit(1)


def compute_sha256(file_path: Path) -> str:
    """Computes SHA-256 hash of a file."""
    h = hashlib.sha256()
    with open(file_path, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def format_bytes(num_bytes: int) -> str:
    """Formats raw bytes into a human-readable size string."""
    for unit in ["B", "KB", "MB", "GB"]:
        if abs(num_bytes) < 1024.0:
            return f"{num_bytes:3.2f} {unit}"
        num_bytes /= 1024.0
    return f"{num_bytes:.2f} TB"


def decdeg_to_dms(dd: float, is_lat: bool = True) -> str:
    """Converts decimal degrees to human-readable DMS string with proper rollover."""
    try:
        sign = -1 if dd < 0 else 1
        abs_dd = abs(dd)
        deg = int(abs_dd)
        rem = (abs_dd - deg) * 60.0
        minute = int(rem)
        second = (rem - minute) * 60.0

        # Handle floating point rounding carrying over to minutes/degrees
        second = round(second, 2)
        if second >= 60.0:
            second -= 60.0
            minute += 1
        if minute >= 60:
            minute -= 60
            deg += 1

        hemi = ("N" if sign >= 0 else "S") if is_lat else ("E" if sign >= 0 else "W")
        return f"{deg}° {minute:02d}' {second:05.2f}\" {hemi}"
    except Exception:
        return f"{dd:.8f}°"


def normalize_value(val: Any) -> Any:
    """Recursively converts NumPy and arbitrary Python objects into standard JSON-serializable types."""
    if isinstance(val, dict):
        return {str(k): normalize_value(v) for k, v in val.items()}
    elif isinstance(val, (list, tuple)):
        return [normalize_value(v) for v in val]
    elif isinstance(val, np.integer):
        return int(val)
    elif isinstance(val, np.floating):
        if np.isnan(val) or np.isinf(val):
            return None
        return float(val)
    elif isinstance(val, np.ndarray):
        # Do not dump full raw arrays directly into metadata JSON; summarize them
        return {"shape": list(val.shape), "dtype": str(val.dtype), "size": int(val.size)}
    elif isinstance(val, (datetime.date, datetime.datetime)):
        return val.isoformat()
    return val


def analyze_layer(layer_name: str, layer_dict: Dict[str, Any]) -> Dict[str, Any]:
    """Dynamically inspects a grid layer and computes diagnostics and statistics without hardcoding schemas."""
    layer_summary: Dict[str, Any] = {}

    # 1. Copy all non-grid properties dynamically (preserving unknown/future fields)
    for k, v in layer_dict.items():
        if k != "grid":
            layer_summary[k] = normalize_value(v)

    # 2. Inspect the grid array if present
    grid_data = layer_dict.get("grid")
    if isinstance(grid_data, np.ndarray):
        shape = list(grid_data.shape)
        dtype_str = str(grid_data.dtype)
        memory_bytes = int(grid_data.nbytes)

        total_nodes = int(grid_data.size)
        # If multidimensional, determine row/col size per band
        num_bands = shape[0] if len(shape) == 3 else 1
        nodes_per_band = int(total_nodes / num_bands) if num_bands > 0 else total_nodes

        layer_summary["array_properties"] = {
            "shape": shape,
            "dimensions_count": len(shape),
            "num_bands": num_bands,
            "dtype": dtype_str,
            "memory_bytes": memory_bytes,
            "memory_formatted": format_bytes(memory_bytes),
            "total_nodes": total_nodes,
            "nodes_per_band": nodes_per_band,
        }

        # Compute statistics per band
        band_stats = []
        for b_idx in range(num_bands):
            band_arr = grid_data[b_idx] if len(shape) == 3 else grid_data
            valid_mask = ~np.isnan(band_arr)
            valid_count = int(np.count_nonzero(valid_mask))
            nan_count = int(band_arr.size - valid_count)
            valid_pct = (
                round((valid_count / band_arr.size) * 100.0, 2) if band_arr.size > 0 else 0.0
            )

            valid_values = band_arr[valid_mask]
            stats = {
                "band_index": b_idx,
                "total_nodes": int(band_arr.size),
                "valid_nodes": valid_count,
                "nan_nodes": nan_count,
                "valid_percentage": valid_pct,
            }

            if len(valid_values) > 0:
                stats.update(
                    {
                        "min": round(float(valid_values.min()), 6),
                        "max": round(float(valid_values.max()), 6),
                        "mean": round(float(valid_values.mean()), 6),
                        "std": round(float(valid_values.std()), 6),
                    }
                )
            else:
                stats.update(
                    {
                        "min": None,
                        "max": None,
                        "mean": None,
                        "std": None,
                    }
                )
            band_stats.append(stats)

        layer_summary["band_statistics"] = band_stats

        # Human-friendly extents and resolution formatting if metadata keys exist
        meta = layer_dict.get("metadata", {})
        if isinstance(meta, dict):
            # Check for geodetic fields (minphi, minla, etc.)
            if "minphi" in meta and "minla" in meta:
                minla = float(meta.get("minla", 0.0))
                maxla = float(meta.get("maxla", 0.0))
                minphi = float(meta.get("minphi", 0.0))
                maxphi = float(meta.get("maxphi", 0.0))
                stepla = float(meta.get("stepla", 0.0))
                stepphi = float(meta.get("stepphi", 0.0))

                layer_summary["formatted_extents"] = {
                    "longitude": {
                        "min_deg": minla,
                        "max_deg": maxla,
                        "span_deg": round(maxla - minla, 8),
                        "min_dms": decdeg_to_dms(minla, is_lat=False),
                        "max_dms": decdeg_to_dms(maxla, is_lat=False),
                        "step_deg": stepla,
                        "step_arcminutes": round(stepla * 60.0, 4),
                    },
                    "latitude": {
                        "min_deg": minphi,
                        "max_deg": maxphi,
                        "span_deg": round(maxphi - minphi, 8),
                        "min_dms": decdeg_to_dms(minphi, is_lat=True),
                        "max_dms": decdeg_to_dms(maxphi, is_lat=True),
                        "step_deg": stepphi,
                        "step_arcminutes": round(stepphi * 60.0, 4),
                    },
                }
            # Check for planimetric fields (mine, maxe, minn, maxn)
            elif "mine" in meta and "minn" in meta:
                mine = float(meta.get("mine", 0.0))
                maxe = float(meta.get("maxe", 0.0))
                minn = float(meta.get("minn", 0.0))
                maxn = float(meta.get("maxn", 0.0))
                stepe = float(meta.get("stepe", 0.0))
                stepn = float(meta.get("stepn", 0.0))

                layer_summary["formatted_extents"] = {
                    "easting": {
                        "min_m": mine,
                        "max_m": maxe,
                        "span_m": round(maxe - mine, 3),
                        "step_m": stepe,
                    },
                    "northing": {
                        "min_m": minn,
                        "max_m": maxn,
                        "span_m": round(maxn - minn, 3),
                        "step_m": stepn,
                    },
                }

    return layer_summary


def extract_full_spg_metadata(spg_path: Path) -> Dict[str, Any]:
    """Exhaustively extracts all metadata, parameters, and grid layer properties from an SPG file."""
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

    file_size_bytes = spg_path.stat().st_size
    sha256_hash = compute_sha256(spg_path)

    # 1. Base File Information
    result: Dict[str, Any] = {
        "_generator": "export_spg_metadata_json.py",
        "_generated_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "file_info": {
            "name": spg_path.name,
            "size_bytes": file_size_bytes,
            "size_formatted": format_bytes(file_size_bytes),
            "sha256": sha256_hash,
        },
        "metadata": {},
        "params": {},
        "grids": {},
    }

    # 2. Extract 100% of metadata attributes dynamically
    raw_metadata = data.get("metadata", {})
    if isinstance(raw_metadata, dict):
        result["metadata"] = normalize_value(raw_metadata)

    # 3. Extract 100% of params attributes dynamically
    raw_params = data.get("params", {})
    if isinstance(raw_params, dict):
        result["params"] = normalize_value(raw_params)

        # Enhance interpolation method descriptions if interpolation dict exists
        interp = raw_params.get("interpolation")
        if isinstance(interp, dict):
            method_names = {
                0: "INTERP_COLOCATE (Nearest-Neighbor)",
                1: "INTERP_BILINEAR (Bilinear)",
                2: "INTERP_BICUBIC (Bicubic Spline)",
            }
            h_val = interp.get("horizontal")
            v_val = interp.get("vertical")
            result["params"]["interpolation_descriptions"] = {
                "horizontal": method_names.get(h_val, f"Method {h_val}"),
                "vertical": method_names.get(v_val, f"Method {v_val}"),
            }

    # 4. Extract and analyze each layer in grids dynamically
    raw_grids = data.get("grids", {})
    if isinstance(raw_grids, dict):
        for layer_name, layer_dict in raw_grids.items():
            if isinstance(layer_dict, dict):
                result["grids"][layer_name] = analyze_layer(layer_name, layer_dict)

    # 5. Capture any other top-level keys that might appear in future releases
    for root_key, root_val in data.items():
        if root_key not in ["metadata", "params", "grids"]:
            result[root_key] = normalize_value(root_val)

    return result


def export_spg_metadata(spg_path: Path, output_json_path: Path) -> Path:
    """Extracts metadata and writes to a structured, pretty JSON file."""
    metadata = extract_full_spg_metadata(spg_path)
    output_json_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_json_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)
    print(f"Full SPG metadata exported successfully to: {output_json_path}")
    return output_json_path


def main():
    grids_dir = PROJECT_ROOT / "pytransdatro" / "grids"
    spg_candidates = list(grids_dir.glob("*.spg"))

    if len(sys.argv) > 1:
        spg_file = Path(sys.argv[1])
    elif len(spg_candidates) == 1:
        spg_file = spg_candidates[0]
    else:
        spg_file = grids_dir / "rom_grid3d_25.09.spg"

    default_output = grids_dir / "spg_metadata.json"
    output_file = Path(sys.argv[2]) if len(sys.argv) > 2 else default_output

    export_spg_metadata(spg_file, output_file)


if __name__ == "__main__":
    main()
