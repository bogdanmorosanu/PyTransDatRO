import os
import sys
import csv
import math
import json
from pathlib import Path

# Add project root to sys path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytransdatro

def get_flag_m(diff):
    abs_diff = abs(diff)
    if abs_diff <= 0.001:
        return 0
    elif abs_diff <= 0.002:
        return 1
    elif abs_diff <= 0.003:
        return 2
    else:
        return 4

def get_flag_sec(diff):
    abs_diff = abs(diff)
    if abs_diff <= 0.00001:
        return 0
    elif abs_diff <= 0.00002:
        return 1
    elif abs_diff <= 0.00003:
        return 2
    else:
        return 4

def get_color_for_flag(flag):
    if flag == 0:
        return 'green'
    elif flag == 1:
        return 'yellow'
    elif flag == 2:
        return 'orange'
    else:
        return 'red'

def read_csv(filepath, has_sexa_latlon=False):
    data = {}
    with open(filepath, 'r', encoding='ANSI') as f:
        for line in f:
            if not line.strip():
                continue
            vals = [v.strip() for v in line.split(',')]
            pid = vals[0]
            if has_sexa_latlon:
                lat = pytransdatro.utils.sexa_dms_chars_to_rad(vals[1])
                lon = pytransdatro.utils.sexa_dms_chars_to_rad(vals[2])
                z = float(vals[3])
                data[pid] = (lat, lon, z, vals[1], vals[2]) # Also store raw string
            else:
                n = float(vals[1])
                e = float(vals[2])
                z = float(vals[3])
                data[pid] = (n, e, z)
    return data

def write_summary(filepath, direction, total, errs, flags_e, flags_n, flags_z):
    max_e, max_n, max_z = errs['max_e'], errs['max_n'], errs['max_z']
    sum_e, sum_n, sum_z = errs['sum_e'], errs['sum_n'], errs['sum_z']
    
    mean_e = sum_e / total if total > 0 else 0
    mean_n = sum_n / total if total > 0 else 0
    mean_z = sum_z / total if total > 0 else 0
    
    # Check if pass
    passed = (flags_e.get(4, 0) == 0) and (flags_n.get(4, 0) == 0) and (flags_z.get(4, 0) == 0)
    
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(f"Transformation Report Summary: {direction}\n")
        f.write("="*50 + "\n")
        f.write(f"Overall Status: {'PASS' if passed else 'FAIL'}\n\n")
        
        f.write("Maximum Absolute Errors:\n")
        if direction == "ETRS89 to Stereo70":
            f.write(f"  Easting:   {max_e:.6f} m\n")
            f.write(f"  Northing:  {max_n:.6f} m\n")
        else:
            f.write(f"  Longitude: {max_e:.8f} sec\n")
            f.write(f"  Latitude:  {max_n:.8f} sec\n")
        f.write(f"  Elevation: {max_z:.6f} m\n\n")
        
        f.write("Mean Absolute Errors:\n")
        if direction == "ETRS89 to Stereo70":
            f.write(f"  Easting:   {mean_e:.6f} m\n")
            f.write(f"  Northing:  {mean_n:.6f} m\n")
        else:
            f.write(f"  Longitude: {mean_e:.8f} sec\n")
            f.write(f"  Latitude:  {mean_n:.8f} sec\n")
        f.write(f"  Elevation: {mean_z:.6f} m\n\n")
        
        f.write("Precision Distribution (Count of points per flag):\n")
        labels = ["Easting/Longitude", "Northing/Latitude", "Elevation"]
        stats = [flags_e, flags_n, flags_z]
        for label, stat in zip(labels, stats):
            f.write(f"  {label}:\n")
            f.write(f"    Flag 0: {stat.get(0, 0)}\n")
            f.write(f"    Flag 1: {stat.get(1, 0)}\n")
            f.write(f"    Flag 2: {stat.get(2, 0)}\n")
            f.write(f"    Flag 4: {stat.get(4, 0)}\n")
            f.write("\n")

def write_html_map(filepath, points, direction):
    # points format: [{'lat': lat_deg, 'lon': lon_deg, 'color_2d': color, 'color_1d': color, 'popup': html_string}]
    points_json = json.dumps(points)
    
    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>Transformation Precision Map - {direction}</title>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <link rel="stylesheet" href="https://unpkg.com/leaflet/dist/leaflet.css" />
    <script src="https://unpkg.com/leaflet/dist/leaflet.js"></script>
    <style>
        body {{ margin: 0; padding: 0; }}
        #map {{ width: 100%; height: 100vh; }}
        .leaflet-popup-content {{ font-family: monospace; }}
    </style>
</head>
<body>
    <div id="map"></div>
    <script>
        var map = L.map('map').setView([45.9, 24.9], 7); // Center of Romania

        var baseLayer = L.tileLayer('https://{{s}}.tile.opentopomap.org/{{z}}/{{x}}/{{y}}.png', {{
            maxZoom: 17,
            attribution: 'Map data: &copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors, <a href="http://viewfinderpanoramas.org">SRTM</a> | Map style: &copy; <a href="https://opentopomap.org">OpenTopoMap</a> (<a href="https://creativecommons.org/licenses/by-sa/3.0/">CC-BY-SA</a>)'
        }}).addTo(map);

        var points = {points_json};
        
        // Custom marker colors using simple SVG icons
        function getIcon(color) {{
            return L.divIcon({{
                className: 'custom-div-icon',
                html: "<div style='background-color:" + color + "; width:12px; height:12px; border-radius:50%; border:1px solid black;'></div>",
                iconSize: [12, 12],
                iconAnchor: [6, 6]
            }});
        }}

        var layer2D = L.layerGroup().addTo(map);
        var layer1D = L.layerGroup();

        points.forEach(function(p) {{
            L.marker([p.lat, p.lon], {{icon: getIcon(p.color_2d)}})
             .bindPopup(p.popup)
             .addTo(layer2D);
             
            L.marker([p.lat, p.lon], {{icon: getIcon(p.color_1d)}})
             .bindPopup(p.popup)
             .addTo(layer1D);
        }});
        
        var overlayMaps = {{
            "2D Precision (Lat/Lon or E/N)": layer2D,
            "1D Precision (Elevation Z)": layer1D
        }};
        L.control.layers(null, overlayMaps, {{collapsed: false}}).addTo(map);
    </script>
</body>
</html>
"""
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(html)

def generate_reports():
    root_dir = Path(__file__).parent.parent
    data_dir = root_dir / "tests" / "data"
    reports_dir = root_dir / "reports"
    reports_dir.mkdir(exist_ok=True)
    
    t = pytransdatro.TransRO()
    
    # 1. ETRS89 to Stereo70
    print("Generating ETRS89 to Stereo70 reports...")
    in_file = data_dir / "etrs89_to_st70_input.csv"
    exp_file = data_dir / "etrs89_to_st70_expected.csv"
    
    if in_file.exists() and exp_file.exists():
        inputs = read_csv(in_file, has_sexa_latlon=True)
        expecteds = read_csv(exp_file, has_sexa_latlon=False)
        
        errs = {'max_e': 0, 'max_n': 0, 'max_z': 0, 'sum_e': 0, 'sum_n': 0, 'sum_z': 0}
        flags_e, flags_n, flags_z = {}, {}, {}
        map_points = []
        
        with open(reports_dir / "etrs89_to_st70_report.csv", "w", newline="", encoding='utf-8') as f_csv:
            writer = csv.writer(f_csv)
            writer.writerow(["Point_ID", "Input_Lat_sexa", "Input_Lon_sexa", "Input_Z", 
                             "Expected_E", "Expected_N", "Expected_Z", 
                             "Computed_E", "Computed_N", "Computed_Z", 
                             "Diff_E", "Diff_N", "Diff_Z", 
                             "Flag_E", "Flag_N", "Flag_Z"])
            
            for pid, (lat_rad, lon_rad, z_in, lat_sexa, lon_sexa) in inputs.items():
                if pid not in expecteds:
                    continue
                exp_n, exp_e, exp_z = expecteds[pid]
                
                try:
                    comp_n, comp_e, comp_z = t.etrs89_to_st70(lat_rad, lon_rad, z_in)
                except Exception as e:
                    print(f"Error computing point {pid}: {e}")
                    continue
                
                diff_e = comp_e - exp_e
                diff_n = comp_n - exp_n
                diff_z = comp_z - exp_z
                
                errs['max_e'] = max(errs['max_e'], abs(diff_e))
                errs['max_n'] = max(errs['max_n'], abs(diff_n))
                errs['max_z'] = max(errs['max_z'], abs(diff_z))
                errs['sum_e'] += abs(diff_e)
                errs['sum_n'] += abs(diff_n)
                errs['sum_z'] += abs(diff_z)
                
                fl_e = get_flag_m(diff_e)
                fl_n = get_flag_m(diff_n)
                fl_z = get_flag_m(diff_z)
                
                flags_e[fl_e] = flags_e.get(fl_e, 0) + 1
                flags_n[fl_n] = flags_n.get(fl_n, 0) + 1
                flags_z[fl_z] = flags_z.get(fl_z, 0) + 1
                
                writer.writerow([pid, lat_sexa, lon_sexa, f"{z_in:.3f}",
                                 f"{exp_e:.3f}", f"{exp_n:.3f}", f"{exp_z:.3f}",
                                 f"{comp_e:.3f}", f"{comp_n:.3f}", f"{comp_z:.3f}",
                                 f"{diff_e:.6f}", f"{diff_n:.6f}", f"{diff_z:.6f}",
                                 fl_e, fl_n, fl_z])
                
                # For map
                flag_2d = max(fl_e, fl_n)
                flag_1d = fl_z
                popup = (f"<b>Point ID:</b> {pid}<br>"
                         f"<b>Computed:</b> E {comp_e:.3f}, N {comp_n:.3f}, Z {comp_z:.3f}<br>"
                         f"<b>Expected:</b> E {exp_e:.3f}, N {exp_n:.3f}, Z {exp_z:.3f}<br>"
                         f"<b>Diffs:</b> E {diff_e:.3f}, N {diff_n:.3f}, Z {diff_z:.3f}")
                map_points.append({
                    'lat': pytransdatro.utils.rad_to_deg(lat_rad),
                    'lon': pytransdatro.utils.rad_to_deg(lon_rad),
                    'color_2d': get_color_for_flag(flag_2d),
                    'color_1d': get_color_for_flag(flag_1d),
                    'popup': popup
                })
        
        write_summary(reports_dir / "etrs89_to_st70_summary.txt", "ETRS89 to Stereo70", len(map_points), errs, flags_e, flags_n, flags_z)
        write_html_map(reports_dir / "etrs89_to_st70_map.html", map_points, "ETRS89 to Stereo70")
        print("Done ETRS89 to Stereo70.")
    
    # 2. Stereo70 to ETRS89
    print("Generating Stereo70 to ETRS89 reports...")
    in_file = data_dir / "st70_to_etrs89_input.csv"
    exp_file = data_dir / "st70_to_etrs89_expected.csv"
    
    if in_file.exists() and exp_file.exists():
        inputs = read_csv(in_file, has_sexa_latlon=False)
        expecteds = read_csv(exp_file, has_sexa_latlon=True)
        
        errs = {'max_e': 0, 'max_n': 0, 'max_z': 0, 'sum_e': 0, 'sum_n': 0, 'sum_z': 0}
        flags_e, flags_n, flags_z = {}, {}, {}
        map_points = []
        
        with open(reports_dir / "st70_to_etrs89_report.csv", "w", newline="", encoding='utf-8') as f_csv:
            writer = csv.writer(f_csv)
            writer.writerow(["Point_ID", "Input_E", "Input_N", "Input_Z", 
                             "Expected_Lat_sexa", "Expected_Lon_sexa", "Expected_Z", 
                             "Computed_Lat_sexa", "Computed_Lon_sexa", 
                             "Computed_Lat_deg", "Computed_Lon_deg", "Computed_Z", 
                             "Diff_Lat_sec", "Diff_Lon_sec", "Diff_Z", 
                             "Flag_Lat", "Flag_Lon", "Flag_Z"])
            
            for pid, (n_in, e_in, z_in) in inputs.items():
                if pid not in expecteds:
                    continue
                exp_lat_rad, exp_lon_rad, exp_z, exp_lat_sexa, exp_lon_sexa = expecteds[pid]
                
                try:
                    comp_lat_rad, comp_lon_rad, comp_z = t.st70_to_etrs89(n_in, e_in, z_in)
                except Exception as e:
                    print(f"Error computing point {pid}: {e}")
                    continue
                
                # compute difference in radians, then convert to arc seconds
                # 1 rad = 360/(2pi) degrees = 360 * 3600 / (2pi) arc seconds
                rad_to_sec = 3600 * 180 / math.pi
                
                diff_lat_rad = comp_lat_rad - exp_lat_rad
                diff_lon_rad = comp_lon_rad - exp_lon_rad
                diff_z = comp_z - exp_z
                
                diff_lat_sec = diff_lat_rad * rad_to_sec
                diff_lon_sec = diff_lon_rad * rad_to_sec
                
                errs['max_e'] = max(errs['max_e'], abs(diff_lon_sec))
                errs['max_n'] = max(errs['max_n'], abs(diff_lat_sec))
                errs['max_z'] = max(errs['max_z'], abs(diff_z))
                errs['sum_e'] += abs(diff_lon_sec)
                errs['sum_n'] += abs(diff_lat_sec)
                errs['sum_z'] += abs(diff_z)
                
                fl_e = get_flag_sec(diff_lon_sec)
                fl_n = get_flag_sec(diff_lat_sec)
                fl_z = get_flag_m(diff_z)
                
                flags_e[fl_e] = flags_e.get(fl_e, 0) + 1
                flags_n[fl_n] = flags_n.get(fl_n, 0) + 1
                flags_z[fl_z] = flags_z.get(fl_z, 0) + 1
                
                comp_lat_deg = pytransdatro.utils.rad_to_deg(comp_lat_rad)
                comp_lon_deg = pytransdatro.utils.rad_to_deg(comp_lon_rad)
                
                # rad_to_sexa produces string like "47 42 56.40000"
                comp_lat_sexa = pytransdatro.utils.rad_to_sexa(comp_lat_rad)
                comp_lon_sexa = pytransdatro.utils.rad_to_sexa(comp_lon_rad)
                
                writer.writerow([pid, f"{e_in:.3f}", f"{n_in:.3f}", f"{z_in:.3f}",
                                 exp_lat_sexa, exp_lon_sexa, f"{exp_z:.3f}",
                                 comp_lat_sexa, comp_lon_sexa, 
                                 f"{comp_lat_deg:.8f}", f"{comp_lon_deg:.8f}", f"{comp_z:.3f}",
                                 f"{diff_lat_sec:.8f}", f"{diff_lon_sec:.8f}", f"{diff_z:.6f}",
                                 fl_n, fl_e, fl_z])
                
                # For map
                flag_2d = max(fl_e, fl_n)
                flag_1d = fl_z
                popup = (f"<b>Point ID:</b> {pid}<br>"
                         f"<b>Computed:</b> Lat {comp_lat_sexa}, Lon {comp_lon_sexa}, Z {comp_z:.3f}<br>"
                         f"<b>Expected:</b> Lat {exp_lat_sexa}, Lon {exp_lon_sexa}, Z {exp_z:.3f}<br>"
                         f"<b>Diffs:</b> Lat {diff_lat_sec:.6f}\", Lon {diff_lon_sec:.6f}\", Z {diff_z:.3f}m")
                map_points.append({
                    'lat': comp_lat_deg,
                    'lon': comp_lon_deg,
                    'color_2d': get_color_for_flag(flag_2d),
                    'color_1d': get_color_for_flag(flag_1d),
                    'popup': popup
                })
        
        write_summary(reports_dir / "st70_to_etrs89_summary.txt", "Stereo70 to ETRS89", len(map_points), errs, flags_e, flags_n, flags_z)
        write_html_map(reports_dir / "st70_to_etrs89_map.html", map_points, "Stereo70 to ETRS89")
        print("Done Stereo70 to ETRS89.")
        print(f"\nAll reports generated successfully in {reports_dir.absolute()}")

if __name__ == "__main__":
    generate_reports()
