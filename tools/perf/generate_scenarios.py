"""Generates deterministic benchmark coordinate datasets for scenarios S1 to S4
and creates an interactive Leaflet HTML map visualizing all scenario layers.

Scenarios:
  - S1: Micro Clustered (10 points) - Cadastre / single construction site (single grid cell).
  - S2: Small Clustered (250 points) - City / municipality cadastre (1-2 grid cells).
  - S3: Medium Corridor (2,500 points) - Infrastructure corridor crossing Romania (dozens of cells).
  - S4: Large Dispersed (25,000 points) - Dispersed across Romanian territory (countrywide grid stress test).
"""

import os
import sys
import csv
import math
import json
import random
from pathlib import Path

# Ensure pytransdatro is importable from repo root
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))

import pytransdatro

def get_trans_validator():
    """Initializes and returns TransRO instance."""
    return pytransdatro.TransRO()

def generate_s1_micro(t, count=10, seed=42):
    """S1: 10 points within a single grid cell (Brașov / Prejmer area, N=455,000, E=565,000)."""
    random.seed(seed)
    center_n, center_e = 455000.0, 565000.0
    radius = 350.0  # within 350 meters
    pts = []
    pid = 1
    while len(pts) < count:
        angle = random.uniform(0, 2 * math.pi)
        r = random.uniform(0, radius)
        n = center_n + r * math.cos(angle)
        e = center_e + r * math.sin(angle)
        z = round(random.uniform(500.0, 550.0), 3)
        n = round(n, 3)
        e = round(e, 3)
        try:
            t.st70_to_etrs89(n, e, z)
            pts.append((f"S1_{pid:03d}", n, e, z))
            pid += 1
        except Exception:
            continue
    return pts

def generate_s2_small(t, count=250, seed=43):
    """S2: 250 points in a localized urban cluster (Bucharest metropolitan zone, N=335,000, E=585,000)."""
    random.seed(seed)
    center_n, center_e = 335000.0, 585000.0
    radius = 6000.0  # 6 km radius
    pts = []
    pid = 1
    while len(pts) < count:
        # Gaussian distribution centered on city center
        n = random.gauss(center_n, radius / 2.5)
        e = random.gauss(center_e, radius / 2.5)
        z = round(random.uniform(70.0, 110.0), 3)
        n = round(n, 3)
        e = round(e, 3)
        try:
            t.st70_to_etrs89(n, e, z)
            pts.append((f"S2_{pid:04d}", n, e, z))
            pid += 1
        except Exception:
            continue
    return pts

def generate_s3_medium(t, count=2500, seed=44):
    """S3: 2,500 points along a major transport corridor (Timișoara -> Deva -> Sibiu -> Brașov -> Ploiești -> Bucharest)."""
    random.seed(seed)
    # Waypoints in Stereo70 (N, E)
    waypoints = [
        (475000.0, 205000.0),  # Near Timișoara
        (485000.0, 310000.0),  # Near Deva
        (460000.0, 405000.0),  # Near Sibiu
        (455000.0, 560000.0),  # Near Brașov
        (395000.0, 585000.0),  # Near Ploiești
        (335000.0, 585000.0)   # Near Bucharest
    ]
    
    # Precompute segment cumulative lengths
    segs = []
    total_len = 0.0
    for i in range(len(waypoints) - 1):
        p1 = waypoints[i]
        p2 = waypoints[i + 1]
        dist = math.hypot(p2[0] - p1[0], p2[1] - p1[1])
        segs.append((p1, p2, dist, total_len))
        total_len += dist

    pts = []
    pid = 1
    while len(pts) < count:
        target_dist = random.uniform(0, total_len)
        # Find segment
        for p1, p2, seg_dist, cum_dist in segs:
            if cum_dist <= target_dist <= cum_dist + seg_dist:
                frac = (target_dist - cum_dist) / seg_dist
                base_n = p1[0] + frac * (p2[0] - p1[0])
                base_e = p1[1] + frac * (p2[1] - p1[1])
                # Jitter perpendicular to corridor (+/- 1.5 km)
                jitter_n = random.gauss(0, 750.0)
                jitter_e = random.gauss(0, 750.0)
                n = round(base_n + jitter_n, 3)
                e = round(base_e + jitter_e, 3)
                z = round(random.uniform(60.0, 600.0), 3)
                try:
                    t.st70_to_etrs89(n, e, z)
                    pts.append((f"S3_{pid:05d}", n, e, z))
                    pid += 1
                except Exception:
                    pass
                break
    return pts

def generate_s4_large(t, count=25000, seed=45):
    """S4: 25,000 points uniformly dispersed across the valid national territory of Romania."""
    random.seed(seed)
    
    # Valid grid extent bounds (from SPG metadata with a safe 25 km inward buffer)
    # E_min: 109,783, E_max: 890,783
    # N_min: 213,634, N_max: 785,634
    e_min_safe = 150000.0
    e_max_safe = 820000.0
    n_min_safe = 250000.0
    n_max_safe = 750000.0

    pts = []
    pid = 1
    attempts = 0
    while len(pts) < count:
        attempts += 1
        n = round(random.uniform(n_min_safe, n_max_safe), 3)
        e = round(random.uniform(e_min_safe, e_max_safe), 3)
        z = round(random.uniform(10.0, 1500.0), 3)
        try:
            t.st70_to_etrs89(n, e, z)
            pts.append((f"S4_{pid:06d}", n, e, z))
            pid += 1
        except Exception:
            continue
            
    print(f"Generated {count} dispersed points (acceptance rate: {count/attempts*100:.1f}%)")
    return pts

def save_csv(filepath, points):
    """Saves points to CSV: Point_ID,Northing,Easting,Z."""
    with open(filepath, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(["Point_ID", "Northing", "Easting", "Z"])
        for pt in points:
            writer.writerow([pt[0], f"{pt[1]:.3f}", f"{pt[2]:.3f}", f"{pt[3]:.3f}"])

def build_leaflet_map(output_path, scenario_data):
    """Builds a rich, interactive Leaflet HTML map displaying layers S1-S4 using Canvas rendering."""
    # scenario_data: dict of scenario_name -> list of (pid, n, e, z, lat_deg, lon_deg, h)
    layers_json = {}
    for sc_name, pts in scenario_data.items():
        layers_json[sc_name] = [
            {
                "id": p[0],
                "n": p[1],
                "e": p[2],
                "z": p[3],
                "lat": round(p[4], 6),
                "lon": round(p[5], 6),
                "h": round(p[6], 3)
            }
            for p in pts
        ]

    layers_payload = json.dumps(layers_json)

    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>PyTransDatRO - Benchmark Scenarios Map (S1 - S4)</title>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css" />
    <script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
    <style>
        body {{
            margin: 0;
            padding: 0;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
            background: #111827;
            color: #f3f4f6;
        }}
        #map {{
            width: 100vw;
            height: 100vh;
        }}
        .info-panel {{
            position: absolute;
            top: 16px;
            right: 16px;
            z-index: 1000;
            background: rgba(17, 24, 39, 0.92);
            backdrop-filter: blur(8px);
            padding: 16px 20px;
            border-radius: 10px;
            border: 1px solid rgba(255, 255, 255, 0.12);
            box-shadow: 0 10px 25px -5px rgba(0, 0, 0, 0.5);
            max-width: 320px;
        }}
        .info-panel h2 {{
            margin: 0 0 8px 0;
            font-size: 16px;
            font-weight: 600;
            color: #60a5fa;
        }}
        .info-panel p {{
            margin: 0 0 12px 0;
            font-size: 12px;
            color: #9ca3af;
            line-height: 1.4;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            margin-bottom: 6px;
            font-size: 12px;
        }}
        .legend-color {{
            width: 14px;
            height: 14px;
            border-radius: 50%;
            margin-right: 8px;
            flex-shrink: 0;
        }}
        .leaflet-popup-content-wrapper {{
            background: #1f2937;
            color: #f9fafb;
            border-radius: 8px;
            box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.4);
            border: 1px solid rgba(255, 255, 255, 0.1);
        }}
        .leaflet-popup-content {{
            font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
            font-size: 12px;
            line-height: 1.5;
            margin: 12px 14px;
        }}
        .leaflet-popup-tip {{
            background: #1f2937;
        }}
    </style>
</head>
<body>
    <div id="map"></div>
    <div class="info-panel">
        <h2>Benchmark Scenarios</h2>
        <p>Testing datasets for PyTransDatRO performance evaluation. Use the layer control in top-right to toggle layers.</p>
        <div class="legend-item">
            <div class="legend-color" style="background: #ef4444; border: 1px solid #ffffff;"></div>
            <span><b>S1: Micro Clustered</b> (10 pts)</span>
        </div>
        <div class="legend-item">
            <div class="legend-color" style="background: #f59e0b; border: 1px solid #ffffff;"></div>
            <span><b>S2: Small Clustered</b> (250 pts)</span>
        </div>
        <div class="legend-item">
            <div class="legend-color" style="background: #3b82f6; border: 1px solid #ffffff;"></div>
            <span><b>S3: Medium Corridor</b> (2,500 pts)</span>
        </div>
        <div class="legend-item">
            <div class="legend-color" style="background: #10b981; border: 1px solid #ffffff;"></div>
            <span><b>S4: Large Dispersed</b> (25,000 pts)</span>
        </div>
    </div>

    <script>
        var map = L.map('map', {{
            center: [45.9432, 24.9668],
            zoom: 7,
            preferCanvas: true
        }});

        // Base maps
        var cartoDark = L.tileLayer('https://{{s}}.basemaps.cartocdn.com/dark_all/{{z}}/{{x}}/{{y}}{{r}}.png', {{
            attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> &copy; <a href="https://carto.com/attributions">CARTO</a>',
            subdomains: 'abcd',
            maxZoom: 19
        }}).addTo(map);

        var topoLayer = L.tileLayer('https://{{s}}.tile.opentopomap.org/{{z}}/{{x}}/{{y}}.png', {{
            maxZoom: 17,
            attribution: '&copy; OpenStreetMap, SRTM | OpenTopoMap'
        }});

        var rawData = {layers_payload};
        var canvasRenderer = L.canvas({{ padding: 0.5 }});

        var styleConfigs = {{
            "S1": {{ color: "#ef4444", fillColor: "#f87171", radius: 7, weight: 2, fillOpacity: 0.95 }},
            "S2": {{ color: "#f59e0b", fillColor: "#fbbf24", radius: 5, weight: 1.5, fillOpacity: 0.85 }},
            "S3": {{ color: "#3b82f6", fillColor: "#60a5fa", radius: 3.5, weight: 1, fillOpacity: 0.75 }},
            "S4": {{ color: "#10b981", fillColor: "#34d399", radius: 2.0, weight: 0.5, fillOpacity: 0.60 }}
        }};

        var overlayGroups = {{}};

        for (var scKey in rawData) {{
            var cfg = styleConfigs[scKey];
            var pts = rawData[scKey];
            var grp = L.layerGroup();

            for (var i = 0; i < pts.length; i++) {{
                var p = pts[i];
                var marker = L.circleMarker([p.lat, p.lon], {{
                    renderer: canvasRenderer,
                    radius: cfg.radius,
                    color: cfg.color,
                    fillColor: cfg.fillColor,
                    fillOpacity: cfg.fillOpacity,
                    weight: cfg.weight
                }});

                var popupHtml = "<b>" + p.id + " (" + scKey + ")</b><br>" +
                                "<b>Stereo70:</b><br>" +
                                "N: " + p.n.toFixed(3) + " m<br>" +
                                "E: " + p.e.toFixed(3) + " m<br>" +
                                "Z: " + p.z.toFixed(3) + " m<br><br>" +
                                "<b>ETRS89:</b><br>" +
                                "Lat: " + p.lat.toFixed(6) + "°<br>" +
                                "Lon: " + p.lon.toFixed(6) + "°<br>" +
                                "h: " + p.h.toFixed(3) + " m";

                marker.bindPopup(popupHtml);
                grp.addLayer(marker);
            }}

            overlayGroups[scKey] = grp;
            // Add all layers to map by default
            grp.addTo(map);
        }}

        var baseMaps = {{
            "Carto Dark": cartoDark,
            "OpenTopoMap": topoLayer
        }};

        var layerControls = {{
            "<span style='color:#ef4444;'>●</span> S1: Micro Clustered (10)": overlayGroups["S1"],
            "<span style='color:#f59e0b;'>●</span> S2: Small Clustered (250)": overlayGroups["S2"],
            "<span style='color:#3b82f6;'>●</span> S3: Medium Corridor (2,500)": overlayGroups["S3"],
            "<span style='color:#10b981;'>●</span> S4: Large Dispersed (25,000)": overlayGroups["S4"]
        }};

        L.control.layers(baseMaps, layerControls, {{ collapsed: false, position: 'topleft' }}).addTo(map);
    </script>
</body>
</html>
"""
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)

def main():
    perf_dir = REPO_ROOT / "tools" / "perf"
    data_dir = perf_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    print("Initializing PyTransDatRO validator...")
    t = get_trans_validator()

    print("\n1. Generating Scenario S1 (Micro Clustered, 10 pts)...")
    s1_pts = generate_s1_micro(t, count=10)
    save_csv(data_dir / "s1_micro_clustered.csv", s1_pts)

    print("\n2. Generating Scenario S2 (Small Clustered, 250 pts)...")
    s2_pts = generate_s2_small(t, count=250)
    save_csv(data_dir / "s2_small_clustered.csv", s2_pts)

    print("\n3. Generating Scenario S3 (Medium Corridor, 2,500 pts)...")
    s3_pts = generate_s3_medium(t, count=2500)
    save_csv(data_dir / "s3_medium_corridor.csv", s3_pts)

    print("\n4. Generating Scenario S4 (Large Dispersed, 25,000 pts)...")
    s4_pts = generate_s4_large(t, count=25000)
    save_csv(data_dir / "s4_large_dispersed.csv", s4_pts)

    print("\nValidating all points and computing ETRS89 coordinates for map visualization...")
    scenario_map_data = {}
    for sc_name, pts in [("S1", s1_pts), ("S2", s2_pts), ("S3", s3_pts), ("S4", s4_pts)]:
        transformed = []
        for pid, n, e, z in pts:
            lat_rad, lon_rad, h = t.st70_to_etrs89(n, e, z)
            lat_deg = pytransdatro.utils.rad_to_deg(lat_rad)
            lon_deg = pytransdatro.utils.rad_to_deg(lon_rad)
            transformed.append((pid, n, e, z, lat_deg, lon_deg, h))
        scenario_map_data[sc_name] = transformed
        print(f"  {sc_name}: {len(transformed)} points validated 100% successfully.")

    map_file = perf_dir / "scenarios_map.html"
    print(f"\nGenerating interactive Leaflet map at: {map_file}...")
    build_leaflet_map(map_file, scenario_map_data)
    print("Map generated successfully!")
    print("\nAll datasets generated and stored in tools/perf/data/.")

if __name__ == "__main__":
    main()
