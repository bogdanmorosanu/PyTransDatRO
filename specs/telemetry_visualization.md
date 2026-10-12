# PyTransDatRO: Telemetry Visualization Dashboard Specification

> [!NOTE]
> **Specification Status: IN PROGRESS (Draft)**
> This document currently details **Component 1: Usage Timeline Visual**. It is deliberately incomplete as further visual components (Spatial Heatmap, Modality Matrix, Overall Analytical Dashboard) are being designed and will be added in subsequent iterations.

---

## 1. Overview & Objectives

The goal of the Telemetry Visualization Dashboard is to provide a dedicated, interactive window/view within the PyTransDatRO ecosystem to visualize operational usage data collected in `logs/usage.db`.

### Key Design Principles:
- **Separation of Concerns**: The telemetry dashboard is a standalone data visualization window, completely decoupled from the coordinate conversion workbench and web map tools.
- **Data-Display First**: The primary objective is faithful, clear, and interactive presentation of raw telemetry dimensions (time-series, platforms, operations, spatial distribution) before performing secondary statistical analysis.
- **REST API Driven**: All visuals query the backend REST API, which reads the SQLite database (`mode=ro`) to serve aggregated telemetry metrics dynamically.
- **Ergonomic Interactivity**: Interactive time window scrubbing, custom date range selection, multi-attribute filtering, dynamic KPI recalculation, and hover data inspection.

---

## 2. Component Specification: Usage Timeline Visual

The Usage Timeline is the primary temporal visual of the dashboard, displaying the volume and intensity of service utilization over time.

### 2.1 Coordinate Axes & Scaling
- **X-Axis (Time)**:
  - Formatted chronologically (ISO dates `YYYY-MM-DD` or month labels `YYYY-MM`).
  - Supports continuous scrubbing and zooming across 3+ years of logged history.
- **Dual Y-Axis (Handling Scale Discrepancies)**:
  - **Left Y-Axis (Calls / Invocations)**:
    - Scales to request volume (typically hundreds to low thousands per day).
    - Clear integer formatting (e.g., `0`, `500`, `1,000`, `2,500`).
  - **Right Y-Axis (Coordinates / Points Transformed)**:
    - Scales to coordinate throughput (typically tens of thousands to hundreds of thousands per day, heavily driven by batch/file operations).
    - Clear formatted intervals (e.g., `0`, `50k`, `100k`, `250k`).
  - Both metrics share the exact same temporal X-axis and hover cursor, enabling instant visual correlation between request frequency and data volume.

---

### 2.2 Data Attributes & Grouping Categories

The chart must support grouping, stacking, and filtering by the following dimensions:

1. **Client Source (`source`)**:
   - `web_map` (Interactive web application)
   - `rest_api` (Backend services and automated pipelines)
   - `desktop_cli` (Desktop GIS tools and batch conversions)
   - `unknown` (Unspecified client)

2. **Service Type (`service_name`)**:
   - `point` (Single coordinate conversion)
   - `batch` (Array/list of coordinates)
   - `file` (Multipart file conversion: CSV, TXT, etc.)
   *(Note: Requires upcoming logging logic enhancement to record `service_name` in SQLite)*.

3. **Operation Parameter (`operation` / Transformation Direction)**:
   - `Stereo70ToETRS89` (Forward transformation)
   - `ETRS89ToStereo70` (Inverse transformation)

4. **Dimensionality (`is_3d`)**:
   - `2D` (Planimetric)
   - `3D` (Elevation / Quasigeoid height included)

---

### 2.3 Time Window Navigation & Filters

1. **Quick Preset Toggles**:
   - `[ All Time (3 Years) ]`
   - `[ Past 12 Months ]`
   - `[ Year to Date (YTD) ]`
   - `[ Past 30 Days ]`
   - `[ Past 7 Days ]`

2. **Custom Date Range Selector**:
   - Start Date (`YYYY-MM-DD`) and End Date (`YYYY-MM-DD`) input pickers.
   - Allows users to isolate arbitrary periods (e.g., specific project deadlines or seasonal campaigns).

3. **Aggregation Granularity**:
   - `Daily`: Individual day bars/points (ideal for ranges $\le 90$ days).
   - `Monthly`: Aggregated monthly totals (ideal for multi-year macro trends).
   - *Behavior*: Can manually toggle or automatically adapt based on active range width.

4. **Mini-Timeline Scrubber (Brush)**:
   - A miniature timeline ribbon positioned beneath the main chart.
   - Dual draggable handles allow smooth dragging, panning, and zooming across the historical timeline.

---

### 2.4 Dynamic KPI Scope Summary (Header Card)

Directly above the chart area, a responsive summary banner displays dynamic totals reflecting the **currently filtered view**:

- **Total Invocations**: Sum of visible `count_calls`.
- **Total Points Processed**: Sum of visible `count_points`.
- **Batch Intensity Ratio**: Average coordinates per invocation ($\text{Points} / \text{Calls}$).
- **Active Filter Period**: e.g., `2025-03-01 to 2025-06-30 (122 days)`.

Whenever the user changes the time window or toggles categories on/off, these figures update immediately.

---

### 2.5 Interactivity & Hover Inspection

1. **Interactive Legend & Category Filter Chips**:
   - Clickable pill buttons/checkboxes for active categories (`Source`, `Service`, `Operation`).
   - Clicking toggles series visibility on the chart.
   - Alt-click / double-click isolates a single category.

2. **Crosshair & Multi-Series Tooltip**:
   - Hovering over any point on the timeline displays a synchronized vertical crosshair guide.
   - A high-contrast floating tooltip appears containing:
     - **Date**: e.g., `Friday, 14 Nov 2025`.
     - **Calls Breakdown**: Calls per active category and total daily calls.
     - **Points Breakdown**: Coordinates per active category and total daily points.
     - **Points/Call Ratio** for that specific day.

---

## 3. Backend REST API Contract Requirements (Draft)

To support this visual, the backend API should expose a dedicated timeline query endpoint:

### Endpoint: `GET /api/v1/telemetry/timeline`

#### Query Parameters:
- `from_date` (optional, `YYYY-MM-DD`): Start date filter.
- `to_date` (optional, `YYYY-MM-DD`): End date filter.
- `granularity` (`daily` | `monthly`, default `daily`): Aggregation bucket.
- `source` (optional, list of integers): Filter by client source IDs (`1`, `2`, `3`).
- `service` (optional, string): Filter by service (`point`, `batch`, `file`).
- `operation` (optional, string): Filter by operation (`Stereo70ToETRS89`, `ETRS89ToStereo70`).

#### Response Schema:
```json
{
  "summary": {
    "total_calls": 142500,
    "total_points": 4820100,
    "avg_points_per_call": 33.8,
    "from_date": "2025-01-01",
    "to_date": "2025-06-30",
    "days_count": 181
  },
  "series": [
    {
      "date": "2025-01-01",
      "calls": {
        "web_map": 12,
        "rest_api": 45,
        "desktop_cli": 18,
        "total": 75
      },
      "points": {
        "web_map": 28,
        "rest_api": 840,
        "desktop_cli": 4200,
        "total": 5068
      }
    }
  ]
}
```

---

## 4. Upcoming Visual Components (To Be Defined)

The following components will be specified in upcoming revisions:

- **Component 2: Spatial Density & Romanian Geographic Heatmap**
  - Full-canvas Leaflet/Mapbox rendering of ~1 km binned cells (`lat_deg + lat_grad/100.0`, `lon_deg + lon_grad/100.0`).
  - Top hotspots ranking table (Bucharest, Cluj, Timișoara, Iași, Constanța, etc.).
- **Component 3: Transformation Modality & Operations Matrix**
  - Forward vs. Reverse distribution.
  - 2D vs. 3D height usage distribution.
  - Cross-tabulation by client platform.
- **Component 4: Overall Analysis & Trends Dashboard**
  - Workday vs. weekend distribution ratios.
  - Romanian seasonal cadastral variation indices.
  - Multi-year adoption and growth curves.
