# PyTransDatRO: Telemetry Visualization Dashboard Specification

> [!NOTE]
> **Specification Status: IN PROGRESS (Draft)**
> This document details **Component 1: Usage Timeline Visual** and **Component 2: Spatial Overview Matrix Visual**. It will continue to expand as subsequent visual components (Modality Matrix, Overall Analytical Dashboard) are designed.

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
   *(Note: Planned logging schema addition)*.

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

## 3. Component Specification: Spatial Overview Matrix Visual

The Spatial Overview Matrix is the primary geographical visual of the dashboard. It conveys a clear macroscopic picture of where transformations occur across Romania without claiming false sub-kilometer precision or requiring map-navigation overhead.

### 3.1 Core Principles: GitHub Punchcard Style Macroscopic Matrix
- **Standalone Matrix (No Map Required)**: The visual is strictly an informative, punchcard-style overview inspired by GitHub activity grids. It does not require basemaps, administrative boundaries, or vector overlays; as usage accumulates over time, the silhouette and regional activity clusters naturally emerge from the populated cells.
- **Reduced Resolution ($0.1^\circ$ Decimal Degrees)**:
  - Raw database records ($0.01^\circ$ / ~1 km) are aggregated via integer division (`FLOOR(lat_grad / 10)`, `FLOOR(lon_grad / 10)`) to **1 decimal place** ($\approx 10\text{ km} \times 7.5\text{ km}$ bins).
  - High enough resolution to distinguish regional clusters and corridors; broad enough to remain fast, compact, and macro-oriented.
  - Romania's bounding box spans roughly $43.6^\circ - 48.3^\circ\text{ N}$ by $20.2^\circ - 29.7^\circ\text{ E}$.
- **Matrix Cell Appearance (Uniform Square Tiles)**:
  - Cells are rendered as uniform **square tiles** across the matrix layout.
  - **Faint Ghost Cells**: Empty/unvisited cells are explicitly rendered as subtle ghost squares (`rgba(..., 0.04)` / faint outline). This ensures clear structure even during early rollout stages or narrow date filters where only few cells are populated, preventing visual ambiguity between unvisited areas and rendering errors.
  - Populated cells take a color shade/intensity based on the aggregate hit count.
- **Independence from Component 1 (Timeline Visual)**:
  - The Spatial Overview Matrix is completely decoupled and independent from the Usage Timeline visual. It maintains its own independent temporal controls and filter state.

---

### 3.2 Coordinate Axes & Layout
- **X-Axis (Horizontal)**: Longitude in Decimal Degrees (span $20.0^\circ\text{ E}$ to $30.0^\circ\text{ E}$ with labels at $1.0^\circ$ or $0.5^\circ$ intervals).
- **Y-Axis (Vertical)**: Latitude in Decimal Degrees (span $43.5^\circ\text{ N}$ to $48.5^\circ\text{ N}$ with labels at $0.5^\circ$ intervals).
- **Cell Geometry**: Uniform square grid cells.

---

### 3.3 Metric Toggle & Scaling Modes

1. **Metric Selection Toggle**:
   - A toggle control allowing the user to select what metric the matrix displays:
     - **`[ Coordinates Transformed (Points) ]`**: Shading reflects cumulated count of transformed coordinates (`count_points`).
     - **`[ Service Invocations (Calls) ]`**: Shading reflects cumulated count of transformation calls (`count_calls`).
   - Switching the toggle instantly updates cell values and recalibrates color intensity.

2. **Scale Mode Toggle (Logarithmic vs. Linear)**:
   - **Logarithmic (Default)**: Normalizes extreme variations between high-density metropolitan/infrastructure hubs and rural cadastre cells.
   - **Linear**: Direct proportional scaling, useful for uniform inspection or balanced activity periods.
   - Both modes map values onto a clean 5-tier discrete color intensity scale (Level 0 = ghost empty, Levels 1–4 = increasing saturation).

3. **Accompanying Color Legend**: A clean 5-step discrete legend adjacent to the matrix indicating value intervals corresponding to each intensity tier.

---

### 3.4 Mathematical Quadrant Macro-Aggregator (Header Cards)

Directly above or alongside the matrix, four summary cards display aggregated metrics partitioned strictly by **mathematical quadrants**:
- **Reference Center**: Defined by the midpoint coordinates of the national bounding envelope (approx. $45.9^\circ\text{ N}, 24.9^\circ\text{ E}$).
- **Quadrants**:
  - **NW**: $\text{Lat} \ge 45.9^\circ\text{ N}, \text{Lon} < 24.9^\circ\text{ E}$
  - **NE**: $\text{Lat} \ge 45.9^\circ\text{ N}, \text{Lon} \ge 24.9^\circ\text{ E}$
  - **SW**: $\text{Lat} < 45.9^\circ\text{ N}, \text{Lon} < 24.9^\circ\text{ E}$
  - **SE**: $\text{Lat} < 45.9^\circ\text{ N}, \text{Lon} \ge 24.9^\circ\text{ E}$
- **Card Content**: Displays total volume for the active metric (Points or Calls) and its relative percentage of national volume (e.g. `NW: 342,100 pts (28%)`).
- **Rationale**: Strict mathematical notation avoids confusion or false correlation with administrative/county boundaries while providing immediate regional balance insights.

---

### 3.5 Dynamic Source Toggles & Filters

1. **Dynamic Source Toggle Buttons**:
   - Available client sources are retrieved dynamically from the backend/database (e.g., `web_map`, `rest_api`, `desktop_cli`).
   - Rendered as interactive toggle buttons/chips.
   - Users can toggle sources on/off; matrix and quadrant cards recalculate dynamically.

2. **Independent Temporal Controls**:
   - Own independent time horizon selector (Presets: `All-Time`, `Past 12 Months`, `Year to Date`, or custom date bounds).
   - Aggregates records whose `period_start` and `period_end` intersect the active filter.

3. **Operation & Dimensionality Filters**:
   - Filter by **Operation**: All, `Stereo70ToETRS89`, or `ETRS89ToStereo70`.
   - Filter by **Dimension**: All, `2D (Planimetric)`, or `3D (Elevation)`.

---

### 3.6 Hover Inspection
- Hovering the mouse over any cell highlights its border and reveals a lightweight tooltip:
  - **Coordinates**: e.g., `45.6° N, 25.6° E` (Center of the $0.1^\circ$ cell).
  - **Cell Bounds**: e.g., `[45.6° – 45.7° N, 25.6° – 25.7° E]`.
  - **Selected Metric Value**: e.g., `4,820 points` (or `142 calls`).
  - **Quadrant**: e.g., `SW`.
  - **Breakdown by Active Sources**: e.g., `web_map: 320 | rest_api: 2,100 | desktop_cli: 2,400`.
  - **Dominant Mode**: e.g., `3D (78%) | Stereo70 -> ETRS89 (64%)`.

---

## 4. Backend REST API Contract Requirements (Draft)

### Endpoint: `GET /api/v1/telemetry/spatial-matrix`

#### Query Parameters:
- `from_date` (optional, `YYYY-MM-DD`): Start date boundary.
- `to_date` (optional, `YYYY-MM-DD`): End date boundary.
- `metric` (`points` | `calls`, default `points`): Metric to aggregate per cell.
- `sources` (optional, comma-separated integer list): Active client source IDs (e.g., `1,2,3`).
- `operation` (optional, integer): `1` (st70_to_etrs89) or `2` (etrs89_to_st70).
- `is_3d` (optional, integer): `0` or `1`.

#### Response Schema:
```json
{
  "resolution": 0.1,
  "bounds": {
    "min_lat": 43.6,
    "max_lat": 48.3,
    "min_lon": 20.2,
    "max_lon": 29.7
  },
  "summary": {
    "active_cells_count": 1845,
    "total_metric_value": 49200000,
    "max_cell_hits": 352000,
    "quadrants": {
      "center": {"lat": 45.9, "lon": 24.9},
      "NW": {"value": 13776000, "percentage": 28.0},
      "NE": {"value": 9840000, "percentage": 20.0},
      "SW": {"value": 8856000, "percentage": 18.0},
      "SE": {"value": 16728000, "percentage": 34.0}
    }
  },
  "cells": [
    {
      "lat": 44.4,
      "lon": 26.1,
      "hits": 352000
    },
    {
      "lat": 46.7,
      "lon": 23.6,
      "hits": 128500
    }
  ]
}
```

---

## 5. Upcoming Visual Components (To Be Defined)

The following components will be specified in upcoming iterations:

- **Component 3: Transformation Modality & Operations Matrix**
  - Forward vs. Reverse distribution.
  - 2D vs. 3D height usage distribution.
  - Cross-tabulation by client platform and service type.
- **Component 4: Overall Analysis & Trends Dashboard**
  - Workday vs. weekend distribution ratios.
  - Romanian seasonal cadastral variation indices.
  - Multi-year adoption and growth curves.
