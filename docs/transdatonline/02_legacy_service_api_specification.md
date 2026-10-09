# TransDatOnline Legacy Service API Specification

## 1. Overview & Service Contract

This document provides the definitive specification of the legacy **TransDatOnline Coordinate Operation Web Service**. This service was released over 10 years ago to allow automated coordinate transformations between Romania's national projected system (**Stereo 70**) and European geographic coordinates (**ETRS89**).

To guarantee that existing clients (desktop scripts, GIS extensions, CAD plugins, web mashups) continue to function without modifications, the new REST API must maintain **100% backward compatibility** with this endpoint specification.

---

## 2. Endpoint Addressing

### Historical & Public Service URLs:
* **Historical Endpoint (University of Bucharest)**:
  `http://earth.unibuc.ro:8080/transdatonline/cooOpService`
* **Current Active Public Endpoint (geo-spatial.org)**:
  `https://www.geo-spatial.org/transdatonline/cooOpService`
* **New Service Compatibility Mounts**:
  * `/transdatonline/cooOpService` (Full path alias matching existing deployments)
  * `/cooOpService` (Convenience root alias)

---

## 3. Coordinate Conventions & Units

| Property | Rectangular Coordinates (Stereo 70) | Geographic Coordinates (ETRS89) |
| :--- | :--- | :--- |
| **Components** | Northing ($N$), Easting ($E$), Normal Height ($H$) | Latitude ($B$ / Lat), Longitude ($L$ / Lon), Ellipsoidal Height ($h$) |
| **Units** | **Meters** ($m$) | **Radians** ($rad$) for $B, L$; **Meters** ($m$) for $h$ |
| **Component Order** | `[N, E]` (2D) or `[N, E, H]` (3D) | `[Lat, Lon]` (2D) or `[Lat, Lon, h]` (3D) |

> [!IMPORTANT]
> The legacy service operates **strictly in radians** for geographic angles and **meters** for planar and height coordinates. It does **not** take or return degrees or DMS on this legacy endpoint.

---

## 4. HTTP GET Specification (Single Coordinate Transformation)

### Request
* **HTTP Method**: `GET`
* **Path**: `/transdatonline/cooOpService`
* **Query Parameters**:

| Parameter | Type | Required | Description | Example Values |
| :--- | :--- | :--- | :--- | :--- |
| `cooOp` | `string` | **Yes** | Coordinate operation name (case-sensitive) | `Stereo70ToETRS89`<br/>`ETRS89ToStereo70`<br/>`Stereo30ToETRS89`<br/>`ETRS89ToStereo30` |
| `coos` | `string` | **Yes** | Semicolon-separated (`";"`) coordinate numbers | `500000;500000;100` (3D)<br/>`500000;500000` (2D) |

#### Example GET Requests:
```http
GET /transdatonline/cooOpService?cooOp=Stereo70ToETRS89&coos=500000;500000;100 HTTP/1.1
Host: www.geo-spatial.org
```
```http
GET /transdatonline/cooOpService?cooOp=ETRS89ToStereo70&coos=0.8028465450500996;0.43630521911977493;139.7825537763764 HTTP/1.1
Host: www.geo-spatial.org
```

---

### Response Schemas (GET)

* **Content-Type**: `application/json` (historical servlet sometimes served `text/html;charset=UTF-8`, but payload is strict JSON)
* **Status Code**: `200 OK`

#### Scenario 1: Successful 3D Transformation
The output JSON object contains a single `coos` array with the transformed coordinates. The `warning` field is **omitted**.
```json
{
  "coos": [
    0.8028465450500997,
    0.43630521911977527,
    139.60844039916992
  ]
}
```

#### Scenario 2: Successful 2D Transformation
When 2 coordinates are provided in the request (`500000;500000`), the response contains a 2-element array:
```json
{
  "coos": [
    0.8028465450500997,
    0.43630521911977527
  ]
}
```

#### Scenario 3: Point Outside the Official Grid (`OutOfGridErr`)
When input coordinates fall outside the bounding envelope of the Romanian grid:
* The original input coordinates are **echoed back** in the `coos` array.
* The `warning` attribute is set to `"Out of grid"`.
```json
{
  "coos": [
    5.0,
    1.0,
    1.0
  ],
  "warning": "Out of grid"
}
```

#### Scenario 4: Point on No-Data Grid Cell (`NoDataGridErr`)
When coordinates are within grid extent, but required subgrid interpolation nodes have missing values:
```json
{
  "coos": [
    500000.0,
    500000.0,
    100.0
  ],
  "warning": "No data on grid"
}
```

#### Scenario 5: Non-numeric / Malformed Input
When input parameters cannot be parsed as floating-point numbers (e.g. `coos=a;b;c`):
* The `coos` array is **empty** (`[]`).
* The `warning` attribute is set to `"Invalid coordinate data"`.
```json
{
  "coos": [],
  "warning": "Invalid coordinate data"
}
```

---

## 5. HTTP POST Specification (Batch Coordinate Transformation)

### Request
* **HTTP Method**: `POST`
* **Path**: `/transdatonline/cooOpService`
* **Content-Type**: `application/x-www-form-urlencoded` or `multipart/form-data`
* **Parameters**:

| Parameter | Type | Required | Description | Example Values |
| :--- | :--- | :--- | :--- | :--- |
| `cooOp` | `string` | **Yes** | Coordinate operation name | `ETRS89ToStereo70`<br/>`Stereo70ToETRS89` |
| `coosArray` | `string` | **Yes** | Stringified JSON array of coordinate objects | `[{"coos":[0.802846, 0.436305, 89.439]}, {"coos":[0.902846, 0.536305, 90.439]}]` |

#### Example POST Body:
```http
POST /transdatonline/cooOpService HTTP/1.1
Host: www.geo-spatial.org
Content-Type: application/x-www-form-urlencoded

cooOp=ETRS89ToStereo70&coosArray=%5B%7B%22coos%22%3A%5B0.8028465450500996%2C0.43630521911977493%2C89.43911984920247%5D%7D%2C%7B%22coos%22%3A%5B0.9028465450500996%2C0.53630521911977493%2C90.43911984920247%5D%7D%5D
```

---

### Response Schema (POST)

The response for `POST` is a **JSON array of coordinate objects**, preserving the exact order of the input array.

```json
[
  {
    "coos": [
      500000.0001279776,
      499999.99995336385,
      49.83067945003255
    ]
  },
  {
    "coos": [
      0.9028465450500996,
      0.53630521911977493,
      90.43911984920247
    ],
    "warning": "Out of grid"
  }
]
```

---

## 6. Numerical Verification: Legacy vs Modern Model

A side-by-side comparison was conducted for point $(N = 500000.0\text{ m}, E = 500000.0\text{ m}, H = 100.0\text{ m})$:

| Field | Legacy Service (2012 EGG97 Model) | Modern `PyTransDatRO` (SPG Model) | Difference / Delta |
| :--- | :--- | :--- | :--- |
| **Latitude ($B$)** | `0.8028465450500996 rad` | `0.8028465450500997 rad` | **$1.1 \times 10^{-16}\text{ rad}$ ($\approx 0.0007\text{ mm}$)** |
| **Longitude ($L$)** | `0.43630521911977493 rad` | `0.43630521911977527 rad` | **$3.4 \times 10^{-16}\text{ rad}$ ($\approx 0.0017\text{ mm}$)** |
| **Height ($h$)** | `139.7825537763764 m` | `139.60844039916992 m` | **$-0.1741\text{ m}$ ($\approx 17.4\text{ cm}$)** |

> [!NOTE]
> The planimetric 2D coordinates are in near-perfect numerical agreement (sub-micrometer precision). The elevation difference of $\approx 17.4\text{ cm}$ is expected and caused by the transition from the legacy EGG97 quasigeoid model to the modern official Romanian quasigeoid grid in `.spg`.

---

## 7. Stereo 30 Handling Strategy

In the legacy Java implementation, operations `Stereo30ToETRS89` and `ETRS89ToStereo30` were supported via dedicated MySQL tables (`etrs89_stereo30_2d`).

In modern Romanian geodesy and `pytransdatro`, Stereo 30 is retired and not included in the modern `.spg` package.

### Recommended Compatibility Behavior:
* If a client calls `Stereo30ToETRS89` or `ETRS89ToStereo30`:
  * Return the input coordinates intact in `"coos"`.
  * Return `"warning": "Stereo30 transformation is obsolete and unsupported"`.
  * Maintain HTTP 200 and standard JSON structure so the client parser does not crash.
