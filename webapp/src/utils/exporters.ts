/**
 * Client-Side Multi-Format Export Generator (CSV, GeoJSON, AutoCAD DXF)
 */

export interface ExportPointItem {
  pointId: string | null;
  coords: number[]; // [North/Lat, East/Long, (Elevation)]
  warning?: string | null;
  isGeo: boolean;
}

function triggerDownload(blob: Blob, filename: string): void {
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

/**
 * Exports transformed points as a CSV file.
 */
export function exportToCSV(
  items: ExportPointItem[],
  delimiter: string = ',',
  order: 'NE' | 'EN' = 'NE',
  filename: string = 'transformed_coordinates.csv'
): void {
  const sep = delimiter === ' ' ? ' ' : delimiter;
  const lines: string[] = [];

  for (const item of items) {
    if (!item.coords || item.coords.length < 2) continue;
    const parts: string[] = [];
    if (item.pointId) parts.push(item.pointId);

    const c1 = item.coords[0]; // North or Lat
    const c2 = item.coords[1]; // East or Long
    const c3 = item.coords[2]; // Elevation

    if (order === 'EN') {
      parts.push(item.isGeo ? c2.toFixed(8) : c2.toFixed(3));
      parts.push(item.isGeo ? c1.toFixed(8) : c1.toFixed(3));
    } else {
      parts.push(item.isGeo ? c1.toFixed(8) : c1.toFixed(3));
      parts.push(item.isGeo ? c2.toFixed(8) : c2.toFixed(3));
    }

    if (c3 !== undefined) {
      parts.push(c3.toFixed(3));
    }

    if (item.warning) {
      parts.push(`# ${item.warning}`);
    }

    lines.push(parts.join(sep));
  }

  const blob = new Blob([lines.join('\r\n')], { type: 'text/csv;charset=utf-8;' });
  triggerDownload(blob, filename);
}

/**
 * Exports transformed points as a GeoJSON FeatureCollection.
 */
export function exportToGeoJSON(
  items: ExportPointItem[],
  filename: string = 'transformed_coordinates.geojson'
): void {
  const features = items
    .filter((it) => it.coords && it.coords.length >= 2 && !it.warning)
    .map((it, idx) => {
      // GeoJSON expects [Longitude/X, Latitude/Y, (Elevation/Z)]
      const x = it.isGeo ? it.coords[1] : it.coords[1]; // Easting or Longitude
      const y = it.isGeo ? it.coords[0] : it.coords[0]; // Northing or Latitude
      const z = it.coords[2] !== undefined ? it.coords[2] : null;

      const coordinates = z !== null ? [x, y, z] : [x, y];

      return {
        type: 'Feature',
        id: it.pointId || idx + 1,
        geometry: {
          type: 'Point',
          coordinates,
        },
        properties: {
          id: it.pointId || String(idx + 1),
          elevation: z,
          system: it.isGeo ? 'ETRS89' : 'Stereo70',
        },
      };
    });

  const geojson = {
    type: 'FeatureCollection',
    features,
  };

  const blob = new Blob([JSON.stringify(geojson, null, 2)], {
    type: 'application/geo+json;charset=utf-8;',
  });
  triggerDownload(blob, filename);
}

/**
 * Exports transformed points as an AutoCAD ASCII DXF file.
 */
export function exportToDXF(
  items: ExportPointItem[],
  filename: string = 'transformed_coordinates.dxf'
): void {
  const dxfLines: string[] = [
    '0', 'SECTION',
    '2', 'HEADER',
    '0', 'ENDSEC',
    '0', 'SECTION',
    '2', 'ENTITIES',
  ];

  for (const it of items) {
    if (!it.coords || it.coords.length < 2 || it.warning) continue;

    // CAD coordinate convention: X = Easting, Y = Northing, Z = Elevation
    const x = it.isGeo ? it.coords[1] : it.coords[1];
    const y = it.isGeo ? it.coords[0] : it.coords[0];
    const z = it.coords[2] !== undefined ? it.coords[2] : 0.0;

    // POINT entity
    dxfLines.push(
      '0', 'POINT',
      '8', 'POINTS',
      '10', x.toFixed(4),
      '20', y.toFixed(4),
      '30', z.toFixed(4)
    );

    // If pointId is present, write a TEXT label near the point
    if (it.pointId) {
      dxfLines.push(
        '0', 'TEXT',
        '8', 'LABELS',
        '10', (x + 0.5).toFixed(4),
        '20', (y + 0.5).toFixed(4),
        '30', z.toFixed(4),
        '40', '1.0', // Text height
        '1', it.pointId
      );
    }
  }

  dxfLines.push('0', 'ENDSEC', '0', 'EOF');

  const blob = new Blob([dxfLines.join('\r\n')], {
    type: 'application/dxf;charset=utf-8;',
  });
  triggerDownload(blob, filename);
}
