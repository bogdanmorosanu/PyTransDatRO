import { describe, it, expect } from 'vitest';
import {
  generateCSVContent,
  generateGeoJSONObject,
  generateDXFContent,
  type ExportPointItem,
} from './exporters';

describe('exporters', () => {
  const samplePlanarItems: ExportPointItem[] = [
    {
      pointId: 'P1',
      coords: [324123.4567, 588654.3214],
      isGeo: false,
    },
    {
      pointId: 'P2',
      coords: [324200.0, 588700.0, 85.1234],
      isGeo: false,
    },
    {
      pointId: 'P3_OOG',
      coords: [],
      warning: 'Out of grid',
      isGeo: false,
    },
    {
      pointId: null,
      coords: [324300.111, 588800.222],
      warning: 'Interpolation edge warning',
      isGeo: false,
    },
  ];

  const sampleGeoItems: ExportPointItem[] = [
    {
      pointId: 'GPS1',
      coords: [44.43231234, 26.10631234, 90.5],
      isGeo: true,
    },
    {
      pointId: 'GPS2',
      coords: [46.77123456, 23.62361234],
      isGeo: true,
    },
  ];

  describe('generateCSVContent', () => {
    it('generates CSV with comma delimiter and NE order for planar coordinates', () => {
      const csv = generateCSVContent(samplePlanarItems, ',', 'NE');
      const lines = csv.split('\r\n');

      expect(lines).toHaveLength(4);
      // P1: North, East (3 decimals)
      expect(lines[0]).toBe('P1,324123.457,588654.321');
      // P2: North, East, Elevation (3 decimals)
      expect(lines[1]).toBe('P2,324200.000,588700.000,85.123');
      // P3: Out of grid
      expect(lines[2]).toBe('P3_OOG, # Out of grid');
      // P4: No pointId, warning attached
      expect(lines[3]).toBe('324300.111,588800.222,# Interpolation edge warning');
    });

    it('generates CSV with EN order and custom semicolon delimiter', () => {
      const csv = generateCSVContent(samplePlanarItems, ';', 'EN');
      const lines = csv.split('\r\n');

      // P1: East, North
      expect(lines[0]).toBe('P1;588654.321;324123.457');
      // P2: East, North, Elevation
      expect(lines[1]).toBe('P2;588700.000;324200.000;85.123');
    });

    it('formats geographic coordinates to 8 decimal places', () => {
      const csv = generateCSVContent(sampleGeoItems, ',', 'NE');
      const lines = csv.split('\r\n');

      expect(lines[0]).toBe('GPS1,44.43231234,26.10631234,90.500');
      expect(lines[1]).toBe('GPS2,46.77123456,23.62361234');
    });

    it('skips items with fewer than 2 coordinates when not an out-of-grid warning', () => {
      const invalidItems: ExportPointItem[] = [
        { pointId: 'Bad', coords: [100], isGeo: false },
      ];
      const csv = generateCSVContent(invalidItems);
      expect(csv).toBe('');
    });

    it('preserves No data on grid warning lines with and without pointId', () => {
      const noDataItems: ExportPointItem[] = [
        { pointId: 'P_ND1', coords: [], warning: 'No data on grid', isGeo: true },
        { pointId: null, coords: [], warning: 'No data on grid', isGeo: true },
      ];
      const csv = generateCSVContent(noDataItems, ',', 'NE');
      expect(csv.split('\r\n')).toEqual([
        'P_ND1, # No data on grid',
        '# No data on grid',
      ]);
    });
  });

  describe('generateGeoJSONObject', () => {
    it('creates a valid GeoJSON FeatureCollection', () => {
      const geojson = generateGeoJSONObject(sampleGeoItems);

      expect(geojson.type).toBe('FeatureCollection');
      expect(geojson.features).toHaveLength(2);

      const f1 = geojson.features[0];
      expect(f1.type).toBe('Feature');
      expect(f1.id).toBe('GPS1');
      // GeoJSON standard: [Longitude, Latitude, Elevation]
      expect(f1.geometry.type).toBe('Point');
      expect(f1.geometry.coordinates).toEqual([26.10631234, 44.43231234, 90.5]);
      expect(f1.properties).toEqual({
        id: 'GPS1',
        elevation: 90.5,
        system: 'ETRS89',
      });

      const f2 = geojson.features[1];
      expect(f2.geometry.coordinates).toEqual([23.62361234, 46.77123456]);
      expect(f2.properties.elevation).toBeNull();
    });

    it('filters out warning items and malformed items', () => {
      const geojson = generateGeoJSONObject(samplePlanarItems);
      // P3 (warning: Out of grid) and P4 (warning: Interpolation edge warning) are excluded
      expect(geojson.features).toHaveLength(2);
      expect(geojson.features[0].id).toBe('P1');
      expect(geojson.features[1].id).toBe('P2');
      expect(geojson.features[0].properties.system).toBe('Stereo70');
    });
  });

  describe('generateDXFContent', () => {
    it('creates standard AutoCAD DXF structure with POINT and TEXT entities', () => {
      const dxf = generateDXFContent(samplePlanarItems);

      // Verify essential DXF section headers
      expect(dxf).toContain('SECTION\r\n2\r\nHEADER\r\n0\r\nENDSEC');
      expect(dxf).toContain('SECTION\r\n2\r\nENTITIES');
      expect(dxf).toContain('ENDSEC\r\n0\r\nEOF');

      // POINT entity for P1 on layer POINTS
      // X = 588654.3214, Y = 324123.4567, Z = 0.0000
      expect(dxf).toContain('POINT\r\n8\r\nPOINTS\r\n10\r\n588654.3214\r\n20\r\n324123.4567\r\n30\r\n0.0000');

      // TEXT entity for P1 label on layer LABELS
      expect(dxf).toContain('TEXT\r\n8\r\nLABELS');
      expect(dxf).toContain('1\r\nP1');

      // POINT entity for P2 with elevation Z = 85.1234
      expect(dxf).toContain('10\r\n588700.0000\r\n20\r\n324200.0000\r\n30\r\n85.1234');
      expect(dxf).toContain('1\r\nP2');

      // Warning items (P3 and P4) should not produce POINT entities
      expect(dxf).not.toContain('P3_OOG');
    });

    it('handles items without point IDs by creating POINT without TEXT', () => {
      const items: ExportPointItem[] = [
        {
          pointId: null,
          coords: [324000.0, 588000.0],
          isGeo: false,
        },
      ];
      const dxf = generateDXFContent(items);
      expect(dxf).toContain('POINT\r\n8\r\nPOINTS');
      expect(dxf).not.toContain('TEXT\r\n8\r\nLABELS');
    });
  });
});
