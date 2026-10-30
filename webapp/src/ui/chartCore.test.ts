import { describe, it, expect } from 'vitest';
import {
  buildSmoothSvgPath,
  calculateNearestXIndex,
  calculateXPosition,
} from './chartCore';

describe('chartCore', () => {
  describe('buildSmoothSvgPath', () => {
    it('returns empty string for empty points array', () => {
      expect(buildSmoothSvgPath([])).toBe('');
    });

    it('returns M command for single point', () => {
      expect(buildSmoothSvgPath([{ x: 10, y: 20 }])).toBe('M 10 20');
    });

    it('generates cubic Bézier commands for multiple points', () => {
      const points = [
        { x: 0, y: 100 },
        { x: 50, y: 50 },
        { x: 100, y: 0 },
      ];
      const path = buildSmoothSvgPath(points);
      expect(path.startsWith('M 0 100')).toBe(true);
      expect(path.includes('C')).toBe(true);
      expect(path.endsWith('100 0')).toBe(true);
    });
  });

  describe('calculateNearestXIndex', () => {
    it('returns 0 for single item list', () => {
      expect(calculateNearestXIndex(250, 50, 100, 1)).toBe(0);
    });

    it('computes correct nearest index on evenly spaced points', () => {
      // paddingLeft: 50, stepX: 100, 5 items (indices 0..4: x = 50, 150, 250, 350, 450)
      expect(calculateNearestXIndex(45, 50, 100, 5)).toBe(0);
      expect(calculateNearestXIndex(140, 50, 100, 5)).toBe(1);
      expect(calculateNearestXIndex(160, 50, 100, 5)).toBe(1);
      expect(calculateNearestXIndex(360, 50, 100, 5)).toBe(3);
    });

    it('clamps to [0, itemCount - 1]', () => {
      expect(calculateNearestXIndex(0, 50, 100, 5)).toBe(0);
      expect(calculateNearestXIndex(999, 50, 100, 5)).toBe(4);
    });
  });

  describe('calculateXPosition', () => {
    it('centers single item in chartWidth', () => {
      expect(calculateXPosition(0, 50, 100, 1, 800)).toBe(450); // 50 + 400
    });

    it('positions items proportionally by index and step', () => {
      expect(calculateXPosition(0, 50, 100, 5, 800)).toBe(50);
      expect(calculateXPosition(2, 50, 100, 5, 800)).toBe(250);
      expect(calculateXPosition(4, 50, 100, 5, 800)).toBe(450);
    });
  });
});
