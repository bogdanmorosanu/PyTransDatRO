import type { Page } from '@playwright/test';
import { expect, test } from './fixtures';
import {
  DISPLAY_DECIMALS,
  REFERENCE_POINTS,
  matchesAtDisplayPrecision,
  type ReferencePoint,
} from './referencePoints';

const EMPTY_OUTPUT = '—';

/** Wait for the automatic first calculation (default ETRS89 -> Stereo 70) so later clicks cannot race it. */
async function waitForInitialResult(page: Page): Promise<void> {
  await page.goto('/');
  const results = page.getByRole('group', { name: 'Computed coordinates' });
  await expect(results.getByLabel('Northing')).not.toHaveValue(EMPTY_OUTPUT);
}

async function selectStereoToEtrs(page: Page): Promise<void> {
  await page.getByRole('radio', { name: 'Stereo 70 → ETRS89' }).click();
  await expect(page.getByRole('radio', { name: 'Stereo 70 → ETRS89' })).toBeChecked();
  const results = page.getByRole('group', { name: 'Computed coordinates' });
  await expect(results.getByLabel('Latitude')).not.toHaveValue(EMPTY_OUTPUT);
}

async function enterAndTransform(
  page: Page,
  values: Record<string, string>,
): Promise<void> {
  const inputs = page.getByRole('group', { name: 'Input coordinates (or click on map)' });
  for (const [label, value] of Object.entries(values)) {
    await inputs.getByLabel(label).fill(value);
  }
  await page.getByRole('button', { name: 'Transform', exact: true }).click();
}

/** Formats decimal degrees as D°MM'SS.SSSSS" (the format accepted by dmsParser). */
function toDms(deg: number): string {
  const totalSeconds = Math.round(deg * 3600 * 1e5) / 1e5;
  const d = Math.floor(totalSeconds / 3600);
  const m = Math.floor((totalSeconds - d * 3600) / 60);
  const s = totalSeconds - d * 3600 - m * 60;
  return `${d}°${String(m).padStart(2, '0')}'${s.toFixed(5).padStart(8, '0')}"`;
}

async function readNumber(locator: ReturnType<Page['getByLabel']>): Promise<string> {
  return locator.inputValue();
}

test.describe('Stereo 70 -> ETRS89 (decimal degrees, 3D)', () => {
  for (const ref of REFERENCE_POINTS) {
    test(`reference point ${ref.id} matches ANCPI at display precision, map not moved`, async ({
      page,
      externalRequests,
      mapView,
    }: {
      page: Page;
      externalRequests: string[];
      mapView: () => ReturnType<typeof import('./fixtures').readMapView>;
    }) => {
      void externalRequests;
      await waitForInitialResult(page);
      await selectStereoToEtrs(page);
      const before = await mapView();

      await enterAndTransform(page, {
        Northing: ref.north.toString(),
        Easting: ref.east.toString(),
        Altitude: ref.height.toString(),
      });

      const results = page.getByRole('group', { name: 'Computed coordinates' });
      await expect
        .poll(async () =>
          matchesAtDisplayPrecision(
            await readNumber(results.getByLabel('Latitude')),
            ref.lat,
            DISPLAY_DECIMALS.degrees,
          ),
        )
        .toBe(true);
      await expect
        .poll(async () =>
          matchesAtDisplayPrecision(
            await readNumber(results.getByLabel('Longitude')),
            ref.lon,
            DISPLAY_DECIMALS.degrees,
          ),
        )
        .toBe(true);
      await expect
        .poll(async () =>
          matchesAtDisplayPrecision(
            await readNumber(results.getByLabel('Height', { exact: true })),
            ref.ellipsoidalHeight,
            DISPLAY_DECIMALS.metres,
          ),
        )
        .toBe(true);
      await expect(results.getByText(/In Grid/)).toBeVisible();

      const after = await mapView();
      expect(after.center).toEqual(before.center);
      expect(after.zoom).toBe(before.zoom);
    });
  }
});

test.describe('ETRS89 (DMS input) -> Stereo 70', () => {
  const ref: ReferencePoint = REFERENCE_POINTS[0];

  test('reference point 1 entered as DMS matches ANCPI at metre precision', async ({
    page,
    externalRequests,
  }) => {
    void externalRequests;
    await waitForInitialResult(page);
    await page.getByLabel('Angular unit:').selectOption({ label: 'DMS' });

    await enterAndTransform(page, {
      Latitude: toDms(ref.lat),
      Longitude: toDms(ref.lon),
      Height: ref.ellipsoidalHeight.toString(),
    });

    const results = page.getByRole('group', { name: 'Computed coordinates' });
    await expect
      .poll(async () =>
        matchesAtDisplayPrecision(
          await readNumber(results.getByLabel('Northing')),
          ref.north,
          DISPLAY_DECIMALS.metres,
        ),
      )
      .toBe(true);
    await expect
      .poll(async () =>
        matchesAtDisplayPrecision(
          await readNumber(results.getByLabel('Easting')),
          ref.east,
          DISPLAY_DECIMALS.metres,
        ),
      )
      .toBe(true);
    await expect
      .poll(async () =>
        matchesAtDisplayPrecision(
          await readNumber(results.getByLabel('Altitude')),
          ref.height,
          DISPLAY_DECIMALS.metres,
        ),
      )
      .toBe(true);
  });
});

test.describe('points without a transformation', () => {
  test('point outside the grid envelope warns and presents no result, map unchanged', async ({
    page,
    externalRequests,
    mapView,
  }) => {
    void externalRequests;
    await waitForInitialResult(page);
    await selectStereoToEtrs(page);
    const before = await mapView();

    // Far outside the grid envelope (tests/test_trans_ro.py out-of-grid cases).
    await enterAndTransform(page, { Northing: '100000', Easting: '100000', Altitude: '100' });

    await expect(
      page.getByText('Coordinate is outside the official TransDatRO coverage area.'),
    ).toBeVisible();
    const results = page.getByRole('group', { name: 'Computed coordinates' });
    await expect(results.getByText(/Outside Grid/)).toBeVisible();
    await expect(results.getByLabel('Latitude')).toHaveValue(EMPTY_OUTPUT);
    await expect(results.getByLabel('Longitude')).toHaveValue(EMPTY_OUTPUT);
    await expect(results.getByLabel('Height', { exact: true })).toHaveValue(EMPTY_OUTPUT);

    const after = await mapView();
    expect(after.center).toEqual(before.center);
    expect(after.zoom).toBe(before.zoom);
  });

  // KNOWN APP BUG (expected to fail): PointConverter.ts only treats the warning 'Out of grid' as
  // untransformable (e.g. line ~429 `isOutOfGrid`); for 'No data on grid' the API echoes the
  // untransformed input, which is then displayed as if it were the transformed output.
  test('point in a no-data zone is flagged and the untransformed input is not shown as output', async ({
    page,
    request,
    externalRequests,
  }) => {
    void externalRequests;
    // tests/test_trans_ro.py::test_st70_to_etrs89_just_out_ro_border_no_data
    const north = 549364;
    const east = 147553;

    const api = await request.post('/api/v1/transform/point', {
      data: { op: 'Stereo70ToETRS89', coos: [north, east, 100], unit: 'degrees' },
    });
    expect(api.ok()).toBe(true);
    expect((await api.json()).warning).toBe('No data on grid');

    await waitForInitialResult(page);
    await selectStereoToEtrs(page);
    await enterAndTransform(page, {
      Northing: String(north),
      Easting: String(east),
      Altitude: '100',
    });

    const results = page.getByRole('group', { name: 'Computed coordinates' });
    await expect(results.getByText(/No Grid Data/)).toBeVisible();
    // The echoed input (549364 / 147553) must not be presented as latitude / longitude.
    await expect(results.getByLabel('Latitude')).toHaveValue(EMPTY_OUTPUT);
    await expect(results.getByLabel('Longitude')).toHaveValue(EMPTY_OUTPUT);
    // A point that could not be transformed is not "valid coordinates".
    await expect(results.getByText('Valid coordinates')).toHaveCount(0);
  });
});
