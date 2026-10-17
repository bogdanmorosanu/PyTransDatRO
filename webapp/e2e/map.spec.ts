import type { Page } from '@playwright/test';
import { REFERENCE_POINTS } from './referencePoints';
import { expect, test, type MapViewSnapshot } from './fixtures';
import { openApp, openTab, setDirectionStereoToEtrs } from './helpers/ui';

// Romania's envelope per .agents/rules/webapp-design.md section 7.
const ROMANIA = { south: 43.5, west: 20.2, north: 48.3, east: 29.8 };

async function transformBatch(page: Page): Promise<void> {
  await openTab(page, 'Text');
  await setDirectionStereoToEtrs(page);
  const input = REFERENCE_POINTS.map((p) => `${p.north}, ${p.east}, ${p.height}`).join('\n');
  await page.getByRole('textbox', { name: 'Input Coordinates' }).fill(input);
  await page.getByRole('button', { name: 'Transform', exact: true }).click();
  await expect(page.getByRole('textbox', { name: 'Computed Results' })).toBeVisible();
}

/** Waits until the map is idle: no zoom animation class and identical snapshots across two 250ms intervals. */
async function settled(page: Page, mapView: () => Promise<MapViewSnapshot>): Promise<MapViewSnapshot> {
  let previous = '';
  let stableCount = 0;
  await expect
    .poll(
      async () => {
        const animating = await page.locator('#map.leaflet-zoom-anim').count();
        const snap = JSON.stringify(await mapView());
        stableCount = animating === 0 && snap === previous ? stableCount + 1 : 0;
        previous = snap;
        return stableCount >= 2; // unchanged across two consecutive 250ms intervals
      },
      { intervals: [250], timeout: 15_000 },
    )
    .toBe(true);
  return mapView();
}

/** Waits until the view reaches the target (guards against sampling before the pan starts), then settles. */
async function settledAt(
  page: Page,
  mapView: () => Promise<MapViewSnapshot>,
  reached: (view: MapViewSnapshot) => boolean,
): Promise<MapViewSnapshot> {
  await expect.poll(async () => reached(await mapView()), { timeout: 15_000 }).toBe(true);
  return settled(page, mapView);
}

function expectCenterInRomania(view: MapViewSnapshot): void {
  expect(view.center.lat).toBeGreaterThanOrEqual(ROMANIA.south);
  expect(view.center.lat).toBeLessThanOrEqual(ROMANIA.north);
  expect(view.center.lng).toBeGreaterThanOrEqual(ROMANIA.west);
  expect(view.center.lng).toBeLessThanOrEqual(ROMANIA.east);
}

test('a successful batch transform does not change map center or zoom', async ({
  page,
  mapView,
}) => {
  await openApp(page);
  const before = await mapView();

  await transformBatch(page);
  await expect.poll(async () => (await mapView()).markersByLayer.text).toBe(REFERENCE_POINTS.length);

  const after = await mapView();
  expect(after.center).toEqual(before.center);
  expect(after.zoom).toBe(before.zoom);
});

test('a successful single-point transform does not change map center or zoom', async ({
  page,
  mapView,
}) => {
  await openApp(page);
  const before = await mapView();

  // Point tab is the default tab.
  await setDirectionStereoToEtrs(page);
  const p = REFERENCE_POINTS[0];
  await page.locator('#pt-input-1').fill(String(p.north));
  await page.locator('#pt-input-2').fill(String(p.east));
  await page.getByRole('button', { name: 'Transform', exact: true }).click();
  await expect.poll(async () => (await mapView()).markersByLayer.point).toBeGreaterThan(0);

  const after = await mapView();
  expect(after.center).toEqual(before.center);
  expect(after.zoom).toBe(before.zoom);
});

// Boundary-Safe Guard (webapp-design.md section 7): the app's own Center on Map / Zoom to Extents
// actions are bounded by Romania's envelope and never pan outside it for out-of-grid points.
// Free user navigation (drag / zoom) is intentionally unrestricted and not tested here.
const OUT_OF_GRID = { north: 100000, east: 100000 };
const OOG_ROW = `${OUT_OF_GRID.north}, ${OUT_OF_GRID.east}, 100`;

/**
 * Transforms one point from the Point tab. The app ignores Transform while an earlier transform
 * (initial load, direction switch) is still in flight and shows no busy state, so the click is
 * retried until `arrived` (a short assertion that this point's result is shown) passes.
 */
async function transformSinglePoint(
  page: Page,
  north: number,
  east: number,
  arrived: () => Promise<void>,
): Promise<void> {
  await setDirectionStereoToEtrs(page);
  await page.locator('#pt-input-1').fill(String(north));
  await page.locator('#pt-input-2').fill(String(east));
  await expect(async () => {
    await page.getByRole('button', { name: 'Transform', exact: true }).click();
    await arrived();
  }).toPass({ timeout: 15_000 });
}

async function transformBatchText(page: Page, rows: string[]): Promise<void> {
  await openTab(page, 'Text');
  await setDirectionStereoToEtrs(page);
  await page.getByRole('textbox', { name: 'Input Coordinates' }).fill(rows.join('\n'));
  await page.getByRole('button', { name: 'Transform', exact: true }).click();
  await expect(page.getByRole('textbox', { name: 'Computed Results' })).toBeVisible();
}

test('single out-of-grid point: Center on Map and Zoom to Extents never pan outside Romania', async ({
  page,
  mapView,
}) => {
  await openApp(page);
  const before = await settled(page, mapView);
  await transformSinglePoint(page, OUT_OF_GRID.north, OUT_OF_GRID.east, () =>
    expect(
      page.getByRole('group', { name: 'Computed coordinates' }).getByText(/Outside Grid/),
    ).toBeVisible({ timeout: 2_000 }),
  );

  for (const name of ['Center on Map', 'Zoom to Extents']) {
    // No plottable in-grid point exists, so the action must be unavailable.
    await expect(page.getByRole('button', { name })).toBeDisabled();
    const view = await settled(page, mapView);
    expectCenterInRomania(view);
    expect(view.center).toEqual(before.center);
    expect(view.zoom).toBe(before.zoom);
  }
});

test('batch with out-of-grid points: centering and zoom to extents ignore them', async ({
  page,
  mapView,
}) => {
  await openApp(page);
  const valid = REFERENCE_POINTS[0];
  await transformBatchText(page, [OOG_ROW, `${valid.north}, ${valid.east}, ${valid.height}`, OOG_ROW]);

  for (const name of ['Center on Map', 'Zoom to Extents']) {
    await page.getByRole('button', { name }).click();
    const view = await settled(page, mapView);
    expectCenterInRomania(view);
    // Only the valid point drives the view.
    expect(Math.abs(view.center.lat - valid.lat)).toBeLessThan(0.01);
    expect(Math.abs(view.center.lng - valid.lon)).toBeLessThan(0.01);
  }
});

test('batch with only out-of-grid points: map actions do not move the view', async ({
  page,
  mapView,
}) => {
  await openApp(page);
  const before = await settled(page, mapView);
  await transformBatchText(page, [OOG_ROW, '110000, 120000, 100']);

  for (const name of ['Center on Map', 'Zoom to Extents']) {
    // No plottable in-grid point exists, so the action must be unavailable.
    await expect(page.getByRole('button', { name })).toBeDisabled();
    const view = await settled(page, mapView);
    expectCenterInRomania(view);
    expect(view.center).toEqual(before.center);
    expect(view.zoom).toBe(before.zoom);
  }
});

test('valid in-grid point: Center on Map keeps zoom, Zoom to Extents goes to level 14', async ({
  page,
  mapView,
}) => {
  await openApp(page);
  const p = REFERENCE_POINTS[0];
  const before = await settled(page, mapView);
  // The card already shows the default point's result from load; wait for this point's latitude.
  const latitude = page.getByRole('group', { name: 'Computed coordinates' }).getByLabel('Latitude');
  await transformSinglePoint(page, p.north, p.east, () =>
    expect
      .poll(async () => Math.abs(parseFloat(await latitude.inputValue()) - p.lat), {
        timeout: 2_000,
      })
      .toBeLessThan(1e-4),
  );

  await page.getByRole('button', { name: 'Center on Map' }).click();
  const near = (v: MapViewSnapshot, tol: number): boolean =>
    Math.abs(v.center.lat - p.lat) < tol && Math.abs(v.center.lng - p.lon) < tol;
  const centered = await settledAt(page, mapView, (v) => near(v, 0.02));
  expectCenterInRomania(centered);
  // At the default zoom one pixel is ~0.007 degrees, so allow a couple of pixels of rounding.
  expect(Math.abs(centered.center.lat - p.lat)).toBeLessThan(0.02);
  expect(Math.abs(centered.center.lng - p.lon)).toBeLessThan(0.02);
  expect(centered.zoom).toBe(before.zoom);

  await page.getByRole('button', { name: 'Zoom to Extents' }).click();
  const zoomed = await settledAt(page, mapView, (v) => v.zoom === 14 && near(v, 0.001));
  expectCenterInRomania(zoomed);
  expect(Math.abs(zoomed.center.lat - p.lat)).toBeLessThan(0.001);
  expect(Math.abs(zoomed.center.lng - p.lon)).toBeLessThan(0.001);
  expect(zoomed.zoom).toBe(14);
});
