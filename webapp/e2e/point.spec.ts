import type { Page } from '@playwright/test';
import { expect, test } from './fixtures';
import {
  DISPLAY_DECIMALS,
  REFERENCE_POINTS,
  matchesAtDisplayPrecision,
  type ReferencePoint,
} from './referencePoints';

const EMPTY_OUTPUT = '';

/** Wait for the automatic first calculation (default ETRS89 -> Stereo 70) so later clicks cannot race it. */
async function waitForInitialResult(page: Page): Promise<void> {
  await page.goto('/');
  const results = page.locator('.result-card');
  await expect(results.getByLabel('Northing')).not.toHaveValue(EMPTY_OUTPUT);
}

async function selectStereoToEtrs(page: Page): Promise<void> {
  await page.getByRole('radio', { name: 'Stereo70 to ETRS89' }).click();
  await expect(page.getByRole('radio', { name: 'Stereo70 to ETRS89' })).toBeChecked();
  const results = page.locator('.result-card');
  await expect(results.getByLabel('Latitude')).not.toHaveValue(EMPTY_OUTPUT);
}

async function enterAndTransform(
  page: Page,
  values: Record<string, string>,
): Promise<void> {
  const inputs = page.locator('.input-card');
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

      const results = page.locator('.result-card');
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
      await expect(page.locator('.report-card').getByText(/Valid/)).toBeVisible();

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
    await page.locator('#angular-unit-select').selectOption('dms');

    await enterAndTransform(page, {
      Latitude: toDms(ref.lat),
      Longitude: toDms(ref.lon),
      Height: ref.ellipsoidalHeight.toString(),
    });

    const results = page.locator('.result-card');
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
    const results = page.locator('.result-card');
    await expect(page.locator('.report-card').getByText(/Outside Grid/)).toBeVisible();
    await expect(results.getByLabel('Latitude')).toHaveValue('');
    await expect(results.getByLabel('Longitude')).toHaveValue('');
    await expect(results.getByLabel('Height', { exact: true })).toHaveValue('');

    const after = await mapView();
    expect(after.center).toEqual(before.center);
    expect(after.zoom).toBe(before.zoom);
  });

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

    const results = page.locator('.result-card');
    await expect(page.locator('.report-card').getByText(/No Grid Data/)).toBeVisible();
    // The echoed input (549364 / 147553) must not be presented as latitude / longitude.
    await expect(results.getByLabel('Latitude')).toHaveValue('');
    await expect(results.getByLabel('Longitude')).toHaveValue('');
    // A point that could not be transformed is not "valid coordinates".
    await expect(page.locator('.report-card').getByText('Valid coordinates')).toHaveCount(0);
  });
});

test.describe('URL Deep-Linking Hydration', () => {
  test('hydrates op, unit, order, and coordinates from query parameters', async ({ page }) => {
    await page.goto('/?op=Stereo70ToETRS89&unit=dms&order=EN&c1=500000.000&c2=500000.000&c3=100.000');
    
    // Verify radio direction
    const radio = page.getByRole('radio', { name: 'Stereo70 to ETRS89' });
    await expect(radio).toBeChecked();
    
    // Verify angular unit select
    const unitSelect = page.locator('#angular-unit-select');
    await expect(unitSelect).toHaveValue('dms');
    
    // Verify order control
    const orderBtn = page.locator('button[data-order-val="EN"]');
    await expect(orderBtn).toHaveClass(/active/);
    
    // Verify input values
    const i1 = page.locator('#pt-input-1');
    const i2 = page.locator('#pt-input-2');
    const i3 = page.locator('#pt-input-3');
    await expect(i1).toHaveValue('500000.000');
    await expect(i2).toHaveValue('500000.000');
    await expect(i3).toHaveValue('100.000');
    
    // Verify computed output has DMS symbol
    const results = page.locator('.result-card');
    await expect(results.getByLabel('Latitude')).not.toHaveValue(EMPTY_OUTPUT);
    await expect(results.getByLabel('Latitude')).toHaveValue(/°/);
  });
});

test.describe('Negative inputs & syntax resilience', () => {
  test('typing non-numeric characters displays Invalid Input, empties output, and disables copy/map', async ({
    page,
  }) => {
    await waitForInitialResult(page);
    await selectStereoToEtrs(page);

    await enterAndTransform(page, {
      Northing: '5001a50.000',
      Easting: '500000.000',
      Altitude: '100.000',
    });

    const reportCard = page.locator('.report-card');
    await expect(reportCard.getByText(/Invalid Input/i)).toBeVisible();

    const results = page.locator('.result-card');
    await expect(results.getByLabel('Latitude')).toHaveValue('');
    await expect(results.getByLabel('Longitude')).toHaveValue('');
    await expect(results.getByLabel('Height', { exact: true })).toHaveValue('');

    // Ensure zero inline error banners exist
    await expect(page.locator('.error-banner')).toHaveCount(0);

    // Copy and map buttons must be disabled
    await expect(page.locator('#pt-copy-code-btn')).toBeDisabled();
    await expect(page.locator('#pt-center-map-btn')).toBeDisabled();
    await expect(page.locator('#pt-zoom-extents-btn')).toBeDisabled();
  });

  test('impossible latitude > 90 deg in ETRS89 mode is flagged as Invalid Angle (not Outside Grid)', async ({
    page,
  }) => {
    await waitForInitialResult(page);
    // ETRS89 to Stereo 70 is default direction
    await enterAndTransform(page, {
      Latitude: '95.50000000',
      Longitude: '25.00000000',
      Height: '100.000',
    });

    const reportCard = page.locator('.report-card');
    await expect(reportCard.getByText(/Invalid Angle/i)).toBeVisible();
    await expect(reportCard.getByText(/Outside Grid/i)).toHaveCount(0);

    const results = page.locator('.result-card');
    await expect(results.getByLabel('Northing')).toHaveValue('');
    await expect(results.getByLabel('Easting')).toHaveValue('');
  });

  test('pasting Stereo 70 numbers into ETRS89 Latitude is caught as Invalid Angle', async ({
    page,
  }) => {
    await waitForInitialResult(page);
    await enterAndTransform(page, {
      Latitude: '500000.000',
      Longitude: '500000.000',
      Height: '100.000',
    });

    const reportCard = page.locator('.report-card');
    await expect(reportCard.getByText(/Invalid Angle/i)).toBeVisible();
    const results = page.locator('.result-card');
    await expect(results.getByLabel('Northing')).toHaveValue('');
  });
});

test.describe('2D vs 3D Dimensionality & Quasigeoid height zeta', () => {
  test('empty 3rd coordinate results in 2D mode with empty quasigeoid field', async ({ page }) => {
    await waitForInitialResult(page);
    await selectStereoToEtrs(page);

    await enterAndTransform(page, {
      Northing: '500000.000',
      Easting: '500000.000',
      Altitude: '',
    });

    const reportCard = page.locator('.report-card');
    await expect(reportCard.getByText(/2D/i)).toBeVisible();
    await expect(reportCard.getByText(/Valid/i)).toBeVisible();

    const quasigeoidBadge = page.locator('.output-quasigeoid-row .info-field-badge');
    await expect(quasigeoidBadge).toHaveText('');
  });

  test('populated 3rd coordinate triggers 3D mode and calculates quasigeoid height with sign', async ({
    page,
  }) => {
    await waitForInitialResult(page);
    await selectStereoToEtrs(page);

    await enterAndTransform(page, {
      Northing: '500000.000',
      Easting: '500000.000',
      Altitude: '100.000',
    });

    const reportCard = page.locator('.report-card');
    await expect(reportCard.getByText(/3D/i)).toBeVisible();
    await expect(reportCard.getByText(/Valid/i)).toBeVisible();

    const quasigeoidBadge = page.locator('.output-quasigeoid-row .info-field-badge');
    await expect(quasigeoidBadge).not.toHaveText('');
    await expect(quasigeoidBadge).toHaveText(/^[+-]\d+\.\d{3}$/);
  });

  test('corrupted 3rd coordinate triggers Invalid Input', async ({ page }) => {
    await waitForInitialResult(page);
    await selectStereoToEtrs(page);

    await enterAndTransform(page, {
      Northing: '500000.000',
      Easting: '500000.000',
      Altitude: 'invalid_alt',
    });

    const reportCard = page.locator('.report-card');
    await expect(reportCard.getByText(/Invalid Input/i)).toBeVisible();
    const results = page.locator('.result-card');
    await expect(results.getByLabel('Latitude')).toHaveValue('');
  });
});

test.describe('Smart Delimiter Split, Keyboard Shortcuts & Copy', () => {
  test('pasting multi-coordinate text into Coordinate 1 splits across inputs and auto-calculates', async ({
    page,
  }) => {
    await waitForInitialResult(page);
    await selectStereoToEtrs(page);

    // Focus input 1 and dispatch a paste event with space-separated coordinates
    const i1 = page.locator('#pt-input-1');
    await i1.focus();
    await page.evaluate(() => {
      const input = document.querySelector('#pt-input-1') as HTMLInputElement;
      const pasteEvent = new ClipboardEvent('paste', {
        bubbles: true,
        cancelable: true,
        clipboardData: new DataTransfer(),
      });
      pasteEvent.clipboardData?.setData('text/plain', '500200.000 500150.000 102.500');
      input.dispatchEvent(pasteEvent);
    });

    // Check all 3 inputs are populated
    await expect(page.locator('#pt-input-1')).toHaveValue('500200.000');
    await expect(page.locator('#pt-input-2')).toHaveValue('500150.000');
    await expect(page.locator('#pt-input-3')).toHaveValue('102.500');

    // Auto-calculation should produce valid results
    const results = page.locator('.result-card');
    await expect(results.getByLabel('Latitude')).not.toHaveValue('');
    await expect(page.locator('.report-card').getByText(/Valid/i)).toBeVisible();
  });

  test('pressing Enter key in input 2 triggers immediate transformation', async ({ page }) => {
    await waitForInitialResult(page);
    await selectStereoToEtrs(page);

    await page.locator('#pt-input-1').fill('501000.000');
    const i2 = page.locator('#pt-input-2');
    await i2.fill('499300.000');
    await page.locator('#pt-input-3').fill('98.300');

    // Press Enter inside input 2
    await i2.press('Enter');

    const results = page.locator('.result-card');
    await expect(results.getByLabel('Latitude')).not.toHaveValue('');
    await expect(page.locator('.report-card').getByText(/Valid/i)).toBeVisible();
  });

  test('copy button copies comma-separated coordinates with visual checkmark feedback', async ({
    page,
    context,
  }) => {
    await context.grantPermissions(['clipboard-read', 'clipboard-write']);
    await waitForInitialResult(page);
    await selectStereoToEtrs(page);

    await enterAndTransform(page, {
      Northing: '500000.000',
      Easting: '500000.000',
      Altitude: '100.000',
    });

    const copyBtn = page.locator('#pt-copy-code-btn');
    await expect(copyBtn).toBeEnabled();
    await copyBtn.click();

    // Verify feedback checkmark icon or tooltip
    await expect(page.locator('.toast').getByText(/Copied/i)).toBeVisible();

    const clipboardText = await page.evaluate(() => navigator.clipboard.readText());
    // Should be comma-separated coordinates
    const tokens = clipboardText.split(',').map((s) => s.trim());
    expect(tokens.length).toBeGreaterThanOrEqual(2);
    expect(Number(tokens[0])).not.toBeNaN();
    expect(Number(tokens[1])).not.toBeNaN();
  });
});

test.describe('Coordinate Order and Angular Unit Synchronization', () => {
  test('switching coordinate order updates input/output labels and order', async ({ page }) => {
    await waitForInitialResult(page);
    await selectStereoToEtrs(page);

    // Initial NE order: input 1 is Northing, input 2 is Easting
    await expect(page.locator('#pt-input-1')).toHaveAttribute('aria-label', 'Northing');
    await expect(page.locator('#pt-input-2')).toHaveAttribute('aria-label', 'Easting');

    // Switch to EN order
    await page.locator('button[data-order-val="EN"]').click();

    // In EN order: input 1 is Easting, input 2 is Northing
    await expect(page.locator('#pt-input-1')).toHaveAttribute('aria-label', 'Easting');
    await expect(page.locator('#pt-input-2')).toHaveAttribute('aria-label', 'Northing');

    // Output labels also swap to Longitude, Latitude
    const results = page.locator('.result-card');
    const outInputs = results.locator('.output-grid-3 input');
    await expect(outInputs.nth(0)).toHaveAttribute('aria-label', 'Longitude');
    await expect(outInputs.nth(1)).toHaveAttribute('aria-label', 'Latitude');
  });

  test('switching angular unit to DMS or Radians reformats output coordinates and updates example badge', async ({
    page,
  }) => {
    await waitForInitialResult(page);
    await selectStereoToEtrs(page);

    await enterAndTransform(page, {
      Northing: '500000.000',
      Easting: '500000.000',
      Altitude: '100.000',
    });

    const unitSelect = page.locator('#angular-unit-select');

    // Switch to DMS
    await unitSelect.selectOption('dms');
    const results = page.locator('.result-card');
    await expect(results.getByLabel('Latitude')).toHaveValue(/°/);
    await expect(results.getByLabel('Longitude')).toHaveValue(/°/);

    // Example badge in config bar reflects DMS
    const exampleBadge = page.locator('#angle-format-example-badge');
    await expect(exampleBadge).toContainText('°');

    // Switch to Radians
    await unitSelect.selectOption('radians');
    await expect(results.getByLabel('Latitude')).toHaveValue(/^0\.\d{10}$/);
  });
});
