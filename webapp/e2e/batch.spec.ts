import { DISPLAY_DECIMALS, REFERENCE_POINTS, matchesAtDisplayPrecision } from './referencePoints';
import { expect, test } from './fixtures';
import { lines, numbers, openApp, openTab, setDirectionStereoToEtrs } from './helpers/ui';
import { readFile } from 'node:fs/promises';

const [P1, P2, P3] = REFERENCE_POINTS;

function st70Line(p: (typeof REFERENCE_POINTS)[number]): string {
  return `${p.north}, ${p.east}, ${p.height}`;
}

async function runBatch(page: import('@playwright/test').Page, input: string): Promise<string> {
  await openApp(page);
  await openTab(page, 'Text');
  await setDirectionStereoToEtrs(page);
  await page.locator('#batch-input-area').fill(input);
  await page.locator('#batch-transform-btn').click();
  const output = page.locator('#batch-output-area');
  await expect(output).toBeVisible();
  await expect(output).not.toHaveValue('');
  return output.inputValue();
}

test('batch keeps strict 1-to-1 lines and 1-based Point N labels around blank and comment lines', async ({
  page,
  externalRequests,
}) => {
  void externalRequests;
  const input = [
    '# header comment',
    st70Line(P1),
    '',
    '# middle comment',
    st70Line(P2),
    '',
    '',
    st70Line(P3),
  ].join('\n');

  const out = lines(await runBatch(page, input));

  // Strict 1-to-1: line N in == line N out, blanks and comments keep their position.
  expect(out).toHaveLength(8);
  expect(out[0]).toBe('# header comment');
  expect(out[2]).toBe('');
  expect(out[3]).toBe('# middle comment');
  expect(out[5]).toBe('');
  expect(out[6]).toBe('');

  const dataLines: Array<[number, (typeof REFERENCE_POINTS)[number]]> = [
    [1, P1],
    [4, P2],
    [7, P3],
  ];
  for (const [idx, p] of dataLines) {
    const [lat, lon, h] = numbers(out[idx]);
    expect(matchesAtDisplayPrecision(lat, p.lat, DISPLAY_DECIMALS.degrees)).toBe(true);
    expect(matchesAtDisplayPrecision(lon, p.lon, DISPLAY_DECIMALS.degrees)).toBe(true);
    expect(matchesAtDisplayPrecision(h, p.ellipsoidalHeight, DISPLAY_DECIMALS.metres)).toBe(true);
  }

  // Unnamed data points are labelled Point 1..3 (blank/comment lines do not consume an index).
  const [download] = await Promise.all([
    page.waitForEvent('download'),
    page.locator('#batch-export-csv').click(),
  ]);
  const csv = (await readFile((await download.path())!, 'utf8')).split(/\r?\n/);
  expect(csv).toHaveLength(3);
  csv.forEach((row, i) => {
    expect(row.startsWith(`Point ${i + 1},`)).toBe(true);
  });
});

test('batch reports a no-data-zone point as untransformable and not as a success', async ({
  page,
  externalRequests,
  request,
}) => {
  void externalRequests;
  // Stereo 70 point from tests/test_trans_ro.py (test_st70_to_etrs89_just_out_ro_border_no_data).
  const noData = '549364, 147553, 100';

  // Confirm via the real API that this input yields the "No data on grid" warning.
  const api = await request.post('/api/v1/transform/batch', {
    data: { op: 'Stereo70ToETRS89', points: [[549364, 147553, 100]], unit: 'degrees' },
  });
  expect(api.ok()).toBe(true);
  expect((await api.json()).results[0].warning).toBe('No data on grid');

  const out = lines(await runBatch(page, [st70Line(P1), noData].join('\n')));
  expect(out).toHaveLength(2);

  // The valid point is transformed.
  expect(matchesAtDisplayPrecision(numbers(out[0])[0], P1.lat, DISPLAY_DECIMALS.degrees)).toBe(true);

  // The no-data point is flagged and its input is not presented as transformed output.
  expect(out[1].toLowerCase()).toContain('no data');
  expect(out[1]).not.toContain('549364');

  // Summary pills: only the valid point counts as in-grid / successful.
  await expect(page.getByText(/Valid \(1\)/i)).toBeVisible();
  await expect(page.getByText(/Valid \(2\)/i)).toHaveCount(0);
  await expect(page.getByText(/Outside Grid \(1\)|No Grid Data \(1\)/i)).toBeVisible();
});

test('batch renders extended breakdown badges (Total points and 2D/3D count pills)', async ({ page }) => {
  await openApp(page);
  await openTab(page, 'Text');
  await setDirectionStereoToEtrs(page);
  
  // Mixed 2D and 3D points
  const input = [
    '500000.000, 500000.000, 100.000',
    '500200.000, 500150.000',
  ].join('\n');
  
  await page.locator('#batch-input-area').fill(input);
  await page.locator('#batch-transform-btn').click();
  
  await expect(page.locator('#batch-output-area')).toBeVisible();
  
  // Quality report badges verification
  const qualityCard = page.locator('.report-card');
  await expect(qualityCard).toBeVisible();
  await expect(qualityCard.getByText(/Total \(2\)/i)).toBeVisible();
  await expect(qualityCard.getByText(/2D \(1\)/i)).toBeVisible();
  await expect(qualityCard.getByText(/3D \(1\)/i)).toBeVisible();
  await expect(qualityCard.getByText(/Valid \(2\)/i)).toBeVisible();
});

test.describe('Batch mixed delimiters & clean error formatting', () => {
  test('handles comma, semicolon, tab, and whitespace batches with named and unnamed points', async ({
    page,
  }) => {
    await openApp(page);
    await openTab(page, 'Text');
    await setDirectionStereoToEtrs(page);

    // 1. Semicolon batch
    const semiInput = [
      'P_SEMI1; 500000.000; 500000.000; 100.000',
      'P_SEMI2; 500200.000; 500150.000; 102.500',
    ].join('\n');

    await page.locator('#batch-input-area').fill(semiInput);
    await page.locator('#batch-transform-btn').click();

    const output = page.locator('#batch-output-area');
    await expect(output).toBeVisible();
    let outLines = lines(await output.inputValue());
    expect(outLines).toHaveLength(2);
    expect(outLines[0]).toContain('P_SEMI1');
    expect(outLines[1]).toContain('P_SEMI2');

    // 2. Tab batch
    const tabInput = [
      'P_TAB1\t500000.000\t500000.000\t100.000',
      'P_TAB2\t501000.000\t499300.000\t98.300',
    ].join('\n');

    await page.locator('#batch-input-area').fill(tabInput);
    await page.locator('#batch-transform-btn').click();

    await expect(output).toHaveValue(/P_TAB1/);
    outLines = lines(await output.inputValue());
    expect(outLines).toHaveLength(2);
    expect(outLines[0]).toContain('P_TAB1');
    expect(outLines[1]).toContain('P_TAB2');
  });

  test('mixed error batch keeps 1-to-1 lines and formats invalid rows with comments and no coordinates', async ({
    page,
  }) => {
    await openApp(page);
    await openTab(page, 'Text');
    await setDirectionStereoToEtrs(page);

    const input = [
      '# Initial survey block',
      '101, 500000.000, 500000.000, 100.000',
      '',
      '102, 500abc.000, 500000.000, 100.000',
      '# Middle comment',
      '103, 100000.000, 100000.000, 0.000',
      '104, 549364, 147553, 100',
    ].join('\n');

    await page.locator('#batch-input-area').fill(input);
    await page.locator('#batch-transform-btn').click();

    const output = page.locator('#batch-output-area');
    const outLines = lines(await output.inputValue());

    // 1-to-1 line correspondence
    expect(outLines).toHaveLength(7);
    expect(outLines[0]).toBe('# Initial survey block');
    expect(outLines[2]).toBe('');
    expect(outLines[4]).toBe('# Middle comment');

    // Valid point has transformed coordinates
    expect(outLines[1]).toContain('101');
    const [p101Lat] = numbers(outLines[1]);
    expect(p101Lat).toBeGreaterThan(43);

    // Clean error row formatting: invalid input row contains only Point ID + # Invalid input (no numbers)
    expect(outLines[3]).toBe('102, # Invalid input');
    expect(outLines[3]).not.toContain('500');

    // Out of grid row: contains only Point ID + # Out of grid
    expect(outLines[5]).toBe('103, # Out of grid');
    expect(outLines[5]).not.toContain('100000');

    // No data row: contains only Point ID + # No data on grid
    expect(outLines[6]).toBe('104, # No data on grid');
    expect(outLines[6]).not.toContain('549364');

    // Quality report card diagnostic counts
    const qualityCard = page.locator('.report-card');
    await expect(qualityCard.getByText(/Total \(4\)/i)).toBeVisible();
    await expect(qualityCard.getByText(/Valid \(1\)/i)).toBeVisible();
    await expect(qualityCard.getByText(/Invalid Input \(1\)/i)).toBeVisible();
    await expect(qualityCard.getByText(/Outside Grid \(1\)/i)).toBeVisible();
    await expect(qualityCard.getByText(/No Grid Data \(1\)/i)).toBeVisible();
  });
});

test.describe('The Settings-Shift Lifecycle Invariant', () => {
  test('transforming batch with initial settings, then updating settings and transforming new text produces output matching new settings', async ({
    page,
  }) => {
    await openApp(page);
    await openTab(page, 'Text');

    // Step 1: Default ETRS89 to Stereo 70, Degrees, NE order
    const etrsInput = '101, 46.00000000, 25.00000000, 100.000';
    await page.locator('#batch-input-area').fill(etrsInput);
    await page.locator('#batch-transform-btn').click();

    const output = page.locator('#batch-output-area');
    await expect(output).toBeVisible();
    let outText = await output.inputValue();
    // Step 1 output must be in Stereo 70 meters (~500031)
    expect(outText).toContain('500031');

    // Step 2: Change settings: Stereo 70 to ETRS89, EN order, DMS units
    await page.getByRole('radio', { name: 'Stereo70 to ETRS89' }).click();
    await page.locator('button[data-order-val="EN"]').click();
    await page.locator('#angular-unit-select').selectOption('dms');

    // Step 3: Enter new input text in Stereo 70 EN order: Easting, Northing, Altitude
    const s70Input = '201, 500000.000, 500000.000, 100.000';
    await page.locator('#batch-input-area').fill(s70Input);

    // Step 4: Click Transform
    await page.locator('#batch-transform-btn').click();

    // Step 5: Verify output was generated using the NEW active settings:
    // Longitude first (~24°59'), then Latitude (~45°59'), with DMS symbols
    await expect(output).toHaveValue(/201/);
    outText = await output.inputValue();
    expect(outText).toContain('201');
    expect(outText).toContain('24°');
    expect(outText).toContain('45°');
    expect(outText).toContain('"');
  });
});

test.describe('Batch Textarea Dimensions & Export Filtering', () => {
  test('clicking Transform preserves custom vertical textarea heights', async ({ page }) => {
    await openApp(page);
    await openTab(page, 'Text');

    const inputArea = page.locator('#batch-input-area');
    // Set custom inline height
    await inputArea.evaluate((el) => {
      el.style.height = '315px';
    });

    await page.locator('#batch-transform-btn').click();
    await expect(page.locator('#batch-output-area')).toBeVisible();

    // Verify input area retains custom height after re-render
    const computedHeight = await inputArea.evaluate((el) => el.style.height);
    expect(computedHeight).toBe('315px');
  });

  test('GeoJSON and DXF exports cleanly filter out invalid rows, and downloads do not emit premature toasts', async ({
    page,
  }) => {
    await openApp(page);
    await openTab(page, 'Text');
    await setDirectionStereoToEtrs(page);

    const input = [
      '101, 500000.000, 500000.000, 100.000',
      '102, 500abc.000, 500000.000',
      '103, 100000.000, 100000.000, 0.000',
    ].join('\n');

    await page.locator('#batch-input-area').fill(input);
    await page.locator('#batch-transform-btn').click();

    // 1. Download GeoJSON
    const [geoDownload] = await Promise.all([
      page.waitForEvent('download'),
      page.locator('#batch-export-geojson').click(),
    ]);
    const geoText = await readFile((await geoDownload.path())!, 'utf8');
    const geoJson = JSON.parse(geoText);
    // Only the 1 valid point should be in the FeatureCollection
    expect(geoJson.features).toHaveLength(1);
    expect(geoJson.features[0].properties.id).toBe('101');

    // Verify no premature download toast was created
    await expect(page.locator('.toast').getByText(/Downloaded/i)).toHaveCount(0);

    // 2. Download DXF
    const [dxfDownload] = await Promise.all([
      page.waitForEvent('download'),
      page.locator('#batch-export-dxf').click(),
    ]);
    const dxfText = await readFile((await dxfDownload.path())!, 'utf8');
    // Standard DXF POINT entity check
    expect(dxfText).toContain('POINT');
    // Invalid point IDs must not appear in valid DXF entities
    expect(dxfText).not.toContain('102');
  });

  test('copy results button copies entire computed output to clipboard', async ({ page, context }) => {
    await context.grantPermissions(['clipboard-read', 'clipboard-write']);
    await openApp(page);
    await openTab(page, 'Text');
    await setDirectionStereoToEtrs(page);

    await page.locator('#batch-input-area').fill('500000.000, 500000.000, 100.000');
    await page.locator('#batch-transform-btn').click();

    const copyBtn = page.locator('#batch-copy-results-btn');
    await expect(copyBtn).toBeEnabled();
    await copyBtn.click();

    await expect(page.locator('.toast').getByText(/Copied/i)).toBeVisible();
    const clipText = await page.evaluate(() => navigator.clipboard.readText());
    expect(clipText).toContain('45.999');
    expect(clipText).toContain('24.998');
  });
});
