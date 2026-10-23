import { readFile } from 'node:fs/promises';
import { DISPLAY_DECIMALS, REFERENCE_POINTS, matchesAtDisplayPrecision } from './referencePoints';
import { expect, test } from './fixtures';
import { openApp, openTab } from './helpers/ui';
import { uploadDiskFixtures, uploadMemoryFiles, readDownloadContent } from './helpers/files';

test('file upload transforms a CSV and the downloaded file is 1-to-1 with the input', async ({
  page,
  externalRequests,
}) => {
  void externalRequests;
  // ETRS89 -> Stereo 70 (default direction), id + lat + lon + ellipsoidal height, degrees.
  const input = REFERENCE_POINTS.map(
    (p) => `${p.id},${p.lat.toFixed(12)},${p.lon.toFixed(12)},${p.ellipsoidalHeight}`,
  );

  await openApp(page);
  await openTab(page, 'Files');
  await page.getByLabel('Drag & drop files here or click to browse').setInputFiles({
    name: 'points.csv',
    mimeType: 'text/csv',
    buffer: Buffer.from(input.join('\n') + '\n', 'utf8'),
  });
  await expect(page.locator('.file-queue-name')).toHaveText('points.csv');

  await page.getByRole('button', { name: 'Transform', exact: true }).click();
  await expect(page.locator('.file-status-badge.completed')).toBeVisible();
  await expect(page.getByText(/Valid \(3\)/i)).toBeVisible();

  const [download] = await Promise.all([
    page.waitForEvent('download'),
    page.locator('.file-download-btn').click(),
  ]);
  expect(download.suggestedFilename()).toBe('transformed_points.csv');
  const content = await readFile((await download.path())!, 'utf8');
  const out = content.split(/\r?\n/).filter((l) => l.length > 0);

  expect(out).toHaveLength(input.length);
  out.forEach((row, i) => {
    const p = REFERENCE_POINTS[i];
    const cols = row.split(',');
    // API file endpoint format: "<id>,<N .4f>,<E .4f>,<H .4f>"
    expect(cols).toHaveLength(4);
    expect(cols[0]).toBe(String(p.id));
    for (const c of cols.slice(1)) expect(c).toMatch(/^-?\d+\.\d{4}$/);
    expect(matchesAtDisplayPrecision(cols[1], p.north, DISPLAY_DECIMALS.metres)).toBe(true);
    expect(matchesAtDisplayPrecision(cols[2], p.east, DISPLAY_DECIMALS.metres)).toBe(true);
    expect(matchesAtDisplayPrecision(cols[3], p.height, DISPLAY_DECIMALS.metres)).toBe(true);
  });
});

test('file upload with ignore header preserves header on line 1 and counts valid data rows', async ({
  page,
  externalRequests,
}) => {
  void externalRequests;
  const header = 'PointID,Latitude,Longitude,Height';
  const dataRows = REFERENCE_POINTS.map(
    (p) => `${p.id},${p.lat.toFixed(12)},${p.lon.toFixed(12)},${p.ellipsoidalHeight}`,
  );
  const input = [header, ...dataRows];

  await openApp(page);
  await openTab(page, 'Files');

  // Toggle "Ignore header / first line" ON
  await page.locator('label[for="file-skip-header-toggle"]').click();
  await expect(page.locator('#file-skip-header-toggle')).toBeChecked();

  await page.getByLabel('Drag & drop files here or click to browse').setInputFiles({
    name: 'headed_points.csv',
    mimeType: 'text/csv',
    buffer: Buffer.from(input.join('\n') + '\n', 'utf8'),
  });

  await page.getByRole('button', { name: 'Transform', exact: true }).click();
  await expect(page.locator('.file-status-badge.completed')).toBeVisible();

  // Exactly 3 valid rows counted, 0 errors
  await expect(page.getByText(/Valid \(3\)/i)).toBeVisible();
  await expect(page.locator('.status-pill.danger')).toHaveCount(0);

  const [download] = await Promise.all([
    page.waitForEvent('download'),
    page.locator('.file-download-btn').click(),
  ]);
  const content = await readFile((await download.path())!, 'utf8');
  const out = content.split(/\r?\n/).filter((l) => l.length > 0);

  expect(out).toHaveLength(input.length);
  // First line in output must match the input header verbatim
  expect(out[0]).toBe(header);
  // Subsequent lines are transformed coordinates
  expect(out[1]).toContain(String(REFERENCE_POINTS[0].id));
});

test('multiple file upload shows per-file report breakdown with total line, and Remove All clears output and map', async ({
  page,
  externalRequests,
}) => {
  void externalRequests;
  const p1 = REFERENCE_POINTS[0];
  const p2 = REFERENCE_POINTS[1];
  const file1Content = `${p1.id},${p1.lat.toFixed(12)},${p1.lon.toFixed(12)},${p1.ellipsoidalHeight}\n`;
  const file2Content = `${p2.id},${p2.lat.toFixed(12)},${p2.lon.toFixed(12)},${p2.ellipsoidalHeight}\n`;

  await openApp(page);
  await openTab(page, 'Files');

  await page.getByLabel('Drag & drop files here or click to browse').setInputFiles([
    {
      name: 'survey_alpha.csv',
      mimeType: 'text/csv',
      buffer: Buffer.from(file1Content, 'utf8'),
    },
    {
      name: 'survey_beta.csv',
      mimeType: 'text/csv',
      buffer: Buffer.from(file2Content, 'utf8'),
    },
  ]);

  await expect(page.locator('.file-queue-item')).toHaveCount(2);

  // Transform both files
  await page.getByRole('button', { name: 'Transform', exact: true }).click();
  await expect(page.locator('.file-status-badge.completed')).toHaveCount(2);

  // Verify multi-line breakdown in Report Card 3
  const subrows = page.locator('.quality-report-subrow');
  await expect(subrows).toHaveCount(3); // 2 files + 1 total row
  await expect(subrows.nth(0)).toContainText('survey_alpha.csv');
  await expect(subrows.nth(1)).toContainText('survey_beta.csv');
  await expect(page.locator('.quality-report-subrow.total-row')).toContainText('Total');
  await expect(page.locator('.quality-report-subrow.total-row')).toContainText('Valid (2)');

  // Verify Remove All button exists and has correct label
  const removeAllBtn = page.locator('#file-clear-output-btn');
  await expect(removeAllBtn).toHaveText(/Remove All/i);

  // Map actions should be enabled before clearing
  await expect(page.locator('#file-center-map-btn')).toBeEnabled();

  // Click Remove All
  await removeAllBtn.click();

  // Transformed files queue must now be cleared and empty state displayed
  await expect(page.locator('.result-card .file-queue-item')).toHaveCount(0);
  await expect(page.locator('.file-empty-state')).toBeVisible();

  // Map actions should now be disabled because all points have been removed
  await expect(page.locator('#file-center-map-btn')).toBeDisabled();
});

test('staged multi-batch file transformation handles initial settings, settings shift, and subsequent file upload independently', async ({
  page,
  externalRequests,
}) => {
  void externalRequests;
  await openApp(page);
  await openTab(page, 'Files');

  // Stage 1: Upload File 1 with default settings (ETRS89 -> Stereo 70, Degrees)
  const p1 = REFERENCE_POINTS[0];
  const file1Content = `${p1.id},${p1.lat.toFixed(12)},${p1.lon.toFixed(12)},${p1.ellipsoidalHeight}\n`;
  await uploadMemoryFiles(page, [{ name: 'staged_etrs.csv', content: file1Content }]);

  // Transform File 1
  await page.getByRole('button', { name: 'Transform', exact: true }).click();
  await expect(page.locator('.file-status-badge.completed')).toHaveCount(1);

  // Stage 2: Shift Settings in UI: Direction Stereo70 to ETRS89, EN order, DMS units
  await page.getByRole('radio', { name: 'Stereo70 to ETRS89' }).click();
  await page.locator('button[data-order-val="EN"]').click();
  await page.locator('#angular-unit-select').selectOption('dms');

  // Stage 3: Upload File 2 with Stereo 70 coordinates: Easting, Northing, Altitude
  const file2Content = `201,500000.000,500000.000,100.000\n`;
  await uploadMemoryFiles(page, [{ name: 'staged_stereo.csv', content: file2Content }]);

  // Verify queue has File 1 completed in output queue, and File 2 pending in input queue
  await expect(page.locator('.result-card .file-queue-item')).toHaveCount(1);
  await expect(page.locator('.studio-card:not(.result-card) .file-queue-item')).toHaveCount(1);

  // Stage 4: Click Transform
  await page.getByRole('button', { name: 'Transform', exact: true }).click();
  await expect(page.locator('.result-card .file-queue-item')).toHaveCount(2);

  // Report card subrows show both files and total row
  const subrows = page.locator('.quality-report-subrow');
  await expect(subrows).toHaveCount(3); // 2 files + 1 total
  await expect(subrows.nth(0)).toContainText('staged_etrs.csv');
  await expect(subrows.nth(1)).toContainText('staged_stereo.csv');

  // Stage 5: Download File 2 (second completed item)
  const downloadBtns = page.locator('.result-card .file-download-btn');
  const [download2] = await Promise.all([
    page.waitForEvent('download'),
    downloadBtns.nth(1).click(),
  ]);
  const content2 = await readDownloadContent(download2);
  // Transformed File 2 must have DMS formatting and longitude first (~24°59')
  expect(content2).toContain('201');
  expect(content2).toContain('24°');
  expect(content2).toContain('45°');
  expect(content2).toContain('"');

  // Stage 6: Download File 1 (first completed item)
  const [download1] = await Promise.all([
    page.waitForEvent('download'),
    downloadBtns.nth(0).click(),
  ]);
  const content1 = await readDownloadContent(download1);
  // Transformed File 1 must have Stereo 70 Northing, Easting in meters
  expect(content1).toContain(String(p1.id));
  expect(content1).toContain(p1.north.toFixed(4).slice(0, 6));
});

test('corrupted file with mixed error rows preserves 1-to-1 row count and omits coordinates on invalid lines', async ({
  page,
  externalRequests,
}) => {
  void externalRequests;
  await openApp(page);
  await openTab(page, 'Files');

  await uploadDiskFixtures(page, ['mixed_errors.csv']);
  await page.getByRole('button', { name: 'Transform', exact: true }).click();

  await expect(page.locator('.file-status-badge.completed')).toBeVisible();

  // Download output
  const [download] = await Promise.all([
    page.waitForEvent('download'),
    page.locator('.file-download-btn').click(),
  ]);
  const content = await readDownloadContent(download);
  const outLines = content.split(/\r?\n/).filter((l) => l.length > 0);

  // Fixture has 8 non-empty lines (comments + valid + errors)
  expect(outLines).toHaveLength(7); // Line 3 was empty so 7 non-empty lines

  // Comment lines preserved
  expect(outLines[0]).toBe('# Header comment describing survey');
  // Valid lines transformed
  expect(outLines[1]).toContain('101');
  expect(outLines[1]).not.toContain('#');

  // Error lines omit coordinates and include comment
  expect(outLines[2]).toContain('102');
  expect(outLines[2]).toContain('#');

  // 103 was impossible angle 99 deg
  expect(outLines[4]).toContain('103');
  expect(outLines[4]).toContain('#');

  // 104 was out of grid
  expect(outLines[5]).toContain('104');
  expect(outLines[5]).toContain('#');

  // 105 was valid
  expect(outLines[6]).toContain('105');
  expect(outLines[6]).not.toContain('#');
});

test('file queue lifecycle: individual item removal before transform and clear queue purge input items', async ({
  page,
  externalRequests,
}) => {
  void externalRequests;
  await openApp(page);
  await openTab(page, 'Files');

  // Upload 2 files
  await uploadMemoryFiles(page, [
    { name: 'queue_a.csv', content: '1,46.0,25.0,100\n' },
    { name: 'queue_b.csv', content: '2,46.0,25.0,100\n' },
  ]);
  await expect(page.locator('.studio-card:not(.result-card) .file-queue-item')).toHaveCount(2);

  // Remove first item using its ✕ button
  await page.locator('.studio-card:not(.result-card) .file-queue-remove').first().click();
  await expect(page.locator('.studio-card:not(.result-card) .file-queue-item')).toHaveCount(1);
  await expect(page.locator('.file-queue-name')).toHaveText('queue_b.csv');

  // Add another file
  await uploadMemoryFiles(page, [
    { name: 'queue_c.csv', content: '3,46.0,25.0,100\n' },
  ]);
  await expect(page.locator('.studio-card:not(.result-card) .file-queue-item')).toHaveCount(2);

  // Click Clear Queue
  await page.locator('#file-clear-queue-btn').click();
  await expect(page.locator('.studio-card:not(.result-card) .file-queue-item')).toHaveCount(0);
});

test('volume warning modal appears when transforming > 2,000 points and handles sample preview', async ({
  page,
  externalRequests,
}) => {
  void externalRequests;
  await openApp(page);
  await openTab(page, 'Files');

  // Generate 2050 points
  const largeBatch = Array.from({ length: 2050 }, (_, i) => `${i + 1},46.00000000,25.00000000,100.000`).join('\n');
  await uploadMemoryFiles(page, [{ name: 'large_volume.csv', content: largeBatch }]);

  await page.getByRole('button', { name: 'Transform', exact: true }).click();

  // Volume Warning Modal should appear
  const modal = page.locator('#modal-root .modal-dialog');
  await expect(modal).toBeVisible();
  await expect(modal.locator('.modal-title')).toContainText('High Volume Point Warning');

  // Click Sample Preview button
  await modal.locator('#vw-sample-btn').click();
  await expect(modal).not.toBeVisible();

  // Status completed
  await expect(page.locator('.file-status-badge.completed')).toBeVisible();
});
