import { readFile } from 'node:fs/promises';
import { DISPLAY_DECIMALS, REFERENCE_POINTS, matchesAtDisplayPrecision } from './referencePoints';
import { expect, test } from './fixtures';
import { openApp, openTab } from './helpers/ui';

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
