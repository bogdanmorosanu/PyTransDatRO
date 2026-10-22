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

  await page.getByRole('button', { name: 'Transform Files', exact: true }).click();
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
