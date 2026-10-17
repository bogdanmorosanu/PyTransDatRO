import { STUB_HEADER, expect, test } from './fixtures';

test('app loads with header, map and e2e handle, without errors or leaked requests', async ({
  page,
  consoleErrors,
  externalRequests,
  mapView,
}) => {
  const leaked: string[] = [];
  page.on('response', (res) => {
    const host = new URL(res.url()).hostname;
    const external = host !== '127.0.0.1' && host !== 'localhost';
    if (external && !res.url().startsWith('data:') && !(STUB_HEADER in res.headers())) {
      leaked.push(res.url());
    }
  });

  await page.goto('/');

  await expect(page).toHaveTitle(/PyTransDatRO/);
  await expect(page.locator('.brand-title-text')).toBeVisible();
  await expect(page.locator('#map')).toBeVisible();
  await expect(page.locator('#map .leaflet-tile-pane')).toBeAttached();

  const snapshot = await mapView();
  expect(snapshot.zoom).toBeGreaterThan(0);
  expect(snapshot.center.lat).toBeGreaterThan(43);
  expect(snapshot.center.lat).toBeLessThan(49);

  // Tiles were requested but only ever answered by the stub.
  await expect.poll(() => externalRequests.length).toBeGreaterThan(0);
  expect(leaked).toEqual([]);
  expect(consoleErrors).toEqual([]);
});
