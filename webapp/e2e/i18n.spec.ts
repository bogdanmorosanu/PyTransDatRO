import { en } from '../src/i18n/en';
import { ro } from '../src/i18n/ro';
import { expect, test } from './fixtures';
import { openApp } from './helpers/ui';

test('switching language en -> ro updates labels and persists across reload', async ({
  page,
  externalRequests,
}) => {
  void externalRequests;
  await openApp(page);

  const tabs = page.locator('.tabs-bar');
  await expect(tabs.getByRole('button', { name: en.tabs.point, exact: true })).toBeVisible();
  await expect(tabs.getByRole('button', { name: en.tabs.file, exact: true })).toBeVisible();
  await expect(page.getByRole('radio', { name: en.point.opETRSToS70 })).toBeVisible();
  await expect(page.getByText(en.point.unitLabel).first()).toBeVisible();

  await page.getByRole('button', { name: en.app.langToggle }).click();
  await page.getByRole('menuitemradio', { name: ro.header.langRo }).click();

  const assertRomanian = async (): Promise<void> => {
    await expect(tabs.getByRole('button', { name: ro.tabs.point, exact: true })).toBeVisible();
    await expect(tabs.getByRole('button', { name: ro.tabs.file, exact: true })).toBeVisible();
    await expect(page.getByRole('radio', { name: ro.point.opETRSToS70 })).toBeVisible();
    await expect(page.getByText(ro.point.unitLabel).first()).toBeVisible();
    await expect(page.getByRole('button', { name: ro.app.langToggle })).toBeVisible();
  };
  await assertRomanian();
  await expect(tabs.getByRole('button', { name: en.tabs.file, exact: true })).toHaveCount(0);

  await page.reload();
  await assertRomanian();
});

test('transformation coordinate boxes strictly remain standard dot decimal while telemetry adapts to locale', async ({ page }) => {
  await openApp(page);
  
  // Switch to Romanian
  await page.getByRole('button', { name: en.app.langToggle }).click();
  await page.getByRole('menuitemradio', { name: ro.header.langRo }).click();
  
  // Coordinate inputs must keep dot decimal separator
  const i1 = page.locator('#pt-input-1');
  const i2 = page.locator('#pt-input-2');
  await expect(i1).toHaveValue(/\d+\.\d+/);
  await expect(i2).toHaveValue(/\d+\.\d+/);
  expect(await i1.inputValue()).not.toContain(',');
  expect(await i2.inputValue()).not.toContain(',');
  
  // Navigate to Usage page
  await page.locator('#nav-btn-telemetry').click();
  await expect(page.locator('.telemetry-page-view')).toBeVisible();
  await expect(page.locator('.telemetry-kpi-banner')).toBeVisible();

  // Spatial tab percentage shares must format with comma in Romanian
  await page.locator('#tab-spatial-btn').click();
  await expect(page.locator('.spatial-quadrants-deck')).toBeVisible();
  const shareTexts = await page.locator('.quadrant-share').allTextContents();
  expect(shareTexts.length).toBeGreaterThan(0);
  for (const st of shareTexts) {
    expect(st).not.toContain('.');
  }
  
  // Return to Studio / Transform via Brand Link (verifying SPA navigation)
  await page.locator('.header-brand-link').click();
  await expect(page.locator('.workbench-layout')).toBeVisible();
});
