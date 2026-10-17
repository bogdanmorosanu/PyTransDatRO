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
  await expect(tabs.getByText(en.tabs.prefix)).toBeVisible();
  await expect(tabs.getByRole('button', { name: en.tabs.file, exact: true })).toBeVisible();
  await expect(page.getByRole('radio', { name: en.point.opETRSToS70 })).toBeVisible();
  await expect(page.getByText(en.point.unitLabel).first()).toBeVisible();

  await page.getByRole('button', { name: en.app.langToggle }).click();
  await page.getByRole('menuitemradio', { name: ro.header.langRo }).click();

  const assertRomanian = async (): Promise<void> => {
    await expect(tabs.getByText(ro.tabs.prefix)).toBeVisible();
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
