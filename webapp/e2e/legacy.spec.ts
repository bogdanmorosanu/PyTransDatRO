import { test, expect } from './fixtures';

test.describe('Legacy Usage Statistics Subpage', () => {
  test('navigates to Usage -> Legacy tab and renders legacy stats', async ({ page }) => {
    await page.goto('/');

    // 1. Click on Usage view button in Header (#nav-btn-telemetry)
    const usageNavBtn = page.locator('#nav-btn-telemetry');
    await expect(usageNavBtn).toBeVisible();
    await usageNavBtn.click();

    // 2. Telemetry / Usage page should be visible
    const telemetryPage = page.locator('.telemetry-page-view');
    await expect(telemetryPage).toBeVisible();

    // 3. Tab buttons should include Computation, Spatial, and Legacy
    const legacyTabBtn = page.locator('#tab-legacy-btn');
    await expect(legacyTabBtn).toBeVisible();
    await expect(legacyTabBtn).toContainText('Legacy');

    // 4. Click on Legacy tab
    await legacyTabBtn.click();
    await expect(legacyTabBtn).toHaveClass(/active/);

    // 5. Verify KPI cards for legacy metrics
    const kpiBanner = page.locator('.telemetry-kpi-banner');
    await expect(kpiBanner).toBeVisible();

    // 6. Verify date inputs show the legacy bounds (starting 2011-01-09, ending 2026-10-14)
    const startDateInput = page.locator('#legacy-date-start');
    const endDateInput = page.locator('#legacy-date-end');
    await expect(startDateInput).toBeVisible();
    await expect(endDateInput).toBeVisible();
    await expect(startDateInput).toHaveValue('2011-01-09');
    await expect(endDateInput).toHaveValue('2026-10-14');

    // 7. Verify source toggles and op toggles
    const appSourceBtn = page.locator('[data-legacy-source="app"]');
    const serviceSourceBtn = page.locator('[data-legacy-source="service"]');
    await expect(appSourceBtn).toBeVisible();
    await expect(serviceSourceBtn).toBeVisible();

    const st70OpBtn = page.locator('[data-legacy-op="st70"]');
    const st30OpBtn = page.locator('[data-legacy-op="st30"]');
    await expect(st70OpBtn).toBeVisible();
    await expect(st30OpBtn).toBeVisible();

    // 8. Verify chart slot exists
    const chartSlot = page.locator('#legacy-chart-slot');
    await expect(chartSlot).toBeAttached();

    // 9. Verify legacy total overlay toggle button exists and responds to clicks
    const totalToggleBtn = page.locator('#legacy-toggle-total');
    await expect(totalToggleBtn).toBeVisible();
    await expect(totalToggleBtn).toHaveClass(/active/);
    await totalToggleBtn.click();
    await expect(totalToggleBtn).not.toHaveClass(/active/);
    await totalToggleBtn.click();
    await expect(totalToggleBtn).toHaveClass(/active/);
  });
});
