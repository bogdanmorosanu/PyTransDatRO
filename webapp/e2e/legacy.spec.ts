import { test, expect } from './fixtures';

test.describe('Legacy Usage Statistics Subpage', () => {
  test('navigates to Usage -> Legacy tab and renders legacy stats with unified archetypes', async ({ page }) => {
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

    // 5. Verify 5 centered Archetype 8.1 KPI cards in five-cards grid
    const kpiBanner = page.locator('.telemetry-kpi-banner.five-cards');
    await expect(kpiBanner).toBeVisible();
    const kpiCards = kpiBanner.locator('.telemetry-kpi-card');
    await expect(kpiCards).toHaveCount(5);

    await expect(kpiCards.nth(0).locator('.telemetry-kpi-label')).toHaveText('Total Points');
    await expect(kpiCards.nth(1).locator('.telemetry-kpi-label')).toHaveText('Web App Points');
    await expect(kpiCards.nth(2).locator('.telemetry-kpi-label')).toHaveText('REST API Points');
    await expect(kpiCards.nth(3).locator('.telemetry-kpi-label')).toHaveText('Stereo 70 Points');
    await expect(kpiCards.nth(4).locator('.telemetry-kpi-label')).toHaveText('Stereo 30 Points');

    // 6. Verify date inputs show the legacy bounds (starting 2011-01-09, ending 2026-10-14)
    const startDateInput = page.locator('#legacy-date-start');
    const endDateInput = page.locator('#legacy-date-end');
    await expect(startDateInput).toBeVisible();
    await expect(endDateInput).toBeVisible();
    await expect(startDateInput).toHaveValue('2011-01-09');
    await expect(endDateInput).toHaveValue('2026-10-14');

    // 7. Verify Date Limits with Info Archetype (Full YYYY-MM-DD Resolution)
    const dateLimits = page.locator('.telemetry-date-limits-cluster');
    await expect(dateLimits).toBeVisible();
    await expect(dateLimits.locator('.info-label').nth(0)).toHaveText('Archive Start');
    await expect(dateLimits.locator('.info-field-readonly').nth(0)).toHaveText('2011-01-09');
    await expect(dateLimits.locator('.info-label').nth(1)).toHaveText('Archive End');
    await expect(dateLimits.locator('.info-field-readonly').nth(1)).toHaveText('2026-10-14');

    // 8. Verify source toggles and op toggles as Enclosed Multi-Select Trays (Archetype 2.3)
    const sourcesTray = page.locator('.filters-left-group .telemetry-btn-toggle-group.multi-select');
    await expect(sourcesTray).toBeVisible();
    const appSourceBtn = page.locator('[data-legacy-source="app"]');
    const serviceSourceBtn = page.locator('[data-legacy-source="service"]');
    await expect(appSourceBtn).toBeVisible();
    await expect(serviceSourceBtn).toBeVisible();

    const opsTray = page.locator('.filters-right-group .telemetry-btn-toggle-group.multi-select');
    await expect(opsTray).toBeVisible();
    const st70OpBtn = page.locator('[data-legacy-op="st70"]');
    const st30OpBtn = page.locator('[data-legacy-op="st30"]');
    await expect(st70OpBtn).toBeVisible();
    await expect(st30OpBtn).toBeVisible();

    // 9. Verify chart slot exists
    const chartSlot = page.locator('#legacy-chart-slot');
    await expect(chartSlot).toBeAttached();

    // 10. Verify stationary top-left legacy total toggle button (Icon-free, 0 layout shift)
    const totalToggleBtn = page.locator('#legacy-toggle-total');
    await expect(totalToggleBtn).toBeVisible();
    await expect(totalToggleBtn).toHaveText('Total');
    await expect(totalToggleBtn).not.toHaveText(/✓/);
    await expect(totalToggleBtn).toHaveClass(/active/);

    const initialBox = await totalToggleBtn.boundingBox();
    expect(initialBox).not.toBeNull();

    // Click to turn OFF total
    await totalToggleBtn.click();
    await expect(totalToggleBtn).not.toHaveClass(/active/);
    const afterOffBox = await totalToggleBtn.boundingBox();
    expect(afterOffBox).not.toBeNull();
    expect(Math.abs((afterOffBox?.x ?? 0) - (initialBox?.x ?? 0))).toBeLessThanOrEqual(1);

    // Click to turn ON total
    await totalToggleBtn.click();
    await expect(totalToggleBtn).toHaveClass(/active/);
    const afterOnBox = await totalToggleBtn.boundingBox();
    expect(afterOnBox).not.toBeNull();
    expect(Math.abs((afterOnBox?.x ?? 0) - (initialBox?.x ?? 0))).toBeLessThanOrEqual(1);
  });
});
