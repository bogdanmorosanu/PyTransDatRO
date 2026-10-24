import { expect, test } from './fixtures';

test.describe('Usage Computation Subpage', () => {
  test('renders refined header, centered KPI cards, info date limits, segmented trays, and stationary total toggle', async ({ page }) => {
    await page.goto('/');

    // 1. Navigate to Usage page
    const usageNavBtn = page.locator('#nav-btn-telemetry');
    await usageNavBtn.click();
    await expect(page.locator('.telemetry-page-container')).toBeVisible();

    // 2. Subtitle removed, tabs flush left
    await expect(page.locator('.telemetry-page-title')).toHaveCount(0);
    const computationTab = page.locator('#tab-timeline-btn');
    await expect(computationTab).toBeVisible();
    await expect(computationTab).toHaveClass(/active/);

    // 3. KPI Cards centered with sentence-case labels
    const kpiCards = page.locator('.telemetry-kpi-card');
    await expect(kpiCards).toHaveCount(3);
    await expect(kpiCards.nth(0).locator('.telemetry-kpi-label')).toHaveText('Points');
    await expect(kpiCards.nth(1).locator('.telemetry-kpi-label')).toHaveText('Calls');
    await expect(kpiCards.nth(2).locator('.telemetry-kpi-label')).toHaveText('Batch Intensity (Points/Call)');

    // 4. Date Limits with Info Archetype (UK standard DD/MM/YY format)
    const dateLimits = page.locator('.telemetry-date-limits-cluster');
    await expect(dateLimits).toBeVisible();
    await expect(dateLimits.locator('.info-label').nth(0)).toHaveText('Earliest Available');
    await expect(dateLimits.locator('.info-field-readonly').nth(0)).toHaveText(/^\d{2}\/\d{2}\/\d{2}$/);
    await expect(dateLimits.locator('.info-label').nth(1)).toHaveText('Current Date');
    await expect(dateLimits.locator('.info-field-readonly').nth(1)).toHaveText(/^\d{2}\/\d{2}\/\d{2}$/);

    const activePeriod = page.locator('#comp-effective-period-group');
    await expect(activePeriod).toBeVisible();
    await expect(activePeriod.locator('.telemetry-active-main-label')).toHaveText('Active Period');
    await expect(activePeriod.locator('.info-label').nth(0)).toHaveText('Start');
    await expect(activePeriod.locator('.info-field-readonly').nth(0)).toHaveText(/^\d{2}\/\d{2}\/\d{2}$/);
    await expect(activePeriod.locator('.info-label').nth(1)).toHaveText('End');
    await expect(activePeriod.locator('.info-field-readonly').nth(1)).toHaveText(/^\d{2}\/\d{2}\/\d{2}$/);
    await expect(activePeriod.locator('#comp-active-period-info-btn')).toBeVisible();
    await expect(activePeriod.locator('#comp-active-period-info-btn .info-tooltip-bubble')).toHaveText(/The active period reflects/);

    // 5. Sources and Metric Enclosed Multi-Select Segmented Trays
    const multiSelectTrays = page.locator('.telemetry-btn-toggle-group.multi-select');
    await expect(multiSelectTrays).toHaveCount(2);

    // Sources tray buttons have individual class with status dots
    const sourceBtns = multiSelectTrays.nth(0).locator('.telemetry-toggle-btn.individual');
    await expect(sourceBtns).toHaveCount(3);

    // Metric tray buttons
    const pointsMetricBtn = page.locator('#comp-metric-points-btn');
    const callsMetricBtn = page.locator('#comp-metric-calls-btn');
    await expect(pointsMetricBtn).toBeVisible();
    await expect(callsMetricBtn).toBeVisible();

    // 6. Stationary Total Toggle Button (Top-Left Anchor without icon)
    const totalBtn = page.locator('#comp-toggle-total');
    await expect(totalBtn).toBeVisible();
    await expect(totalBtn).toHaveText('Total');
    await expect(totalBtn).not.toHaveText(/✓/);

    // Measure button position before toggle
    const initialBox = await totalBtn.boundingBox();
    expect(initialBox).not.toBeNull();

    // Click to turn OFF total
    await totalBtn.click();
    await expect(totalBtn).not.toHaveClass(/active/);
    const afterOffBox = await totalBtn.boundingBox();
    expect(afterOffBox).not.toBeNull();
    // Bounding box x coordinate MUST be stationary (0 layout shift)
    expect(Math.abs((afterOffBox?.x ?? 0) - (initialBox?.x ?? 0))).toBeLessThanOrEqual(1);

    // Click to turn ON total again
    await totalBtn.click();
    await expect(totalBtn).toHaveClass(/active/);
    const afterOnBox = await totalBtn.boundingBox();
    expect(afterOnBox).not.toBeNull();
    expect(Math.abs((afterOnBox?.x ?? 0) - (initialBox?.x ?? 0))).toBeLessThanOrEqual(1);

    // 7. Timeline SVG Axis Titles
    const callsAxisTitle = page.locator('.chart-axis-title.calls-title');
    const pointsAxisTitle = page.locator('.chart-axis-title.points-title');
    await expect(callsAxisTitle).toBeVisible();
    await expect(callsAxisTitle).toHaveText('Calls');
    await expect(pointsAxisTitle).toBeVisible();
    await expect(pointsAxisTitle).toHaveText('Points');
  });
});
