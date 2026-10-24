import { expect, test } from './fixtures';

test.describe('Usage Spatial Subpage', () => {
  test('renders 2 centered Archetype 8.1 KPI cards, sentence-case quadrants, info date limits, multi-select sources, and mutually exclusive metric tray', async ({ page }) => {
    await page.goto('/');

    // 1. Navigate to Usage page
    const usageNavBtn = page.locator('#nav-btn-telemetry');
    await usageNavBtn.click();
    await expect(page.locator('.telemetry-page-container')).toBeVisible();

    // 2. Switch to Spatial tab
    const spatialTab = page.locator('#tab-spatial-btn');
    await expect(spatialTab).toBeVisible();
    await spatialTab.click();
    await expect(spatialTab).toHaveClass(/active/);

    // 3. Top KPI Area: 2 Centered Archetype 8.1 KPI Cards (Points/Calls & Active Cells)
    const kpiBanner = page.locator('.telemetry-kpi-banner.two-cards');
    await expect(kpiBanner).toBeVisible();
    const kpiCards = kpiBanner.locator('.telemetry-kpi-card');
    await expect(kpiCards).toHaveCount(2);

    // Initial default metric is Points
    await expect(kpiCards.nth(0).locator('.telemetry-kpi-label')).toHaveText('Points');
    await expect(kpiCards.nth(0).locator('.telemetry-kpi-val')).toBeVisible();
    await expect(kpiCards.nth(1).locator('.telemetry-kpi-label')).toHaveText('Active Cells');
    await expect(kpiCards.nth(1).locator('.telemetry-kpi-val')).toBeVisible();

    // 4. Mathematical Quadrants Deck: 4 Cards with Natural Sentence-Case Names
    const quadrantDeck = page.locator('.spatial-quadrants-deck');
    await expect(quadrantDeck).toBeVisible();
    const quadrantCards = quadrantDeck.locator('.quadrant-card');
    await expect(quadrantCards).toHaveCount(4);

    await expect(quadrantCards.nth(0).locator('.quadrant-name')).toHaveText('North-West (NW)');
    await expect(quadrantCards.nth(1).locator('.quadrant-name')).toHaveText('North-East (NE)');
    await expect(quadrantCards.nth(2).locator('.quadrant-name')).toHaveText('South-West (SW)');
    await expect(quadrantCards.nth(3).locator('.quadrant-name')).toHaveText('South-East (SE)');

    // Verify text-transform on quadrant-name is none
    const nwTransform = await quadrantCards.nth(0).locator('.quadrant-name').evaluate((el) => {
      return window.getComputedStyle(el).textTransform;
    });
    expect(nwTransform).toBe('none');

    // 5. Date Limits with Info Archetype (UK standard DD/MM/YY format)
    const dateLimits = page.locator('.telemetry-date-limits-cluster');
    await expect(dateLimits).toBeVisible();
    await expect(dateLimits.locator('.info-label').nth(0)).toHaveText('Earliest Available');
    await expect(dateLimits.locator('.info-field-readonly').nth(0)).toHaveText(/^\d{2}\/\d{2}\/\d{2}$/);
    await expect(dateLimits.locator('.info-label').nth(1)).toHaveText('Current Date');
    await expect(dateLimits.locator('.info-field-readonly').nth(1)).toHaveText(/^\d{2}\/\d{2}\/\d{2}$/);

    const activePeriod = page.locator('#spatial-effective-period-group');
    await expect(activePeriod).toBeVisible();
    await expect(activePeriod.locator('.telemetry-active-main-label')).toHaveText('Active Period');
    await expect(activePeriod.locator('.info-label').nth(0)).toHaveText('Start');
    await expect(activePeriod.locator('.info-field-readonly').nth(0)).toHaveText(/^\d{2}\/\d{2}\/\d{2}$/);
    await expect(activePeriod.locator('.info-label').nth(1)).toHaveText('End');
    await expect(activePeriod.locator('.info-field-readonly').nth(1)).toHaveText(/^\d{2}\/\d{2}\/\d{2}$/);
    await expect(activePeriod.locator('#spatial-active-period-info-btn')).toBeVisible();
    await expect(activePeriod.locator('#spatial-active-period-info-btn .info-tooltip-bubble')).toHaveText(/The active period reflects/);

    // 6. Sources Enclosed Multi-Select Tray (Dot Indicators)
    const sourcesTray = page.locator('.filters-left-group .telemetry-btn-toggle-group.multi-select');
    await expect(sourcesTray).toBeVisible();
    const sourceBtns = sourcesTray.locator('.telemetry-toggle-btn.individual');
    await expect(sourceBtns).toHaveCount(3);

    // 7. Metric Mutually Exclusive Segmented Pill Tray (Archetype 2.2)
    const metricTray = page.locator('.filters-right-group .telemetry-btn-toggle-group');
    await expect(metricTray).toBeVisible();
    await expect(metricTray).not.toHaveClass(/multi-select/);

    const pointsMetricBtn = page.locator('#spatial-metric-points-btn');
    const callsMetricBtn = page.locator('#spatial-metric-calls-btn');
    await expect(pointsMetricBtn).toBeVisible();
    await expect(callsMetricBtn).toBeVisible();
    await expect(pointsMetricBtn).toHaveClass(/active/);
    await expect(callsMetricBtn).not.toHaveClass(/active/);

    // Switch metric to Calls: Mutually exclusive transition
    await callsMetricBtn.click();
    await expect(callsMetricBtn).toHaveClass(/active/);
    await expect(pointsMetricBtn).not.toHaveClass(/active/);

    // Card 1 label dynamically updates to Calls
    await expect(kpiCards.nth(0).locator('.telemetry-kpi-label')).toHaveText('Calls');

    // Switch back to Points
    await pointsMetricBtn.click();
    await expect(pointsMetricBtn).toHaveClass(/active/);
    await expect(callsMetricBtn).not.toHaveClass(/active/);
    await expect(kpiCards.nth(0).locator('.telemetry-kpi-label')).toHaveText('Points');

    // 8. Heatmap Controls, Canvas Viewport, and Scaling Selector
    const canvasSlot = page.locator('#spatial-matrix-slot canvas');
    await expect(canvasSlot).toBeVisible();

    const scaleLogBtn = page.locator('#spatial-scale-log-btn');
    const scaleLinearBtn = page.locator('#spatial-scale-linear-btn');
    await expect(scaleLogBtn).toBeVisible();
    await expect(scaleLinearBtn).toBeVisible();
    await expect(scaleLogBtn).toHaveClass(/active/);

    await scaleLinearBtn.click();
    await expect(scaleLinearBtn).toHaveClass(/active/);
    await expect(scaleLogBtn).not.toHaveClass(/active/);

    await scaleLogBtn.click();
    await expect(scaleLogBtn).toHaveClass(/active/);
    await expect(scaleLinearBtn).not.toHaveClass(/active/);
  });
});
