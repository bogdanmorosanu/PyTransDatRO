import { expect, test } from './fixtures';
import { openApp } from './helpers/ui';

test.describe('Application Lifecycle & Global Controls', () => {
  test('brand logo fresh first-load reset completely purges inputs and floating panels while preserving language', async ({
    page,
  }) => {
    await openApp(page);

    // Change Direction, Coordinate Order, and Angular Unit
    await page.getByRole('radio', { name: 'Stereo70 to ETRS89' }).click();
    await page.locator('button[data-order-val="EN"]').click();
    await page.locator('#angular-unit-select').selectOption('dms');

    // Fill point coordinates
    await page.locator('#pt-input-1').fill('500000.000');
    await page.locator('#pt-input-2').fill('500000.000');

    // Open Help panel via Header dropdown
    await page.locator('#header-help-details summary').click();
    await page.locator('#menu-help-app-btn').click();
    await expect(page.locator('.floating-help-panel')).toBeVisible();

    // Click Brand Logo for complete reset
    await page.locator('#header-brand-logo-wrap').click();

    // Verify Direction reset to ETRS89 to Stereo70
    await expect(page.getByRole('radio', { name: 'ETRS89 to Stereo70' })).toBeChecked();

    // Verify Coordinate Order reset to NE(H)
    await expect(page.locator('button[data-order-val="NE"]')).toHaveClass(/active/);

    // Verify Angular Unit reset to degrees
    await expect(page.locator('#angular-unit-select')).toHaveValue('degrees');

    // Verify Point inputs are reset to fresh first-load defaults (46° N, 25° E)
    await expect(page.locator('#pt-input-1')).toHaveValue('46.00000000');
    await expect(page.locator('#pt-input-2')).toHaveValue('25.00000000');

    // Verify Help panel is closed
    await expect(page.locator('.floating-help-panel')).not.toBeVisible();
  });

  test('brand title soft navigation returns to studio without resetting active inputs', async ({
    page,
  }) => {
    await openApp(page);

    // Enter coordinates in Point tab
    await page.locator('#pt-input-1').fill('46.12345678');
    await page.locator('#pt-input-2').fill('25.98765432');

    // Switch to Telemetry view
    await page.locator('#nav-btn-telemetry').click();
    await expect(page.locator('.telemetry-page-view')).toBeVisible();

    // Click Brand Title Text for soft navigation back to studio
    await page.locator('.brand-title-text').click();

    // Studio workbench is visible and inputs are preserved
    await expect(page.locator('#pt-input-1')).toHaveValue('46.12345678');
    await expect(page.locator('#pt-input-2')).toHaveValue('25.98765432');
  });

  test('global offline alert banner and health dot react to API outage and open feedback modal', async ({
    page,
  }) => {
    // Intercept /health to simulate service outage
    await page.route('**/health', (route) => route.fulfill({ status: 503, body: 'Service Unavailable' }));
    await openApp(page);

    // Health dot in header reflects offline state
    const healthDot = page.locator('#header-health-dot');
    await expect(healthDot).toHaveClass(/offline/);

    // Global offline banner appears below header
    const offlineBanner = page.locator('#global-offline-banner');
    await expect(offlineBanner).toBeVisible();

    // Click "Report issue" button in the offline banner
    await page.locator('#global-offline-report-btn').click();

    // Feedback panel opens with report type pre-selected as api_down
    const feedbackPanel = page.locator('.floating-feedback-panel');
    await expect(feedbackPanel).toBeVisible();
    await expect(page.locator('#fb-type-select')).toHaveValue('api_down');
  });

  test('floating panels: grid metadata search filtering and geodetic help live benchmark card', async ({
    page,
  }) => {
    await openApp(page);

    // Open Grid Metadata from Header dropdown
    await page.locator('#header-help-details summary').click();
    await page.locator('#menu-grid-info-btn').click();

    const gridPanel = page.locator('.floating-grid-panel');
    await expect(gridPanel).toBeVisible();

    // Filter tree using search input
    await page.locator('#fgp-search').fill('etrs89');
    await expect(gridPanel.locator('.fgp-tree-branch').first()).toBeVisible();

    // Close Grid Metadata panel
    await page.locator('#fgp-close-btn').click();
    await expect(gridPanel).not.toBeVisible();

    // Open Geodetic Help panel
    await page.locator('#header-help-details summary').click();
    await page.locator('#menu-help-app-btn').click();

    const helpPanel = page.locator('.floating-help-panel');
    await expect(helpPanel).toBeVisible();

    // Live ANCPI benchmark card should be rendered
    const benchSlot = helpPanel.locator('#fhp-benchmark-slot');
    await expect(benchSlot).toBeVisible();
    await expect(benchSlot.locator('.fhp-benchmark-table')).toBeVisible();

    // Close help panel via header close button
    await page.locator('#fhp-close-btn').click();
    await expect(helpPanel).not.toBeVisible();
  });

  test('feedback modal bot challenge rejects incorrect math answers', async ({
    page,
  }) => {
    await openApp(page);

    // Open Feedback modal from Header dropdown
    await page.locator('#header-help-details summary').click();
    await page.locator('#menu-feedback-btn').click();

    const feedbackPanel = page.locator('.floating-feedback-panel');
    await expect(feedbackPanel).toBeVisible();

    // Fill description
    await page.locator('#fb-message-input').fill('Detailed discrepancy description for testing.');

    // Intentionally enter wrong math answer
    await page.locator('#fb-math-input').fill('99');

    // Click Submit
    await page.locator('#fb-submit-btn').click();

    // Warning toast is displayed and panel remains open
    await expect(page.locator('.toast.warning')).toBeVisible();
    await expect(feedbackPanel).toBeVisible();
  });
});
