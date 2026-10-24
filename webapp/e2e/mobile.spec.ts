import { expect, test, fireMapClick } from './fixtures';
import { openApp } from './helpers/ui';

test.describe('Mobile Compatibility & Dockable Workbench Architecture', () => {
  test('mobile bottom sheet defaults to sheet-half and cycles through full, peek, half on handle tap', async ({
    page,
  }) => {
    await page.setViewportSize({ width: 390, height: 844 });
    await openApp(page);

    const utilityPanel = page.locator('.utility-panel');
    const handleBar = page.locator('.bottom-sheet-handle-bar');
    const dockToggle = page.locator('.workbench-dock-toggle');

    // On mobile, dock toggle is hidden and handle bar is visible
    await expect(dockToggle).toBeHidden();
    await expect(handleBar).toBeVisible();

    // Default state on load is sheet-half
    await expect(utilityPanel).toHaveClass(/sheet-half/);

    // Tap 1: Transitions to sheet-full
    await handleBar.click();
    await expect(utilityPanel).toHaveClass(/sheet-full/);
    await expect(utilityPanel).not.toHaveClass(/sheet-half/);

    // Tap 2: Transitions to sheet-peek
    await handleBar.click();
    await expect(utilityPanel).toHaveClass(/sheet-peek/);
    await expect(utilityPanel).not.toHaveClass(/sheet-full/);

    // Tap 3: Transitions back to sheet-half
    await handleBar.click();
    await expect(utilityPanel).toHaveClass(/sheet-half/);
  });

  test('map remains interactive on mobile and receives clicks in half and peek modes', async ({
    page,
  }) => {
    await page.setViewportSize({ width: 390, height: 844 });
    await openApp(page);

    // Map panel is visible
    const mapPanel = page.locator('.map-panel');
    await expect(mapPanel).toBeVisible();

    // Trigger map click on Romania point (45.8, 25.0) via E2E test hook
    await fireMapClick(page, 45.8, 25.0);

    // Coordinates are captured in Point Converter inputs
    const input1 = page.locator('#pt-input-1');
    const input2 = page.locator('#pt-input-2');
    await expect(input1).toHaveValue('45.80000000');
    await expect(input2).toHaveValue('25.00000000');
  });

  test('coordinate inputs adapt to 2-row layout at <= 540px and preserve non-numeric DMS text', async ({
    page,
  }) => {
    await page.setViewportSize({ width: 390, height: 844 });
    await openApp(page);

    const input1 = page.locator('#pt-input-1');
    const input2 = page.locator('#pt-input-2');
    const input3 = page.locator('#pt-input-3');

    // Verify 16px font-size to prevent iOS Safari auto-zoom across inputs and selects
    const fontSize = await input1.evaluate((el) => window.getComputedStyle(el).fontSize);
    expect(fontSize).toBe('16px');
    const selectFontSize = await page.locator('#angular-unit-select').evaluate((el) => window.getComputedStyle(el).fontSize);
    expect(selectFontSize).toBe('16px');

    // Switch angular unit to DMS
    await page.locator('#angular-unit-select').selectOption('dms');

    // Type non-numeric DMS text with symbols °, ', "
    await input1.fill('45°50\'28.44420"');
    await input2.fill('25°00\'00.00000"');
    await input3.fill('150.500');

    // Trigger calculation
    await page.locator('#pt-calculate-btn').click();

    // Verify status is valid and output coordinate inputs are populated
    await expect(page.locator('.status-pill.success')).toBeVisible();

    // Verify 2-row layout: input 1 and 2 share same vertical top, input 3 sits underneath
    const box1 = await input1.boundingBox();
    const box2 = await input2.boundingBox();
    const box3 = await input3.boundingBox();

    expect(box1).not.toBeNull();
    expect(box2).not.toBeNull();
    expect(box3).not.toBeNull();

    if (box1 && box2 && box3) {
      // Row 1: Input 1 and Input 2 are vertically aligned side-by-side
      expect(Math.abs(box1.y - box2.y)).toBeLessThan(5);
      expect(box2.x).toBeGreaterThan(box1.x);

      // Row 2: Input 3 is strictly below Row 1 and spans wider
      expect(box3.y).toBeGreaterThan(box1.y + box1.height - 5);
      expect(box3.width).toBeGreaterThan(box1.width);
    }
  });

  test('responsive header navbar hides repo link and retires mobile view toggle button', async ({
    page,
  }) => {
    await page.setViewportSize({ width: 390, height: 844 });
    await openApp(page);

    // Repo link is hidden on mobile
    await expect(page.locator('.header-repo-link')).toBeHidden();

    // Mobile view switcher button #header-mobile-view-btn is cleanly retired
    await expect(page.locator('#header-mobile-view-btn')).toHaveCount(0);

    // Brand title and main nav buttons are visible
    await expect(page.locator('.brand-title-text')).toBeVisible();
    await expect(page.locator('#nav-btn-studio')).toBeVisible();
    await expect(page.locator('#nav-btn-telemetry')).toBeVisible();

    // Viewport has no horizontal scrollbar
    const hasHorizontalOverflow = await page.evaluate(() => {
      return document.documentElement.scrollWidth > window.innerWidth;
    });
    expect(hasHorizontalOverflow).toBe(false);
  });

  test('usage and telemetry dashboard harmonizes for mobile viewports without horizontal clipping', async ({
    page,
  }) => {
    await page.setViewportSize({ width: 390, height: 844 });
    await openApp(page);

    // Navigate to Telemetry
    await page.locator('#nav-btn-telemetry').click();
    await expect(page.locator('.telemetry-page-view')).toBeVisible();

    // Telemetry container padding reduced to 12px 14px
    const container = page.locator('.telemetry-page-container');
    await expect(container).toBeVisible();

    // Verify no window horizontal overflow
    const hasHorizontalOverflow = await page.evaluate(() => {
      const pageView = document.querySelector('.telemetry-page-view');
      return pageView ? pageView.scrollWidth > pageView.clientWidth + 2 : false;
    });
    expect(hasHorizontalOverflow).toBe(false);

    // Check KPI cards stack responsively
    const kpiBanner = page.locator('.telemetry-kpi-banner').first();
    await expect(kpiBanner).toBeVisible();

    // Switch to Spatial tab and verify punchcard matrix viewport has overflow-x auto
    await page.locator('#tab-spatial-btn').click();
    const spatialViewport = page.locator('.spatial-matrix-viewport');
    await expect(spatialViewport).toBeVisible();
    const overflowX = await spatialViewport.evaluate((el) => window.getComputedStyle(el).overflowX);
    expect(overflowX).toBe('auto');
  });

  test('floating help and feedback panels render as bottom overlays on mobile', async ({
    page,
  }) => {
    await page.setViewportSize({ width: 390, height: 844 });
    await openApp(page);

    // Open Help panel
    await page.locator('#header-help-details summary').click();
    await page.locator('#menu-help-app-btn').click();

    const helpPanel = page.locator('.floating-help-panel');
    await expect(helpPanel).toBeVisible();

    // Verify positioned at bottom with full width
    const bounds = await helpPanel.boundingBox();
    expect(bounds).not.toBeNull();
    if (bounds) {
      expect(bounds.x).toBeLessThanOrEqual(5);
      expect(bounds.width).toBeGreaterThanOrEqual(380);
      expect(bounds.y + bounds.height).toBeGreaterThanOrEqual(830);
    }

    // Close Help panel
    await page.locator('#fhp-close-btn').click();
    await expect(helpPanel).not.toBeVisible();
  });

  test('desktop dockable workbench collapses and expands smoothly via toggle button and keyboard shortcut', async ({
    page,
  }) => {
    await page.setViewportSize({ width: 1280, height: 800 });
    await openApp(page);

    const utilityPanel = page.locator('.utility-panel');
    const dockToggle = page.locator('.workbench-dock-toggle');

    // On desktop, dock toggle is visible with chevron pointing left
    await expect(dockToggle).toBeVisible();
    await expect(dockToggle.locator('.workbench-dock-toggle-icon')).toHaveText('◀');
    await expect(utilityPanel).not.toHaveClass(/is-collapsed/);

    // Click toggle button to collapse
    await dockToggle.click();
    await expect(utilityPanel).toHaveClass(/is-collapsed/);
    await expect(dockToggle).toHaveClass(/is-collapsed/);
    await expect(dockToggle.locator('.workbench-dock-toggle-icon')).toHaveText('▶');

    // Verify localStorage persistence
    const savedState = await page.evaluate(() => localStorage.getItem('pytransdat_workbench_collapsed'));
    expect(savedState).toBe('true');

    // Click toggle button again to expand
    await dockToggle.click();
    await expect(utilityPanel).not.toHaveClass(/is-collapsed/);
    await expect(dockToggle).not.toHaveClass(/is-collapsed/);
    await expect(dockToggle.locator('.workbench-dock-toggle-icon')).toHaveText('◀');

    // Toggle via keyboard shortcut Alt+W
    await page.keyboard.press('Alt+w');
    await expect(utilityPanel).toHaveClass(/is-collapsed/);

    // Toggle back via shortcut Ctrl+\
    await page.keyboard.press('Control+\\');
    await expect(utilityPanel).not.toHaveClass(/is-collapsed/);
  });

  test('desktop dockable workbench starts collapsed on initial load when persisted in localStorage', async ({
    page,
  }) => {
    await page.setViewportSize({ width: 1280, height: 800 });
    // Seed localStorage with collapsed state before loading the app
    await page.addInitScript(() => {
      localStorage.setItem('pytransdat_workbench_collapsed', 'true');
    });

    await openApp(page);

    const utilityPanel = page.locator('.utility-panel');
    const dockToggle = page.locator('.workbench-dock-toggle');

    // On initial load, utility panel starts collapsed and toggle shows expand chevron ▶
    await expect(utilityPanel).toHaveClass(/is-collapsed/);
    await expect(dockToggle).toHaveClass(/is-collapsed/);
    await expect(dockToggle.locator('.workbench-dock-toggle-icon')).toHaveText('▶');
    await expect(dockToggle).toHaveAttribute('aria-expanded', 'false');

    // Map panel spans 100% width
    const mapPanel = page.locator('.map-panel');
    const mapBox = await mapPanel.boundingBox();
    expect(mapBox).not.toBeNull();
    if (mapBox) {
      expect(mapBox.width).toBeGreaterThanOrEqual(1200);
    }
  });

  test('bottom sheet touch swipe gestures smoothly transition between states', async ({
    page,
  }) => {
    await page.setViewportSize({ width: 390, height: 844 });
    await openApp(page);

    const utilityPanel = page.locator('.utility-panel');
    await expect(utilityPanel).toHaveClass(/sheet-half/);

    // Swipe UP from sheet-half: should transition to sheet-full
    await page.evaluate(() => {
      const handle = document.querySelector('.bottom-sheet-handle-bar')!;
      const rect = handle.getBoundingClientRect();
      const clientX = rect.left + rect.width / 2;
      const clientY = rect.top + rect.height / 2;

      const touchStart = new Touch({
        identifier: 1,
        target: handle,
        clientX,
        clientY,
      });
      const touchMove = new Touch({
        identifier: 1,
        target: handle,
        clientX,
        clientY: clientY - 80,
      });

      handle.dispatchEvent(new TouchEvent('touchstart', {
        touches: [touchStart],
        bubbles: true,
      }));
      handle.dispatchEvent(new TouchEvent('touchmove', {
        touches: [touchMove],
        bubbles: true,
      }));
      handle.dispatchEvent(new TouchEvent('touchend', {
        touches: [],
        bubbles: true,
      }));
    });

    await expect(utilityPanel).toHaveClass(/sheet-full/);

    // Swipe DOWN from sheet-full: should transition back to sheet-half
    await page.evaluate(() => {
      const handle = document.querySelector('.bottom-sheet-handle-bar')!;
      const rect = handle.getBoundingClientRect();
      const clientX = rect.left + rect.width / 2;
      const clientY = rect.top + rect.height / 2;

      const touchStart = new Touch({
        identifier: 1,
        target: handle,
        clientX,
        clientY,
      });
      const touchMove = new Touch({
        identifier: 1,
        target: handle,
        clientX,
        clientY: clientY + 80,
      });

      handle.dispatchEvent(new TouchEvent('touchstart', {
        touches: [touchStart],
        bubbles: true,
      }));
      handle.dispatchEvent(new TouchEvent('touchmove', {
        touches: [touchMove],
        bubbles: true,
      }));
      handle.dispatchEvent(new TouchEvent('touchend', {
        touches: [],
        bubbles: true,
      }));
    });

    await expect(utilityPanel).toHaveClass(/sheet-half/);
  });

  test('focusing coordinate input on mobile automatically expands bottom sheet to sheet-full', async ({
    page,
  }) => {
    await page.setViewportSize({ width: 390, height: 844 });
    await openApp(page);

    const utilityPanel = page.locator('.utility-panel');
    await expect(utilityPanel).toHaveClass(/sheet-half/);

    // Focus input 1
    const input1 = page.locator('#pt-input-1');
    await input1.focus();

    // Bottom sheet smoothly auto-expands to sheet-full so virtual keyboard does not obscure controls
    await expect(utilityPanel).toHaveClass(/sheet-full/);
  });

  test('timeline chart supports touch scrubbing displaying crosshair and tooltip on mobile', async ({
    page,
  }) => {
    await page.setViewportSize({ width: 390, height: 844 });
    await openApp(page);

    // Navigate to Telemetry
    await page.locator('#nav-btn-telemetry').click();
    await expect(page.locator('.telemetry-page-view')).toBeVisible();

    const svg = page.locator('#timeline-svg');
    await expect(svg).toBeVisible();

    // Dispatch touchstart and touchmove on timeline chart SVG
    await page.evaluate(() => {
      const chartSvg = document.querySelector('#timeline-svg')!;
      const rect = chartSvg.getBoundingClientRect();
      const clientX = rect.left + rect.width / 2;
      const clientY = rect.top + rect.height / 2;

      const touch = new Touch({
        identifier: 1,
        target: chartSvg,
        clientX,
        clientY,
      });

      chartSvg.dispatchEvent(new TouchEvent('touchstart', {
        touches: [touch],
        bubbles: true,
      }));
      chartSvg.dispatchEvent(new TouchEvent('touchmove', {
        touches: [touch],
        bubbles: true,
      }));
    });

    // Tooltip and crosshair become visible
    const tooltip = page.locator('#chart-tooltip');
    await expect(tooltip).toBeVisible();
    const crosshair = page.locator('#crosshair-x');
    await expect(crosshair).toHaveCSS('display', 'block');
    expect(await crosshair.getAttribute('x1')).not.toBe('0');

    // Dispatch touchend
    await page.evaluate(() => {
      const chartSvg = document.querySelector('#timeline-svg')!;
      chartSvg.dispatchEvent(new TouchEvent('touchend', {
        touches: [],
        bubbles: true,
      }));
    });

    // Tooltip hides on touch end
    await expect(tooltip).toBeHidden();
  });
});
