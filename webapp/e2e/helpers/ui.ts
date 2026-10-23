import type { Page } from '@playwright/test';
import { expect } from '../fixtures';

export type TabName = 'Point' | 'Text' | 'Files';

export async function openApp(page: Page): Promise<void> {
  await page.goto('/');
  await expect(page.locator('#map .leaflet-tile-pane')).toBeAttached();
}

export async function openTab(page: Page, tab: TabName): Promise<void> {
  await page.locator('.tabs-bar').getByRole('button', { name: tab, exact: true }).click();
}

export async function setDirectionStereoToEtrs(page: Page): Promise<void> {
  await page.getByRole('radio', { name: 'Stereo70 to ETRS89' }).click();
  await expect(page.getByRole('radio', { name: 'Stereo70 to ETRS89' })).toBeChecked();
}

export async function setDirectionEtrsToStereo(page: Page): Promise<void> {
  await page.getByRole('radio', { name: 'ETRS89 to Stereo70' }).click();
  await expect(page.getByRole('radio', { name: 'ETRS89 to Stereo70' })).toBeChecked();
}

export async function setCoordinateOrder(page: Page, order: 'NE' | 'EN'): Promise<void> {
  const btn = page.locator(`.sliding-toggle-btn[data-order-val="${order}"]`);
  await btn.click();
  await expect(btn).toHaveClass(/active/);
}

export async function setAngularUnit(page: Page, unit: 'degrees' | 'dms' | 'radians'): Promise<void> {
  await page.locator('#angular-unit-select').selectOption(unit);
  await expect(page.locator('#angular-unit-select')).toHaveValue(unit);
}

/** Splits a textarea value into lines (1-to-1, blanks preserved). */
export function lines(text: string): string[] {
  return text.split(/\r?\n/);
}

/** Parses "a, b, c" style output into numbers; returns [] for non-numeric tokens. */
export function numbers(line: string): number[] {
  return line
    .split(',')
    .map((s) => s.trim())
    .filter((s) => s.length > 0)
    .map(Number);
}
