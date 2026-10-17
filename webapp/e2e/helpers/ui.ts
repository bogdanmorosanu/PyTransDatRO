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
  await page.getByRole('radio', { name: 'Stereo 70 → ETRS89' }).click();
  await expect(page.getByRole('radio', { name: 'Stereo 70 → ETRS89' })).toBeChecked();
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
