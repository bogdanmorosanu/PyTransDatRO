import type { Page } from '@playwright/test';
import { expect, test } from './fixtures';

const POINT_ENDPOINT = '**/api/v1/transform/point';
const EMPTY_OUTPUT = '';

async function loadAndSettle(page: Page): Promise<void> {
  await page.goto('/');
  const results = page.locator('.result-card');
  await expect(results.getByLabel('Northing')).not.toHaveValue(EMPTY_OUTPUT);
}

test.describe('point transform API failures', () => {
  test('HTTP 500 shows an error to the user and no result', async ({ page, externalRequests }) => {
    void externalRequests;
    await loadAndSettle(page);
    await page.route(POINT_ENDPOINT, (route) =>
      route.fulfill({ status: 500, contentType: 'application/json', body: '{"detail":"boom"}' }),
    );

    await page.getByRole('button', { name: 'Transform', exact: true }).click();

    await expect(page.getByText(/API Error: 500/)).toBeVisible();
    const results = page.locator('.result-card');
    await expect(results.getByLabel('Northing')).toHaveValue(EMPTY_OUTPUT);
    await expect(results.getByLabel('Easting')).toHaveValue(EMPTY_OUTPUT);
  });

  test('HTTP 429 shows an error / rate-limit message to the user and no result', async ({
    page,
    externalRequests,
  }) => {
    void externalRequests;
    await loadAndSettle(page);
    await page.route(POINT_ENDPOINT, (route) =>
      route.fulfill({
        status: 429,
        contentType: 'application/json',
        headers: { 'Retry-After': '30' },
        body: '{"detail":"Rate limit exceeded"}',
      }),
    );

    await page.getByRole('button', { name: 'Transform', exact: true }).click();

    await expect(page.locator('.toast').getByText(/429|rate limit|too many requests/i)).toBeVisible();
    const results = page.locator('.result-card');
    await expect(results.getByLabel('Northing')).toHaveValue(EMPTY_OUTPUT);
  });
});
