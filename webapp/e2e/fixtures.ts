import { test as base, expect, type Page, type Route } from '@playwright/test';
import type { BoundsBox, MapViewSnapshot } from '../src/components/mapViewSnapshot';

/** localStorage key holding the UI language (see src/i18n/index.ts). */
export const LANGUAGE_STORAGE_KEY = 'pytransdat_lang';
/** Response header set on every stubbed external request; lets specs prove nothing leaked. */
export const STUB_HEADER = 'x-e2e-stub';

const LOCAL_HOSTS = new Set(['127.0.0.1', 'localhost']);

// 1x1 fully transparent PNG.
const TRANSPARENT_PNG = Buffer.from(
  'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==',
  'base64',
);

/** Plain-data snapshot returned by `window.__e2e.mapView()`; shared with the app, not copied. */
export type { BoundsBox, MapViewSnapshot };

/** Reads the e2e map handle. Throws if the app was not built with VITE_E2E=true. */
export async function readMapView(page: Page): Promise<MapViewSnapshot> {
  return page.evaluate(() => {
    const handle = (window as unknown as { __e2e?: { mapView(): unknown } }).__e2e;
    if (!handle) throw new Error('window.__e2e is missing (build without VITE_E2E=true?)');
    return handle.mapView();
  }) as Promise<MapViewSnapshot>;
}

export async function fireMapClick(page: Page, lat: number, lng: number): Promise<void> {
  await page.evaluate(
    ({ lat, lng }) => {
      const handle = (window as unknown as { __e2e?: { fireMapClick?(lat: number, lng: number): void } }).__e2e;
      if (!handle?.fireMapClick) throw new Error('fireMapClick missing on window.__e2e');
      handle.fireMapClick(lat, lng);
    },
    { lat, lng },
  );
}

async function stubExternal(route: Route): Promise<void> {
  const url = new URL(route.request().url());
  const type = route.request().resourceType();
  const headers = { [STUB_HEADER]: '1' };
  if (/\.(png|jpe?g)$/.test(url.pathname) || type === 'image') {
    await route.fulfill({ status: 200, contentType: 'image/png', body: TRANSPARENT_PNG, headers });
  } else if (type === 'stylesheet') {
    await route.fulfill({ status: 200, contentType: 'text/css', body: '', headers });
  } else {
    await route.fulfill({ status: 204, body: '', headers });
  }
}

interface Fixtures {
  /** URLs of external (non-local) requests that were intercepted and stubbed. */
  externalRequests: string[];
  /** Messages of console errors and uncaught page errors seen so far. */
  consoleErrors: string[];
  /** Reads `window.__e2e.mapView()` as plain data. */
  mapView: () => Promise<MapViewSnapshot>;
}

export const test = base.extend<Fixtures>({
  context: async ({ context }, use) => {
    // Pin language to English unless a test already chose one (survives reloads).
    await context.addInitScript((key: string) => {
      try {
        if (!window.localStorage.getItem(key)) window.localStorage.setItem(key, 'en');
      } catch (err) {
        console.warn('e2e: could not pin language', err);
      }
    }, LANGUAGE_STORAGE_KEY);
    await use(context);
  },
  externalRequests: async ({ context }, use) => {
    const seen: string[] = [];
    await context.route(
      (url) => !LOCAL_HOSTS.has(url.hostname) && url.protocol.startsWith('http'),
      async (route) => {
        seen.push(route.request().url());
        await stubExternal(route);
      },
    );
    await use(seen);
  },
  consoleErrors: async ({ page }, use) => {
    const errors: string[] = [];
    page.on('console', (msg) => {
      if (msg.type() === 'error') errors.push(msg.text());
    });
    page.on('pageerror', (err) => errors.push(err.message));
    await use(errors);
  },
  mapView: async ({ page, externalRequests }, use) => {
    void externalRequests; // ensure external stubbing is active whenever the map is used
    await use(() => readMapView(page));
  },
});

export { expect };
