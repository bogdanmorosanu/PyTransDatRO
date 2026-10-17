import { defineConfig, devices } from '@playwright/test';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

import fs from 'node:fs';

const webappDir = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(webappDir, '..');
const tmpDir = path.join(webappDir, 'e2e', '.tmp');

const E2E_HOST = '127.0.0.1';
const E2E_PORT = 8100;
const BASE_URL = `http://${E2E_HOST}:${E2E_PORT}`;
const UNLIMITED_RATE = '100000';

// Auto-detect project Python: process.env.PYTHON -> standard Conda env on Windows -> .venv -> system python.
function resolvePython(): string {
  if (process.env.PYTHON) return process.env.PYTHON;
  const defaultCondaWindows = 'D:\\anaconda3\\envs\\pytransdat\\python.exe';
  if (fs.existsSync(defaultCondaWindows)) return defaultCondaWindows;
  const localVenvPosix = path.join(repoRoot, '.venv', 'bin', 'python');
  if (fs.existsSync(localVenvPosix)) return localVenvPosix;
  const localVenvWindows = path.join(repoRoot, '.venv', 'Scripts', 'python.exe');
  if (fs.existsSync(localVenvWindows)) return localVenvWindows;
  return 'python';
}

// Quoted because Windows conda paths often contain spaces.
const PYTHON = `"${resolvePython()}"`;

export default defineConfig({
  testDir: './e2e',
  retries: 0,
  reporter: [['list'], ['html', { open: 'never' }]],
  use: {
    baseURL: BASE_URL,
    locale: 'en-US',
    trace: 'retain-on-failure',
    screenshot: 'only-on-failure',
  },
  projects: [{ name: 'chromium', use: { ...devices['Desktop Chrome'] } }],
  webServer: {
    // Build the frontend with the e2e map handle enabled, then serve web + API from one process.
    command:
      'npm --prefix webapp run build:e2e && ' +
      `${PYTHON} -m uvicorn api.main:app --host ${E2E_HOST} --port ${E2E_PORT}`,
    cwd: repoRoot,
    url: `${BASE_URL}/health`,
    reuseExistingServer: false,
    timeout: 300_000,
    env: {
      VITE_E2E: 'true',
      PYTRANSDAT_WEBAPP_DIST_DIR: path.join(webappDir, 'dist-e2e'),
      PYTRANSDAT_LOG_DIR: path.join(tmpDir, 'logs'),
      PYTRANSDAT_FEEDBACK_LOG_PATH: path.join(tmpDir, 'feedback.log'),
      PYTRANSDAT_SMTP_HOST: '',
      PYTRANSDAT_RATE_LIMIT_POINT: UNLIMITED_RATE,
      PYTRANSDAT_RATE_LIMIT_BATCH: UNLIMITED_RATE,
    },
  },
});
