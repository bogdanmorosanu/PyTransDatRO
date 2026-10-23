import type { Page, Download } from '@playwright/test';
import { readFile } from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const helpersDir = path.dirname(fileURLToPath(import.meta.url));
export const FIXTURES_FILES_DIR = path.resolve(helpersDir, '..', 'fixtures', 'files');

export interface MemoryFile {
  name: string;
  content: string;
  mimeType?: string;
}

/**
 * Ingests files into the FileUpload component from in-memory string buffers.
 */
export async function uploadMemoryFiles(page: Page, files: MemoryFile[]): Promise<void> {
  const fileInputs = files.map((f) => ({
    name: f.name,
    mimeType: f.mimeType || 'text/csv',
    buffer: Buffer.from(f.content, 'utf8'),
  }));
  await page.getByLabel('Drag & drop files here or click to browse').setInputFiles(fileInputs);
}

/**
 * Ingests files from the fixtures/files directory on disk.
 */
export async function uploadDiskFixtures(page: Page, fileNames: string[]): Promise<void> {
  const fullPaths = fileNames.map((name) => path.join(FIXTURES_FILES_DIR, name));
  await page.getByLabel('Drag & drop files here or click to browse').setInputFiles(fullPaths);
}

/**
 * Reads downloaded file content as a UTF-8 string.
 */
export async function readDownloadContent(download: Download): Promise<string> {
  const downloadPath = await download.path();
  if (!downloadPath) {
    throw new Error('Download path is null');
  }
  return readFile(downloadPath, 'utf8');
}
