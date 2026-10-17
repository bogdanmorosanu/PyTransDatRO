import { DISPLAY_DECIMALS, REFERENCE_POINTS, matchesAtDisplayPrecision } from './referencePoints';
import { expect, test } from './fixtures';
import { lines, numbers, openApp, openTab, setDirectionStereoToEtrs } from './helpers/ui';
import { readFile } from 'node:fs/promises';

const [P1, P2, P3] = REFERENCE_POINTS;

function st70Line(p: (typeof REFERENCE_POINTS)[number]): string {
  return `${p.north}, ${p.east}, ${p.height}`;
}

async function runBatch(page: import('@playwright/test').Page, input: string): Promise<string> {
  await openApp(page);
  await openTab(page, 'Text');
  await setDirectionStereoToEtrs(page);
  await page.getByRole('textbox', { name: 'Input Coordinates' }).fill(input);
  await page.getByRole('button', { name: 'Transform', exact: true }).click();
  const output = page.getByRole('textbox', { name: 'Computed Results' });
  await expect(output).toBeVisible();
  return output.inputValue();
}

test('batch keeps strict 1-to-1 lines and 1-based Point N labels around blank and comment lines', async ({
  page,
  externalRequests,
}) => {
  void externalRequests;
  const input = [
    '# header comment',
    st70Line(P1),
    '',
    '# middle comment',
    st70Line(P2),
    '',
    '',
    st70Line(P3),
  ].join('\n');

  const out = lines(await runBatch(page, input));

  // Strict 1-to-1: line N in == line N out, blanks and comments keep their position.
  expect(out).toHaveLength(8);
  expect(out[0]).toBe('# header comment');
  expect(out[2]).toBe('');
  expect(out[3]).toBe('# middle comment');
  expect(out[5]).toBe('');
  expect(out[6]).toBe('');

  const dataLines: Array<[number, (typeof REFERENCE_POINTS)[number]]> = [
    [1, P1],
    [4, P2],
    [7, P3],
  ];
  for (const [idx, p] of dataLines) {
    const [lat, lon, h] = numbers(out[idx]);
    expect(matchesAtDisplayPrecision(lat, p.lat, DISPLAY_DECIMALS.degrees)).toBe(true);
    expect(matchesAtDisplayPrecision(lon, p.lon, DISPLAY_DECIMALS.degrees)).toBe(true);
    expect(matchesAtDisplayPrecision(h, p.ellipsoidalHeight, DISPLAY_DECIMALS.metres)).toBe(true);
  }

  // Unnamed data points are labelled Point 1..3 (blank/comment lines do not consume an index).
  const [download] = await Promise.all([
    page.waitForEvent('download'),
    page.getByRole('button', { name: 'CSV', exact: true }).click(),
  ]);
  const csv = (await readFile((await download.path())!, 'utf8')).split(/\r?\n/);
  expect(csv).toHaveLength(3);
  csv.forEach((row, i) => {
    expect(row.startsWith(`Point ${i + 1},`)).toBe(true);
  });
});

test('batch reports a no-data-zone point as untransformable and not as a success', async ({
  page,
  externalRequests,
  request,
}) => {
  void externalRequests;
  // Stereo 70 point from tests/test_trans_ro.py (test_st70_to_etrs89_just_out_ro_border_no_data).
  const noData = '549364, 147553, 100';

  // Confirm via the real API that this input yields the "No data on grid" warning.
  const api = await request.post('/api/v1/transform/batch', {
    data: { op: 'Stereo70ToETRS89', points: [[549364, 147553, 100]], unit: 'degrees' },
  });
  expect(api.ok()).toBe(true);
  expect((await api.json()).results[0].warning).toBe('No data on grid');

  const out = lines(await runBatch(page, [st70Line(P1), noData].join('\n')));
  expect(out).toHaveLength(2);

  // The valid point is transformed.
  expect(matchesAtDisplayPrecision(numbers(out[0])[0], P1.lat, DISPLAY_DECIMALS.degrees)).toBe(true);

  // The no-data point is flagged and its input is not presented as transformed output.
  expect(out[1].toLowerCase()).toContain('no data');
  expect(out[1]).not.toContain('549364');

  // Summary pills: only the valid point counts as in-grid / successful.
  await expect(page.getByText(/In Grid: 1$/)).toBeVisible();
  await expect(page.getByText(/In Grid: 2/)).toHaveCount(0);
  await expect(page.getByText(/Outside Grid: 1|No Grid Data: 1/)).toBeVisible();
});
