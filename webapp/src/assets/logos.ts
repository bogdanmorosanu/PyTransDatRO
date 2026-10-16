export interface LogoPalette {
  bg: string;
  border: string;
  grid: string;
  muted: string;
  accent1: string;
  accent2: string;
  main: string;
}

export const CSS_VAR_PALETTE: LogoPalette = {
  bg: 'var(--logo-bg, #14171c)',
  border: 'var(--logo-border, #30363d)',
  grid: 'var(--logo-grid, #252a32)',
  muted: 'var(--logo-node-muted, #8b949e)',
  accent1: 'var(--logo-accent-1, #8b949e)',
  accent2: 'var(--logo-accent-2, #56b6c2)',
  main: 'var(--logo-node-main, #e6edf3)',
};

export const DARK_STANDALONE_PALETTE: LogoPalette = {
  bg: '#14171c',
  border: '#30363d',
  grid: '#252a32',
  muted: '#8b949e',
  accent1: '#8b949e',
  accent2: '#56b6c2',
  main: '#e6edf3',
};

export const LIGHT_STANDALONE_PALETTE: LogoPalette = {
  bg: '#ffffff',
  border: '#d0d7de',
  grid: '#eaeef2',
  muted: '#6e7781',
  accent1: '#57606a',
  accent2: '#0e7490',
  main: '#1f2328',
};

/**
 * Official PyTransDatRO Brand Logo (C2: Map Projection Graticule)
 * Curvilinear globe parallels and meridians morphing into a flat 2D digital cartographic grid.
 */
export function renderMapProjectionGraticuleInner(p: LogoPalette): string {
  return `
    <rect x="10" y="10" width="20" height="20" rx="3" stroke="${p.grid}" stroke-width="1.3" fill="none" />
    <path d="M10 16.5Q20 13.5 30 16.5" stroke="${p.accent1}" stroke-width="1.6" fill="none" />
    <path d="M10 23.5Q20 26.5 30 23.5" stroke="${p.accent1}" stroke-width="1.6" fill="none" />
    <path d="M16.5 10Q13.5 20 16.5 30" stroke="${p.accent2}" stroke-width="1.6" fill="none" />
    <path d="M23.5 10Q26.5 20 23.5 30" stroke="${p.accent2}" stroke-width="1.6" fill="none" />
    <circle cx="20" cy="20" r="2.6" fill="${p.main}" />
  `;
}

export function getLogoInlineSvg(size: number = 36): string {
  const p = CSS_VAR_PALETTE;
  return `<svg width="${size}" height="${size}" viewBox="0 0 40 40" fill="none" xmlns="http://www.w3.org/2000/svg" class="brand-logo-svg" aria-label="PyTransDatRO Logo">
    <rect x="0.75" y="0.75" width="38.5" height="38.5" rx="7.25" fill="${p.bg}" stroke="${p.border}" stroke-width="1.5"/>
    ${renderMapProjectionGraticuleInner(p)}
  </svg>`;
}

export function getLogoStandaloneSvg(theme: 'dark' | 'light'): string {
  const p = theme === 'light' ? LIGHT_STANDALONE_PALETTE : DARK_STANDALONE_PALETTE;
  return `<svg width="40" height="40" viewBox="0 0 40 40" fill="none" xmlns="http://www.w3.org/2000/svg">
  <rect x="0.75" y="0.75" width="38.5" height="38.5" rx="7.25" fill="${p.bg}" stroke="${p.border}" stroke-width="1.5"/>
${renderMapProjectionGraticuleInner(p).trim()}
</svg>
`;
}

export function syncFaviconWithTheme(theme: 'dark' | 'light'): void {
  if (typeof document === 'undefined') return;
  const svgMarkup = getLogoStandaloneSvg(theme);
  const encoded = `data:image/svg+xml;utf8,${encodeURIComponent(svgMarkup)}`;
  let link = document.querySelector("link[rel~='icon']") as HTMLLinkElement | null;
  if (!link) {
    link = document.createElement('link');
    link.rel = 'icon';
    link.type = 'image/svg+xml';
    document.head.appendChild(link);
  }
  link.href = encoded;
}
