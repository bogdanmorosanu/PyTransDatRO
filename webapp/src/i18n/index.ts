import { ro } from './ro';
import { en } from './en';

export type Language = 'ro' | 'en';

const dictionaries = { ro, en };

let currentLang: Language = (localStorage.getItem('pytransdat_lang') as Language) || 'ro';
const listeners: Array<(lang: Language) => void> = [];

export function getLanguage(): Language {
  return currentLang;
}

export function setLanguage(lang: Language): void {
  if (lang !== currentLang) {
    currentLang = lang;
    localStorage.setItem('pytransdat_lang', lang);
    listeners.forEach((fn) => fn(lang));
  }
}

export function toggleLanguage(): Language {
  const next = currentLang === 'ro' ? 'en' : 'ro';
  setLanguage(next);
  return next;
}

export function onLanguageChange(fn: (lang: Language) => void): () => void {
  listeners.push(fn);
  return () => {
    const idx = listeners.indexOf(fn);
    if (idx !== -1) listeners.splice(idx, 1);
  };
}

/**
 * Access a translation key like 'point.title' or 'app.title'.
 * Supports optional parameter interpolation: t('volumeWarning.message', { count: 5000 })
 */
export function t(path: string, params?: Record<string, string | number>): string {
  const dict = dictionaries[currentLang] || dictionaries.ro;
  const parts = path.split('.');
  let curr: any = dict;

  for (const part of parts) {
    if (curr && typeof curr === 'object' && part in curr) {
      curr = curr[part];
    } else {
      return path; // Fallback to key path if missing
    }
  }

  if (typeof curr !== 'string') return path;

  if (params) {
    let result = curr;
    for (const [k, v] of Object.entries(params)) {
      result = result.replace(new RegExp(`\\{${k}\\}`, 'g'), String(v));
    }
    return result;
  }

  return curr;
}
