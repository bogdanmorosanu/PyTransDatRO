import { ro } from './ro';
import { en } from './en';

export type Language = 'ro' | 'en';

const dictionaries = { ro, en };

const getInitialLanguage = (): Language => {
  if (typeof localStorage !== 'undefined') {
    try {
      return (localStorage.getItem('pytransdat_lang') as Language) || 'ro';
    } catch {
      return 'ro';
    }
  }
  return 'ro';
};

let currentLang: Language = getInitialLanguage();
const listeners: Array<(lang: Language) => void> = [];

export function getLanguage(): Language {
  return currentLang;
}

export function setLanguage(lang: Language): void {
  if (lang !== currentLang) {
    currentLang = lang;
    if (typeof localStorage !== 'undefined') {
      try {
        localStorage.setItem('pytransdat_lang', lang);
      } catch {
        // ignore storage errors
      }
    }
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
 * Translates a key path using the active language dictionary.
 * Supports dot notation: t('app.title'), t('point.inputCoordsHint')
 * Supports optional parameter interpolation: t('volumeWarning.message', { count: 5000 })
 */
export function t(path: string, params?: Record<string, string | number>): string {
  const dict = dictionaries[currentLang] || dictionaries.ro;
  const parts = path.split('.');
  let current: any = dict;

  for (const part of parts) {
    if (current && typeof current === 'object' && part in current) {
      current = current[part];
    } else {
      // Fallback to Romanian dictionary if key is missing in active language
      let fallback: any = dictionaries.ro;
      for (const fbPart of parts) {
        if (fallback && typeof fallback === 'object' && fbPart in fallback) {
          fallback = fallback[fbPart];
        } else {
          return path;
        }
      }
      current = fallback;
      break;
    }
  }

  if (typeof current !== 'string') {
    return path;
  }

  if (params) {
    let result = current;
    for (const [key, value] of Object.entries(params)) {
      result = result.replace(new RegExp(`\\{${key}\\}`, 'g'), String(value));
    }
    return result;
  }

  return current;
}
