import { describe, it, expect, beforeEach } from 'vitest';
import { en } from './en';
import { ro } from './ro';
import { t, setLanguage, getLanguage, toggleLanguage, onLanguageChange } from './index';

function extractLeaves(obj: Record<string, any>, prefix = ''): Map<string, string> {
  const result = new Map<string, string>();
  for (const [key, value] of Object.entries(obj)) {
    const fullKey = prefix ? `${prefix}.${key}` : key;
    if (value !== null && typeof value === 'object' && !Array.isArray(value)) {
      const nested = extractLeaves(value, fullKey);
      for (const [nKey, nVal] of nested.entries()) {
        result.set(nKey, nVal);
      }
    } else if (typeof value === 'string') {
      result.set(fullKey, value);
    }
  }
  return result;
}

function extractPlaceholders(text: string): string[] {
  const matches = text.match(/\{([a-zA-Z0-9_]+)\}/g);
  if (!matches) return [];
  return matches.map((m) => m.slice(1, -1)).sort();
}

describe('i18n Translation Dictionaries Parity', () => {
  const enLeaves = extractLeaves(en);
  const roLeaves = extractLeaves(ro);

  it('contains identical leaf keys across en and ro', () => {
    const enKeys = Array.from(enLeaves.keys()).sort();
    const roKeys = Array.from(roLeaves.keys()).sort();

    const missingInRo = enKeys.filter((k) => !roLeaves.has(k));
    const missingInEn = roKeys.filter((k) => !enLeaves.has(k));

    expect(missingInRo, `Keys missing in Romanian: ${missingInRo.join(', ')}`).toEqual([]);
    expect(missingInEn, `Keys missing in English: ${missingInEn.join(', ')}`).toEqual([]);
  });

  it('maintains identical variable placeholders in matching strings', () => {
    const mismatchedPlaceholders: string[] = [];

    for (const [key, enText] of enLeaves.entries()) {
      const roText = roLeaves.get(key);
      if (!roText) continue;

      const enHolders = extractPlaceholders(enText);
      const roHolders = extractPlaceholders(roText);

      if (enHolders.join(',') !== roHolders.join(',')) {
        mismatchedPlaceholders.push(
          `${key}: EN [${enHolders.join(', ')}] vs RO [${roHolders.join(', ')}]`
        );
      }
    }

    expect(
      mismatchedPlaceholders,
      `Mismatched placeholders: \n${mismatchedPlaceholders.join('\n')}`
    ).toEqual([]);
  });
});

describe('i18n Runtime & Fallback Behavior', () => {
  beforeEach(() => {
    setLanguage('ro');
  });

  it('translates nested keys in both languages', () => {
    setLanguage('ro');
    expect(t('app.title')).toBe('PyTransDatRO');
    expect(t('common.close')).toBe('Închide');

    setLanguage('en');
    expect(t('common.close')).toBe('Close');
  });

  it('interpolates parameters correctly into placeholders', () => {
    setLanguage('en');
    const msg = t('volumeWarning.message', { count: '100' });
    expect(msg).toContain('100');

    setLanguage('ro');
    const roMsg = t('volumeWarning.message', { count: '100' });
    expect(roMsg).toContain('100');
  });

  it('returns key path if key is not found', () => {
    expect(t('nonexistent.deep.key')).toBe('nonexistent.deep.key');
  });

  it('switches languages via toggleLanguage and notifies listeners', () => {
    setLanguage('ro');
    let notifiedLang = '';
    const unsubscribe = onLanguageChange((lang) => {
      notifiedLang = lang;
    });

    const next = toggleLanguage();
    expect(next).toBe('en');
    expect(getLanguage()).toBe('en');
    expect(notifiedLang).toBe('en');

    unsubscribe();
    toggleLanguage();
    expect(notifiedLang).toBe('en'); // Listener was detached
  });
});
