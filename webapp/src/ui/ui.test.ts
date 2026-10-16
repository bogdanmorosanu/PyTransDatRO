import { describe, it, expect, vi, beforeEach } from 'vitest';
import {
  renderDropdownMenu,
  renderSegmentedControl,
  renderSelect,
  renderCopyBadge,
  renderCopyButton,
} from './index';

// Minimal DOM mock environment for Node.js Vitest runner
class MockClassList {
  private classes: Set<string> = new Set();

  constructor(initial = '') {
    initial.split(/\s+/).filter(Boolean).forEach((c) => this.classes.add(c));
  }

  add(...tokens: string[]) {
    tokens.forEach((t) => this.classes.add(t));
  }

  remove(...tokens: string[]) {
    tokens.forEach((t) => this.classes.delete(t));
  }

  contains(token: string) {
    return this.classes.has(token);
  }

  toString() {
    return Array.from(this.classes).join(' ');
  }
}

class MockElement {
  public tagName: string;
  public id: string = '';
  public title: string = '';
  public tabIndex: number = 0;
  public disabled: boolean = false;
  public selected: boolean = false;
  public value: string = '';
  public open: boolean = false;
  public type: string = '';
  public href: string = '';
  public target: string = '';
  public rel: string = '';
  public style: Record<string, any> = {};
  public classList: MockClassList;
  public children: MockElement[] = [];
  public parentNode: MockElement | null = null;
  public listeners: Record<string, ((e?: any) => void)[]> = {};
  private attrs: Record<string, string> = {};
  private _textContent: string = '';
  private _innerHTML: string = '';

  constructor(tagName: string) {
    this.tagName = tagName.toUpperCase();
    this.classList = new MockClassList();
  }

  get className(): string {
    return this.classList.toString();
  }

  set className(val: string) {
    this.classList = new MockClassList(val);
  }

  get textContent(): string {
    return this._textContent;
  }

  set textContent(val: string) {
    this._textContent = val;
  }

  get innerHTML(): string {
    return this._innerHTML;
  }

  set innerHTML(val: string) {
    this._innerHTML = val;
    // Extract any class from simple HTML tags if needed
  }

  setAttribute(name: string, val: string) {
    this.attrs[name] = val;
  }

  getAttribute(name: string): string | null {
    return this.attrs[name] ?? null;
  }

  appendChild<T extends MockElement>(child: T): T {
    child.parentNode = this;
    this.children.push(child);
    return child;
  }

  addEventListener(event: string, callback: (e?: any) => void) {
    if (!this.listeners[event]) this.listeners[event] = [];
    this.listeners[event].push(callback);
  }

  dispatchEvent(event: { type: string; [k: string]: any }) {
    const list = this.listeners[event.type] || [];
    list.forEach((cb) => cb(event));
  }

  click() {
    this.dispatchEvent({ type: 'click', stopPropagation: () => {} });
  }

  focus() {
    (globalThis as any).document.activeElement = this;
    this.dispatchEvent({ type: 'focus' });
  }

  querySelector(selector: string): MockElement | null {
    if (selector.startsWith('.')) {
      const cls = selector.slice(1);
      return this.find((el) => el.classList.contains(cls));
    }
    return null;
  }

  querySelectorAll(selector: string): MockElement[] {
    const results: MockElement[] = [];
    if (selector.startsWith('.')) {
      const cls = selector.slice(1);
      this.findAll((el) => el.classList.contains(cls), results);
    }
    return results;
  }

  private find(predicate: (el: MockElement) => boolean): MockElement | null {
    for (const child of this.children) {
      if (predicate(child)) return child;
      const sub = child.find(predicate);
      if (sub) return sub;
    }
    return null;
  }

  private findAll(predicate: (el: MockElement) => boolean, acc: MockElement[]) {
    for (const child of this.children) {
      if (predicate(child)) acc.push(child);
      child.findAll(predicate, acc);
    }
  }
}

// Setup global mock document before tests
beforeEach(() => {
  const mockDoc = {
    activeElement: null as MockElement | null,
    createElement: (tag: string) => new MockElement(tag),
  };
  (globalThis as any).document = mockDoc;
  try {
    Object.defineProperty(globalThis.navigator, 'clipboard', {
      value: {
        writeText: vi.fn().mockResolvedValue(undefined),
      },
      configurable: true,
      writable: true,
    });
  } catch {
    // fallback if navigator is plain object
  }
  (globalThis as any).window = {
    setTimeout: (fn: Function) => {
      fn();
      return 1;
    },
    clearTimeout: () => {},
  };
});

describe('UI Factory: renderDropdownMenu', () => {
  it('renders a details dropdown with nav-menu-trigger by default', () => {
    const onClickMock = vi.fn();
    const menu = renderDropdownMenu({
      id: 'test-menu',
      triggerLabel: 'Help & Resources',
      items: [
        { type: 'group', title: 'Category' },
        { type: 'action', label: 'Item 1', onClick: onClickMock, id: 'action-1' },
        { type: 'divider' },
        { type: 'link', label: 'Docs', href: '/docs', isExternal: true },
      ],
    }) as any;

    expect(menu.id).toBe('test-menu');
    expect(menu.classList.contains('header-dropdown-wrap')).toBe(true);

    const summary = menu.children[0];
    expect(summary.classList.contains('nav-menu-trigger')).toBe(true);
    expect(summary.getAttribute('role')).toBe('button');
    expect(summary.getAttribute('aria-expanded')).toBe('false');

    const menuPanel = menu.children[1];
    expect(menuPanel.classList.contains('header-dropdown-menu')).toBe(true);
    expect(menuPanel.getAttribute('role')).toBe('menu');

    // Trigger action item click
    const actionBtn = menuPanel.children.find((c: any) => c.id === 'action-1');
    expect(actionBtn).toBeDefined();
    actionBtn.click();
    expect(onClickMock).toHaveBeenCalledTimes(1);
    expect(menu.open).toBe(false);
  });

  it('supports btn-outline triggerVariant when specified', () => {
    const menu = renderDropdownMenu({
      triggerLabel: 'Options',
      triggerVariant: 'btn-outline',
      items: [],
    }) as any;

    const summary = menu.children[0];
    expect(summary.classList.contains('btn-outline')).toBe(true);
  });
});

describe('UI Factory: renderSegmentedControl', () => {
  it('renders radio group buttons and dispatches change events', () => {
    const onChange = vi.fn();
    const control = renderSegmentedControl({
      options: [
        { value: 'NE', label: 'North-East' },
        { value: 'EN', label: 'East-North' },
      ],
      activeValue: 'NE',
      onChange,
      ariaLabel: 'Coordinate Order',
    }) as any;

    expect(control.getAttribute('role')).toBe('radiogroup');
    expect(control.getAttribute('aria-label')).toBe('Coordinate Order');
    expect(control.children.length).toBe(2);

    const [btn1, btn2] = control.children;
    expect(btn1.classList.contains('active')).toBe(true);
    expect(btn1.getAttribute('aria-checked')).toBe('true');
    expect(btn2.classList.contains('active')).toBe(false);
    expect(btn2.getAttribute('aria-checked')).toBe('false');

    // Click button 2
    btn2.click();
    expect(onChange).toHaveBeenCalledWith('EN');
    expect(btn2.classList.contains('active')).toBe(true);
    expect(btn1.classList.contains('active')).toBe(false);
  });
});

describe('UI Factory: renderSelect', () => {
  it('renders select dropdown and sets initial selected value', () => {
    const onChange = vi.fn();
    const select = renderSelect({
      options: [
        { value: 'degrees', label: 'Degrees' },
        { value: 'dms', label: 'DMS' },
        { value: 'radians', label: 'Radians' },
      ],
      selectedValue: 'dms',
      onChange,
      ariaLabel: 'Angular Unit',
    }) as any;

    expect(select.classList.contains('form-select')).toBe(true);
    expect(select.getAttribute('aria-label')).toBe('Angular Unit');
    expect(select.children.length).toBe(3);

    const dmsOption = select.children[1];
    expect(dmsOption.selected).toBe(true);

    select.value = 'radians';
    select.dispatchEvent({ type: 'change' });
    expect(onChange).toHaveBeenCalledWith('radians');
  });
});

describe('UI Factory: renderCopyBadge', () => {
  it('renders copyable badge and triggers clipboard write on click and keydown', () => {
    const onCopied = vi.fn();
    const badge = renderCopyBadge({
      text: '45.84123450',
      tooltip: 'Click to copy',
      onCopied,
    }) as any;

    expect(badge.classList.contains('format-example-badge')).toBe(true);
    expect(badge.getAttribute('role')).toBe('button');
    expect(badge.tabIndex).toBe(0);
    expect(badge.title).toBe('Click to copy');

    // Click
    badge.click();
    expect(navigator.clipboard.writeText).toHaveBeenCalledWith('45.84123450');
    expect(onCopied).toHaveBeenCalledTimes(1);

    // Keydown Enter
    badge.dispatchEvent({ type: 'keydown', key: 'Enter', preventDefault: () => {} });
    expect(onCopied).toHaveBeenCalledTimes(2);
  });
});

describe('UI Factory: renderCopyButton', () => {
  it('renders copy button and handles text copy and temporary state toggle', () => {
    const onCopied = vi.fn();
    const btn = renderCopyButton({
      getCopyText: () => '100.000, 200.000',
      title: 'Copy coordinates',
      onCopied,
    }) as any;

    expect(btn.classList.contains('copy-code-btn')).toBe(true);
    expect(btn.title).toBe('Copy coordinates');

    btn.click();
    expect(navigator.clipboard.writeText).toHaveBeenCalledWith('100.000, 200.000');
    expect(onCopied).toHaveBeenCalledTimes(1);
  });
});
