/**
 * Reusable Dropdown Menu & Menu Item Factory
 * Supports both sleek navbar item triggers (.nav-menu-trigger) and classic outlined triggers.
 */

export type MenuItem =
  | {
      type: 'action';
      label: string;
      iconSvg?: string;
      onClick: () => void;
      id?: string;
      danger?: boolean;
    }
  | {
      type: 'link';
      label: string;
      href: string;
      isExternal?: boolean;
      iconSvg?: string;
      id?: string;
    }
  | {
      type: 'group';
      title: string;
    }
  | {
      type: 'divider';
      dashed?: boolean;
    };

export interface DropdownMenuProps {
  id?: string;
  triggerLabel: string;
  triggerIconSvg?: string;
  triggerVariant?: 'menu-item' | 'btn-outline';
  items: MenuItem[];
  align?: 'left' | 'right';
  className?: string;
  ariaLabel?: string;
  onOpenChange?: (isOpen: boolean) => void;
}

const EXTERNAL_LINK_ICON = `
  <svg class="ext-icon" width="11" height="11" viewBox="0 0 12 12" fill="none">
    <path d="M3.5 2H10V8.5M10 2L2 10" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>
  </svg>
`;

const CHEVRON_ICON = `
  <svg class="chevron-icon" width="12" height="12" viewBox="0 0 12 12" fill="none">
    <path d="M2.5 4.5L6 8L9.5 4.5" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round"/>
  </svg>
`;

export function renderDropdownMenu(props: DropdownMenuProps): HTMLDetailsElement {
  const details = document.createElement('details');
  details.className = `header-dropdown-wrap ${props.className || ''}`.trim();
  if (props.id) details.id = props.id;

  // Trigger Summary
  const summary = document.createElement('summary');
  const isOutline = props.triggerVariant === 'btn-outline';
  summary.className = isOutline ? 'btn btn-outline header-dropdown-trigger' : 'nav-menu-trigger';
  summary.setAttribute('role', 'button');
  summary.setAttribute('aria-haspopup', 'true');
  summary.setAttribute('aria-expanded', 'false');
  if (props.ariaLabel) summary.setAttribute('aria-label', props.ariaLabel);

  let triggerInner = '';
  if (props.triggerIconSvg) {
    triggerInner += props.triggerIconSvg;
  }
  triggerInner += `<span>${props.triggerLabel}</span>`;
  triggerInner += CHEVRON_ICON;
  summary.innerHTML = triggerInner;
  details.appendChild(summary);

  // Dropdown Menu Container
  const menuContainer = document.createElement('div');
  menuContainer.className = `header-dropdown-menu ${props.align === 'left' ? 'align-left' : ''}`.trim();
  menuContainer.setAttribute('role', 'menu');

  props.items.forEach((item) => {
    if (item.type === 'action') {
      const btn = document.createElement('button');
      btn.type = 'button';
      btn.className = `dropdown-item ${item.danger ? 'danger' : ''}`.trim();
      btn.setAttribute('role', 'menuitem');
      if (item.id) btn.id = item.id;

      let btnHtml = '';
      if (item.iconSvg) {
        btnHtml += `<span class="menu-item-left">${item.iconSvg}<span>${item.label}</span></span>`;
      } else {
        btnHtml += `<span>${item.label}</span>`;
      }
      btn.innerHTML = btnHtml;

      btn.addEventListener('click', (e) => {
        e.stopPropagation();
        details.open = false;
        summary.setAttribute('aria-expanded', 'false');
        props.onOpenChange?.(false);
        item.onClick();
      });

      menuContainer.appendChild(btn);
    } else if (item.type === 'link') {
      const link = document.createElement('a');
      link.href = item.href;
      link.className = 'dropdown-item';
      link.setAttribute('role', 'menuitem');
      if (item.id) link.id = item.id;

      if (item.isExternal) {
        link.target = '_blank';
        link.rel = 'noopener noreferrer';
      }

      let linkHtml = '';
      if (item.iconSvg) {
        linkHtml += `<span class="menu-item-left">${item.iconSvg}<span>${item.label}</span></span>`;
      } else {
        linkHtml += `<span>${item.label}</span>`;
      }

      if (item.isExternal) {
        linkHtml += EXTERNAL_LINK_ICON;
      }
      link.innerHTML = linkHtml;

      link.addEventListener('click', () => {
        details.open = false;
        summary.setAttribute('aria-expanded', 'false');
        props.onOpenChange?.(false);
      });

      menuContainer.appendChild(link);
    } else if (item.type === 'group') {
      const title = document.createElement('div');
      title.className = 'dropdown-group-title';
      title.setAttribute('role', 'presentation');
      title.textContent = item.title;
      menuContainer.appendChild(title);
    } else if (item.type === 'divider') {
      const divider = document.createElement('div');
      divider.className = item.dashed ? 'dropdown-divider-dashed' : 'dropdown-divider';
      divider.setAttribute('role', 'separator');
      menuContainer.appendChild(divider);
    }
  });

  details.appendChild(menuContainer);

  // Toggle & State Tracking
  details.addEventListener('toggle', () => {
    summary.setAttribute('aria-expanded', details.open ? 'true' : 'false');
    props.onOpenChange?.(details.open);
  });

  // Keyboard navigation inside dropdown items
  details.addEventListener('keydown', (e) => {
    if (e.key === 'Escape' && details.open) {
      e.preventDefault();
      details.open = false;
      summary.setAttribute('aria-expanded', 'false');
      summary.focus();
      props.onOpenChange?.(false);
      return;
    }

    if (!details.open) return;

    const focusableItems = Array.from(
      menuContainer.querySelectorAll<HTMLElement>('.dropdown-item')
    );
    if (focusableItems.length === 0) return;

    const currentIndex = focusableItems.indexOf(document.activeElement as HTMLElement);

    if (e.key === 'ArrowDown') {
      e.preventDefault();
      const nextIndex = currentIndex < focusableItems.length - 1 ? currentIndex + 1 : 0;
      focusableItems[nextIndex].focus();
    } else if (e.key === 'ArrowUp') {
      e.preventDefault();
      const prevIndex = currentIndex > 0 ? currentIndex - 1 : focusableItems.length - 1;
      focusableItems[prevIndex].focus();
    }
  });

  return details;
}
