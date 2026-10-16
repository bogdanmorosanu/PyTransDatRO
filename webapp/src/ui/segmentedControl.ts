/**
 * Reusable Segmented Control UI Factory
 * Provides uniform radio-switch style segmented buttons with keyboard navigation & ARIA.
 */

export interface SegmentedOption<T extends string = string> {
  value: T;
  label: string;
  title?: string;
  id?: string;
}

export interface SegmentedControlProps<T extends string = string> {
  options: SegmentedOption<T>[];
  activeValue: T;
  onChange: (value: T) => void;
  id?: string;
  ariaLabel?: string;
  className?: string;
  style?: Partial<CSSStyleDeclaration>;
}

export function renderSegmentedControl<T extends string = string>(
  props: SegmentedControlProps<T>
): HTMLElement {
  const container = document.createElement('div');
  container.className = `segmented-control ${props.className || ''}`.trim();
  container.setAttribute('role', 'radiogroup');
  if (props.id) container.id = props.id;
  if (props.ariaLabel) container.setAttribute('aria-label', props.ariaLabel);

  if (props.style) {
    Object.assign(container.style, props.style);
  }

  const buttons: HTMLButtonElement[] = [];

  const updateActiveStates = (selectedValue: T) => {
    buttons.forEach((btn, index) => {
      const opt = props.options[index];
      const isSelected = opt.value === selectedValue;
      if (isSelected) {
        btn.classList.add('active');
        btn.setAttribute('aria-checked', 'true');
        btn.tabIndex = 0;
      } else {
        btn.classList.remove('active');
        btn.setAttribute('aria-checked', 'false');
        btn.tabIndex = -1;
      }
    });
  };

  props.options.forEach((opt) => {
    const btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'segmented-btn';
    btn.setAttribute('role', 'radio');
    btn.textContent = opt.label;
    if (opt.id) btn.id = opt.id;
    if (opt.title) btn.title = opt.title;

    btn.addEventListener('click', () => {
      if (btn.classList.contains('active')) return;
      updateActiveStates(opt.value);
      props.onChange(opt.value);
    });

    buttons.push(btn);
    container.appendChild(btn);
  });

  updateActiveStates(props.activeValue);

  // Keyboard navigation for radio groups (Left/Right, Up/Down)
  container.addEventListener('keydown', (e) => {
    const currentIndex = buttons.indexOf(document.activeElement as HTMLButtonElement);
    if (currentIndex === -1) return;

    let nextIndex = -1;
    if (e.key === 'ArrowRight' || e.key === 'ArrowDown') {
      e.preventDefault();
      nextIndex = currentIndex < buttons.length - 1 ? currentIndex + 1 : 0;
    } else if (e.key === 'ArrowLeft' || e.key === 'ArrowUp') {
      e.preventDefault();
      nextIndex = currentIndex > 0 ? currentIndex - 1 : buttons.length - 1;
    }

    if (nextIndex !== -1) {
      const targetOpt = props.options[nextIndex];
      updateActiveStates(targetOpt.value);
      buttons[nextIndex].focus();
      props.onChange(targetOpt.value);
    }
  });

  return container;
}
