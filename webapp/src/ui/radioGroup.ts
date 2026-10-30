/**
 * Reusable Radio Card / Selection Group UI Factory (Archetype 2.5)
 * Provides accessible, keyboard-navigable single-choice radio option cards.
 */

export interface RadioOption<T extends string = string> {
  value: T;
  label: string;
  description?: string;
  badge?: string;
  id?: string;
}

export interface RadioGroupProps<T extends string = string> {
  name: string;
  options: RadioOption<T>[];
  selectedValue: T;
  onChange: (value: T) => void;
  id?: string;
  ariaLabel?: string;
  className?: string;
}

export function renderRadioGroup<T extends string = string>(
  props: RadioGroupProps<T>
): HTMLElement {
  const container = document.createElement('div');
  container.className = `radio-group ${props.className || ''}`.trim();
  container.setAttribute('role', 'radiogroup');
  if (props.id) container.id = props.id;
  if (props.ariaLabel) container.setAttribute('aria-label', props.ariaLabel);

  const cards: HTMLLabelElement[] = [];
  const inputs: HTMLInputElement[] = [];

  const updateSelectedState = (val: T) => {
    inputs.forEach((input, idx) => {
      const isSelected = input.value === val;
      input.checked = isSelected;
      input.setAttribute('aria-checked', String(isSelected));
      input.tabIndex = isSelected ? 0 : -1;
      cards[idx].classList.toggle('active', isSelected);
    });
  };

  props.options.forEach((opt) => {
    const card = document.createElement('label');
    card.className = 'radio-card';
    if (opt.id) card.setAttribute('for', opt.id);

    const input = document.createElement('input');
    input.type = 'radio';
    input.name = props.name;
    input.value = opt.value;
    input.className = 'radio-input';
    if (opt.id) input.id = opt.id;
    input.setAttribute('role', 'radio');

    const circle = document.createElement('span');
    circle.className = 'radio-circle';
    circle.setAttribute('aria-hidden', 'true');

    const content = document.createElement('div');
    content.className = 'radio-content';

    const titleRow = document.createElement('div');
    titleRow.className = 'radio-title-row';

    const labelSpan = document.createElement('span');
    labelSpan.className = 'radio-label';
    labelSpan.textContent = opt.label;
    titleRow.appendChild(labelSpan);

    if (opt.badge) {
      const badgeSpan = document.createElement('span');
      badgeSpan.className = 'radio-badge';
      badgeSpan.textContent = opt.badge;
      titleRow.appendChild(badgeSpan);
    }
    content.appendChild(titleRow);

    if (opt.description) {
      const descSpan = document.createElement('span');
      descSpan.className = 'radio-desc';
      descSpan.textContent = opt.description;
      content.appendChild(descSpan);
    }

    card.appendChild(input);
    card.appendChild(circle);
    card.appendChild(content);

    input.addEventListener('change', () => {
      if (input.checked) {
        updateSelectedState(opt.value);
        props.onChange(opt.value);
      }
    });

    cards.push(card);
    inputs.push(input);
    container.appendChild(card);
  });

  // Set initial selected state
  updateSelectedState(props.selectedValue);

  // Keyboard navigation for radio group (ArrowUp/ArrowDown, ArrowLeft/ArrowRight)
  container.addEventListener('keydown', (e) => {
    const activeIndex = inputs.findIndex((input) => input === document.activeElement);
    if (activeIndex === -1) return;

    let targetIndex = -1;
    if (e.key === 'ArrowDown' || e.key === 'ArrowRight') {
      e.preventDefault();
      targetIndex = activeIndex < inputs.length - 1 ? activeIndex + 1 : 0;
    } else if (e.key === 'ArrowUp' || e.key === 'ArrowLeft') {
      e.preventDefault();
      targetIndex = activeIndex > 0 ? activeIndex - 1 : inputs.length - 1;
    }

    if (targetIndex !== -1) {
      const targetOpt = props.options[targetIndex];
      updateSelectedState(targetOpt.value);
      inputs[targetIndex].focus();
      props.onChange(targetOpt.value);
    }
  });

  return container;
}
