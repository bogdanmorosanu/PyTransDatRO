/**
 * Reusable Select Dropdown Control UI Factory
 */

export interface SelectOption<T extends string = string> {
  value: T;
  label: string;
}

export interface SelectControlProps<T extends string = string> {
  options: SelectOption<T>[];
  selectedValue: T;
  onChange: (value: T) => void;
  id?: string;
  ariaLabel?: string;
  className?: string;
  style?: Partial<CSSStyleDeclaration>;
}

export function renderSelect<T extends string = string>(
  props: SelectControlProps<T>
): HTMLSelectElement {
  const select = document.createElement('select');
  select.className = `form-select ${props.className || ''}`.trim();
  if (props.id) select.id = props.id;
  if (props.ariaLabel) select.setAttribute('aria-label', props.ariaLabel);

  if (props.style) {
    Object.assign(select.style, props.style);
  }

  props.options.forEach((opt) => {
    const option = document.createElement('option');
    option.value = opt.value;
    option.textContent = opt.label;
    if (opt.value === props.selectedValue) {
      option.selected = true;
    }
    select.appendChild(option);
  });

  select.addEventListener('change', () => {
    props.onChange(select.value as T);
  });

  return select;
}
