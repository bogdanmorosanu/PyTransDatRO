/**
 * Reusable Copyable Badge Component
 * Interactive chip that selects and copies its contents on click or keyboard activation.
 */

export interface CopyBadgeProps {
  text: string;
  copyValue?: string;
  tooltip?: string;
  className?: string;
  id?: string;
  onCopied?: () => void;
}

export function renderCopyBadge(props: CopyBadgeProps): HTMLElement {
  const badge = document.createElement('span');
  badge.className = `format-example-badge ${props.className || ''}`.trim();
  badge.textContent = props.text;
  badge.tabIndex = 0;
  badge.setAttribute('role', 'button');
  if (props.id) badge.id = props.id;
  if (props.tooltip) badge.title = props.tooltip;

  const performCopy = () => {
    const valToCopy = props.copyValue !== undefined ? props.copyValue : props.text;

    // Visual selection range
    try {
      const range = document.createRange();
      range.selectNodeContents(badge);
      const sel = window.getSelection();
      sel?.removeAllRanges();
      sel?.addRange(range);
    } catch {
      // Ignore if DOM selection fails
    }

    // Clipboard write
    if (typeof navigator !== 'undefined' && navigator.clipboard?.writeText) {
      navigator.clipboard.writeText(valToCopy).catch(() => {});
    }

    props.onCopied?.();
  };

  badge.addEventListener('click', performCopy);

  badge.addEventListener('keydown', (e) => {
    if (e.key === 'Enter' || e.key === ' ') {
      e.preventDefault();
      performCopy();
    }
  });

  return badge;
}
