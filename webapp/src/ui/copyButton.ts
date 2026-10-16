/**
 * Reusable Code-Block Style Copy Button
 * Provides icon-driven copy action with temporary animated checkmark feedback.
 */

export interface CopyButtonProps {
  getCopyText: () => string;
  title?: string;
  ariaLabel?: string;
  id?: string;
  className?: string;
  disabled?: boolean;
  feedbackText?: string;
  onCopied?: () => void;
}

const COPY_ICON_SVG = `
  <svg class="copy-icon" width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
    <rect x="9" y="9" width="13" height="13" rx="2" ry="2"></rect>
    <path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"></path>
  </svg>
`;

export function renderCopyButton(props: CopyButtonProps): HTMLButtonElement {
  const btn = document.createElement('button');
  btn.type = 'button';
  btn.className = `copy-code-btn ${props.className || ''}`.trim();
  if (props.id) btn.id = props.id;
  if (props.title) btn.title = props.title;
  if (props.ariaLabel) btn.setAttribute('aria-label', props.ariaLabel);
  btn.disabled = !!props.disabled;

  const feedbackText = props.feedbackText || 'Copied';

  btn.innerHTML = `
    ${COPY_ICON_SVG}
    <span class="copy-feedback hidden">✓ ${feedbackText}</span>
  `;

  let resetTimeout: number | null = null;

  btn.addEventListener('click', (e) => {
    e.stopPropagation();
    if (btn.disabled) return;

    const textToCopy = props.getCopyText();
    if (!textToCopy) return;

    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(textToCopy).catch(() => {});
    }

    const icon = btn.querySelector('.copy-icon');
    const feedback = btn.querySelector('.copy-feedback');

    btn.classList.add('copied');
    icon?.classList.add('hidden');
    feedback?.classList.remove('hidden');

    if (resetTimeout !== null) {
      window.clearTimeout(resetTimeout);
    }

    resetTimeout = window.setTimeout(() => {
      btn.classList.remove('copied');
      icon?.classList.remove('hidden');
      feedback?.classList.add('hidden');
      resetTimeout = null;
    }, 1500);

    props.onCopied?.();
  });

  return btn;
}
