import { t, onLanguageChange } from '../i18n';
import { apiClient, type TransformationOp } from '../utils/apiClient';
import { showToast } from '../utils/toast';

export class FileUpload {
  private container: HTMLElement;
  private op: TransformationOp = 'Stereo70ToETRS89';
  private selectedFile: File | null = null;
  private previewLines: string[] = [];
  private isUploading: boolean = false;

  constructor() {
    this.container = document.createElement('div');
    this.container.className = 'tab-content';
    this.render();

    onLanguageChange(() => this.render());
  }

  public getElement(): HTMLElement {
    return this.container;
  }

  private handleFileSelected(file: File): void {
    if (file.size > 20 * 1024 * 1024) {
      showToast('File size exceeds 20 MB limit', 'danger');
      return;
    }

    this.selectedFile = file;

    // Read first 5 lines for preview
    const reader = new FileReader();
    reader.onload = (e) => {
      const text = (e.target?.result as string) || '';
      const lines = text.split(/\r\n|\r|\n/).filter((l) => l.trim().length > 0);
      this.previewLines = lines.slice(0, 5);
      this.render();
    };
    reader.readAsText(file.slice(0, 4096)); // Read first 4 KB
  }

  private async startStreamingUpload(): Promise<void> {
    if (!this.selectedFile || this.isUploading) return;
    this.isUploading = true;
    this.render();

    try {
      showToast(t('file.streamingProgress'), 'info');
      const response = await apiClient.transformFile(this.selectedFile, this.op, 'degrees', 'auto');

      // Trigger streaming download directly from response blob
      const blob = await response.blob();
      const downloadUrl = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = downloadUrl;
      a.download = `transformed_${this.selectedFile.name}`;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(downloadUrl);

      showToast('File transformed and downloaded successfully!', 'success');
    } catch (err: any) {
      showToast(err.message || 'Error during streaming upload', 'danger');
    } finally {
      this.isUploading = false;
      this.render();
    }
  }

  private render(): void {
    const isS70 = this.op === 'Stereo70ToETRS89';

    this.container.innerHTML = `
      <div class="studio-card">
        <div class="card-title">
          <span>${t('file.title')}</span>
        </div>

        <!-- Direction Switch -->
        <div class="form-group" style="margin-bottom: 14px;">
          <label class="form-label">${t('point.opLabel')}</label>
          <div class="segmented-control">
            <button class="segmented-btn ${isS70 ? 'active' : ''}" id="file-op-s70">${t('point.opS70ToETRS')}</button>
            <button class="segmented-btn ${!isS70 ? 'active' : ''}" id="file-op-etrs">${t('point.opETRSToS70')}</button>
          </div>
        </div>

        <!-- Drag and Drop Dropzone -->
        <div 
          id="file-dropzone" 
          style="
            border: 2px dashed var(--border-medium); 
            border-radius: var(--radius-md); 
            padding: 32px 20px; 
            text-align: center; 
            cursor: pointer;
            background-color: var(--bg-surface-input);
            transition: all var(--transition-fast);
          "
        >
          <div style="font-size: 32px; margin-bottom: 8px;">📂</div>
          <div style="font-weight: 600; font-size: 14px; margin-bottom: 4px; color: var(--text-main);">
            ${this.selectedFile ? this.selectedFile.name : t('file.dropzoneTitle')}
          </div>
          <div style="color: var(--text-muted); font-size: 12px;">
            ${this.selectedFile ? `${(this.selectedFile.size / 1024).toFixed(1)} KB` : t('file.dropzoneHint')}
          </div>
          <input type="file" id="file-input-hidden" accept=".csv,.txt,.xyz" style="display: none;" />
        </div>

        <!-- First 5 lines preview -->
        ${
          this.previewLines.length > 0
            ? `
          <div style="margin-top: 16px;">
            <div class="form-label">${t('file.previewTitle')}</div>
            <div class="form-input mono" style="font-size: 11px; white-space: pre; overflow-x: auto; max-height: 120px;">
${this.previewLines.join('\n')}
            </div>
          </div>
        `
            : ''
        }

        <!-- Upload Action Button -->
        <button 
          class="btn btn-primary" 
          id="file-start-btn" 
          style="width: 100%; margin-top: 16px;" 
          ${!this.selectedFile ? 'disabled' : ''}
        >
          <span>${this.isUploading ? '⏳ ' + t('file.streamingProgress') : '🚀 ' + t('file.uploadAndTransform')}</span>
        </button>
      </div>
    `;

    this.attachEventListeners();
  }

  private attachEventListeners(): void {
    // Op switches
    this.container.querySelector('#file-op-s70')?.addEventListener('click', () => {
      this.op = 'Stereo70ToETRS89';
      this.render();
    });
    this.container.querySelector('#file-op-etrs')?.addEventListener('click', () => {
      this.op = 'ETRS89ToStereo70';
      this.render();
    });

    // Dropzone
    const dropzone = this.container.querySelector('#file-dropzone');
    const hiddenInput = this.container.querySelector('#file-input-hidden') as HTMLInputElement;

    if (dropzone && hiddenInput) {
      dropzone.addEventListener('click', () => hiddenInput.click());

      dropzone.addEventListener('dragover', (e: any) => {
        e.preventDefault();
        dropzone.setAttribute('style', `${dropzone.getAttribute('style')}; border-color: var(--color-primary); background-color: var(--color-primary-glow);`);
      });

      dropzone.addEventListener('dragleave', (e: any) => {
        e.preventDefault();
        this.render();
      });

      dropzone.addEventListener('drop', (e: any) => {
        e.preventDefault();
        if (e.dataTransfer && e.dataTransfer.files.length > 0) {
          this.handleFileSelected(e.dataTransfer.files[0]);
        }
      });

      hiddenInput.addEventListener('change', () => {
        if (hiddenInput.files && hiddenInput.files.length > 0) {
          this.handleFileSelected(hiddenInput.files[0]);
        }
      });
    }

    this.container.querySelector('#file-start-btn')?.addEventListener('click', () => {
      this.startStreamingUpload();
    });
  }
}
