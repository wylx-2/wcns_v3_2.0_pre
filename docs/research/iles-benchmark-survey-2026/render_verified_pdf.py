"""Use the packaged renderer rasterization after a native Word PDF export.

The regular render_docx.py failed because LibreOffice is not installed.
Word's PDF is the verified conversion backend for this Windows host; the
packaged renderer still performs page rasterization through bundled Poppler.
"""
from pathlib import Path
import importlib.util, os, json
from pypdf import PdfReader
R=Path(__file__).resolve().parent
STEM='ILES湍流验证算例调研报告'
SCRIPT=Path(r'C:\Users\15018\.codex\plugins\cache\openai-primary-runtime\documents\26.909.12148\skills\documents\render_docx.py')
os.environ['PATH']=r'C:\Users\15018\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin'+os.pathsep+os.environ['PATH']
spec=importlib.util.spec_from_file_location('packaged_renderer',SCRIPT)
renderer=importlib.util.module_from_spec(spec);spec.loader.exec_module(renderer)
pdf=R/'qa'/(STEM+'.pdf')
assert pdf.is_file() and pdf.stat().st_size>0
renderer.convert_to_pdf=lambda *args,**kwargs:(str(pdf),'Native Word ExportAsFixedFormat backend; LibreOffice unavailable')
renderer.rasterize(str(R/(STEM+'.docx')),str(R/'qa'/'final-pages'),144,False,False)
reader=PdfReader(pdf)
pages=[]
for i,page in enumerate(reader.pages,1):
 text=page.extract_text() or ''
 pages.append({'page':i,'chars':len(text),'beginning':text[:100],'ending':text[-100:]})
(R/'qa'/'page_text_audit.json').write_text(json.dumps(pages,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'pages':len(pages),'min_chars':min(p['chars'] for p in pages),'backend':'Word PDF + packaged render_docx.py rasterize'},ensure_ascii=True))
