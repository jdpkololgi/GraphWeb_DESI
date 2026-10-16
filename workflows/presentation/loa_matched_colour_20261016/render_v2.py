import os
import socket
from run import OUT,render,sha,write
from pathlib import Path
assert os.environ.get('SLURM_JOB_ID') and socket.gethostname().startswith('nid')
render('figures_v2','histograms_v2.csv')
assert sha(OUT/'histograms.csv')==sha(OUT/'histograms_v2.csv')
write(OUT/'PRESENTATION_V2.json',dict(code_sha256=sha(Path(__file__).with_name('run.py')),
    original_manifest_sha256=sha(OUT/'MANIFEST.json'),
    change='Shared density limits cover all panels; annotation headroom; symmetric residual limits',
    files={str(p.relative_to(OUT)):sha(p) for p in sorted((OUT/'figures_v2').iterdir())}))
