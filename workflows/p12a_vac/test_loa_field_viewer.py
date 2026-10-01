"""Exercise full-atlas rendering on real smoke tiles in a disposable fixture.

The fixture is explicitly incomplete, never a science atlas or delivered HTML.
"""
import argparse
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
from workflows.p12a_vac import loa_field_viewer as viewer
from workflows.p12a_vac import check_loa_viewer as qa


def main(a):
    plan=json.loads((a.root/'PLAN.json').read_text())
    smoke=json.loads((a.root/'SMOKE_COMPLETE.json').read_text())
    names={r['case']['id'] for r in smoke['results']}
    plan['tiles']=[c for c in plan['tiles'] if c['id'] in names]
    assert len(plan['tiles'])==3 and {c['cap'] for c in plan['tiles']}=={'NGC','SGC'}
    with tempfile.TemporaryDirectory(prefix='loa-viewer-INCOMPLETE-TEST-') as temp:
        root=Path(temp)
        for name in ['vac_properties.npz','details','tiles']:(root/name).symlink_to(a.root/name)
        (root/'PLAN.json').write_text(json.dumps(plan))
        (root/'ATLAS_COMPLETE.json').write_text(json.dumps({'passed':True,'tiles':3,'TEST_FIXTURE_ONLY':True}))
        viewer.build(SimpleNamespace(root=root,out=root/'viewer',details_only=False,surface_library=a.root/'viewer_dependencies'))
        qa.main(SimpleNamespace(html=root/'viewer'/'Loa_3D.html',out=a.out,playwright_path=a.root/'playwright_runtime'))
    (a.out/'SMOKE_SCOPE.json').write_text(json.dumps({'test_only':True,'real_smoke_tiles':sorted(names),
        'qualification':'Disposable three-tile test of atlas viewer code; no complete atlas claimed. Test HTML removed.'},indent=2)+'\n')
    print('Both-cap atlas rendering smoke passed; disposable HTML removed.',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    main(p.parse_args())
