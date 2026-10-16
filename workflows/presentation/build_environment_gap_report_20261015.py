"""Package reviewed local evidence into one offline text-and-images document."""
from pathlib import Path
import base64,hashlib,json,html,re
import markdown
from PIL import Image
H=Path('/global/u2/d/dkololgi'); R=H/'fmsu_tokyo_20261010'; G=H/'GraphWeb_DESI'; I=H/'TNG/Illustris'
D=G/'docs/figures/environment_gap_three_tests_20261015'; O=R/'environment_gap_three_tests_20261015';O.mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def data(p,mime):return 'data:'+mime+';base64,'+base64.b64encode(p.read_bytes()).decode()
def md(s):return markdown.markdown(s,extensions=['tables','fenced_code','toc'])
assets=[]
def fig(p,caption=None):
    with Image.open(p) as im:im.verify()
    assets.append(dict(path=str(p),sha256=sha(p),kind='figure'))
    return '<figure><img loading="lazy" src="'+data(p,'image/png')+'" alt="'+html.escape(caption or p.stem)+'"><figcaption>'+html.escape(caption or p.stem.replace('_',' '))+'</figcaption></figure>'
def pdf(p):
    b=p.read_bytes();assert b.startswith(b'%PDF-') and b.rstrip().endswith(b'%%EOF')
    assets.append(dict(path=str(p),sha256=sha(p),kind='pdf'))
    return '<p><a download="'+p.name+'" href="'+data(p,'application/pdf')+'">Download embedded PDF: '+html.escape(p.name)+'</a> (complete original figure supplement)</p>'
a=json.loads((D/'test1.json').read_text());b=json.loads((D/'test3.json').read_text())
s=(G/'docs/environment_gap_three_tests_20261015.md').read_text()
rows=['| Population | N | R2 contrast [95% interval] (pp) | R7 contrast [95% interval] (pp) |','|---|---:|---:|---:|']
for high in [False,True]:
    rr=[next(r for r in a['arms'] if r['radius']==rad and r['high_mass']==high and r['estimator']=='mass_controlled_trace_quartiles') for rad in [2,7]]
    rows.append('| '+('High mass: 10.5≤log M<11.5' if high else 'Full: 9≤log M<11.5')+f" | {rr[0]['n']:,} | "+' | '.join(f"{100*r['contrast'][0]:.2f} [{100*r['ci95'][0][0]:.2f}, {100*r['ci95'][0][1]:.2f}]" for r in rr)+' |')
s=s.replace('<!-- TEST1_TABLE -->','\n'.join(rows))
rows=['| z range | Common parent N | All eligible (pp) | Low-width half (pp) |','|---|---:|---:|---:|']
for r in b['windows']:
    if r['boundary_over_20_mpc']:continue
    rows.append(f"| {r['z'][0]:.2f}–{r['z'][1]:.2f} | {r['common_parent_rows']:,} | "+' | '.join(f"{x['contrast_pp']:.2f} [{x['ci95_pp'][0]:.2f}, {x['ci95_pp'][1]:.2f}]" for x in r['arms'])+' |')
s=s.replace('<!-- TEST3_TABLE -->','\n'.join(rows))
# Persist the numerical tables in the readable source report too.
(G/'docs/environment_gap_three_tests_20261015.md').write_text(s)
body=md(s)
for token,path,caption in [('FIG_TEST1',D/'test1_ground_truth.png','True-field mass and scale experiment; 95% spatial intervals within one TNG volume.'),('FIG_TEST2',R/'loa_actual_gat_overlay_20261015/legacy_targetfit_soft_controlled_quenched.png','Actual GAT and Abacus on shared observed galaxies; see original report for each panel’s reference mixture.'),('FIG_TEST3',D/'test3_uncertainty.png','Posterior-width selection at fixed mass/redshift; selection weakens the observed association.')]:
    body=body.replace('<!-- '+token+' -->',fig(path,caption))
groups=[('Test 1 evidence',[(I/'docs/tng_density_control_20261015.md')],['tng_density_control_20261015']),('Test 2 evidence',[G/'docs/loa_actual_gat_overlay_20261015.md',I/'docs/legacy_gat_density_20261015.md',G/'docs/loa_tng_abacus_overlay_20261015.md',G/'docs/loa_selection_sensitivity_20261015.md'],['loa_actual_gat_overlay_20261015','loa_tng_abacus_overlay_20261015','loa_selection_sensitivity_20261015','legacy_gat_density_20261015']),('Test 3 evidence',[G/'docs/legacy_gat_comparison_20261014.md',R/'p12a_cg5_check_20261014/REPORT.md'],['legacy_gat_comparison_20261014','p12a_cg5_check_20261014']),('Additional test 4 evidence',[I/'docs/tracer_density_causality_20261015.md',G/'docs/raw_environment_checks_20261014.md',G/'docs/density_property_checks_20261014.md'],['abacus_density_followup_20261015','raw_environment_checks_20261014','density_checks_20261014']),('Complete chronological ledger and original numeric receipts',[G/'docs/legacy_gat_sequential_tests_20261014.md',R/'legacy_selection_tests_20261014/REPORT.md'],['legacy_selection_tests_20261014'])]
body+='<h2 id="archive">Complete source evidence and figure archive</h2><p>These sections preserve the original reports. Their chronological status is not silently rewritten: the synthesis above controls current conclusions. No blank or failed-support panel has been removed.</p>'
seen={str(x['path']) for x in assets}
for title,reports,folders in groups:
    body+='<details><summary>'+title+'</summary>'
    for p in reports:
        if not p.exists():continue
        assets.append(dict(path=str(p),sha256=sha(p),kind='source_report'))
        body+='<article><h3>'+html.escape(p.name)+'</h3><p class="source">Original source: '+html.escape(str(p))+'</p>'+md(p.read_text())+'</article>'
    for folder in folders:
        root=R/folder;body+='<h3>'+html.escape(folder)+'</h3>'
        pngs=[p for p in sorted(root.glob('*.png')) if not p.stem.endswith('_preview')]
        for p in pngs:
            if str(p) not in seen:body+=fig(p);seen.add(str(p))
        for p in sorted(root.glob('*.svg')):
            if p.with_suffix('.png').exists():continue
            assets.append(dict(path=str(p),sha256=sha(p),kind='figure'))
            body+='<figure><img loading="lazy" alt="'+html.escape(p.stem)+'" src="'+data(p,'image/svg+xml')+'"><figcaption>'+html.escape(p.stem)+'</figcaption></figure>'
        # Embed all PDFs: multi-page property panels are not lost when only a PNG summary exists.
        for p in sorted(root.glob('*.pdf')):body+=pdf(p)
        for name in ['RESULTS.json','VALIDATION.json','REPORT.md']:
            p=root/name
            if p.exists() and p not in reports:
                assets.append(dict(path=str(p),sha256=sha(p),kind='numeric_receipt' if p.suffix=='.json' else 'source_report'))
                body+='<details><summary>'+html.escape(folder+'/'+name)+'</summary><pre>'+html.escape(p.read_text())+'</pre></details>'
    body+='</details>'
body+='<h2>New machine-readable receipts</h2>'
for p in sorted(D.glob('*.json')):
    if p.name in ['MANIFEST.json','DOCUMENT_VALIDATION.json']:continue
    assets.append(dict(path=str(p),sha256=sha(p),kind='new_receipt'))
    body+='<details><summary>'+p.name+'</summary><pre>'+html.escape(p.read_text())+'</pre></details>'
body+='<h2>Embedded evidence inventory</h2>'+md('| Kind | Source | SHA-256 |\n|---|---|---|\n'+'\n'.join('| '+a['kind']+' | '+a['path'].replace(str(H)+'/','~/')+' | '+a['sha256']+' |' for a in assets))
style='''body{font:17px/1.6 system-ui,sans-serif;color:#20303e;background:#f6f8fa;margin:0}main{max-width:1120px;margin:auto;padding:3rem;background:white}h1{font-size:2.5rem;line-height:1.15}h2{border-top:2px solid #d9e5ed;padding-top:1.6rem;margin-top:2.5rem;color:#12647a}h3{color:#284e64}table{border-collapse:collapse;width:100%;font-size:.9rem;margin:1.5rem 0}th,td{padding:.55rem;border:1px solid #dbe3e9;text-align:left;overflow-wrap:anywhere}th{background:#edf4f7}tr:nth-child(even){background:#f8fafb}img{max-width:100%;height:auto}figure{margin:2rem 0}figcaption,.source{font-size:.85rem;color:#566876}details{border:1px solid #cfdee6;border-radius:6px;margin:1rem 0;padding:1rem}summary{cursor:pointer;font-weight:650}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:12px;max-height:40rem;overflow:auto}a{color:#006b92}nav{background:#eaf3f7;padding:1rem;border-radius:8px}article{border-bottom:3px solid #dae6eb;padding-bottom:2rem}@media print{body,main{background:white;padding:0}details{break-inside:auto}figure,table{break-inside:avoid}pre{max-height:none}nav{display:none}}'''
page='<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Environment relation: three diagnostic tests</title><style>'+style+'</style></head><body><main><nav><b>Offline evidence report</b> · Text, figures and original PDFs are embedded. Expand the evidence sections for full records. Use your browser’s Find to locate a test or filename. <a href="#archive">Jump to evidence archive</a>.</nav>'+body+'</main></body></html>'
p=O/'environment_gap_three_tests.html';p.write_text(page)
manifest=dict(document=str(p),document_sha256=sha(p),document_bytes=p.stat().st_size,assets=assets)
(D/'MANIFEST.json').write_text(json.dumps(manifest,indent=2))
assert '<!-- TEST' not in page and '<!-- FIG' not in page
assert page.count('<details>')==page.count('</details>')
assert page.count('<figure>')==page.count('</figure>')
(D/'DOCUMENT_VALIDATION.json').write_text(json.dumps(dict(html_complete=page.endswith('</html>'),embedded_figures=page.count('<figure>'),embedded_pdf_files=sum(a['kind']=='pdf' for a in assets),source_reports=sum(a['kind']=='source_report' for a in assets),asset_count=len(assets),all_images_verified=True,all_pdfs_header_eof_checked=True,all_assets_sha256_recorded=True,external_image_dependencies=0,rendered_browser_validation=False,limitation='No browser rendering tool available; image inspection, markup structure and asset validation performed.'),indent=2))
print(str(p),p.stat().st_size,'bytes',len(assets),'assets')
