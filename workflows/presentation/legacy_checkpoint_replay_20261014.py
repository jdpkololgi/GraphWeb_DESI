"""Identify old GAT weights using eight galaxies and exact four-hop context."""
import ast
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.nn import functional as F
from torch_geometric.nn import GATv2Conv
from torch_geometric.utils import k_hop_subgraph

ROOT=Path('/global/homes/d/dkololgi/GraphWeb_DESI')
OUT=Path('/global/homes/d/dkololgi/fmsu_tokyo_20261010/legacy_selection_tests_20261014')
ARCHIVE=Path('/pscratch/sd/d/dkololgi/tng_illustris/cache/from_repo_root')


def main():
    torch.set_num_threads(1)
    # Load only the canonical model class, avoiding catalogue/graph construction.
    source=ROOT/'workflows/graph_inference/graph_catalog.py';tree=ast.parse(source.read_text())
    node=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='SimpleGAT')
    namespace=dict(torch=torch,nn=nn,F=F,GATv2Conv=GATv2Conv)
    exec(compile(ast.Module(body=[node],type_ignores=[]),str(source),'exec'),namespace)
    model=namespace['SimpleGAT'](10,4,num_heads=4).eval()
    graph=torch.load(ROOT/'cache/DESI_alpha_geom.pt',map_location='cpu',weights_only=False)
    catalogue=pd.read_pickle(ROOT/'cache/DESI_NETWORKalpha_zcat.pt')
    vac=pd.read_pickle(ROOT/'DESI_BGS_PRERELEASE_VAC.pkl')
    assert np.array_equal(catalogue.TARGETID,vac.TARGETID)
    assert np.array_equal(catalogue.Z,vac.Z)
    selected=[]
    for c in range(4):
        rows=np.flatnonzero(vac.GAT_ENV.to_numpy()==c);selected.extend(rows[[len(rows)//3,2*len(rows)//3]].tolist())
    selected=np.sort(selected);target=torch.as_tensor(selected)
    subset,edges,mapping,emask=k_hop_subgraph(target,4,graph.edge_index,relabel_nodes=True,num_nodes=graph.num_nodes)
    assert len(subset)<100000 and edges.shape[1]<1500000
    x=graph.x[subset];attr=graph.edge_attr[emask]
    observed=vac[['GAT_'+s+'_PROB' for s in ['VOID','WALL','FILAMENT','CLUSTER']]].to_numpy()[selected]
    results=[]
    with torch.inference_mode():
        for path in sorted(ARCHIVE.glob('*gat*.pth')):
            state=torch.load(path,map_location='cpu',weights_only=False)
            model.load_state_dict(state)
            pred=torch.softmax(model(x,edges,attr)[mapping],dim=1).numpy()
            results.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                max_abs_probability_error=float(np.max(abs(pred-observed))),
                                mean_abs_probability_error=float(np.mean(abs(pred-observed))),
                                labels_matching=int(np.sum(pred.argmax(1)==observed.argmax(1)))))
    results.sort(key=lambda r:r['max_abs_probability_error'])
    best=results[0];model.load_state_dict(torch.load(best['path'],map_location='cpu',weights_only=False))
    # One additional hop must leave these targets unchanged; tests context completeness.
    subset5,edges5,mapping5,emask5=k_hop_subgraph(target,5,graph.edge_index,relabel_nodes=True,num_nodes=graph.num_nodes)
    assert len(subset5)<150000 and edges5.shape[1]<2000000
    with torch.inference_mode():
        p4=torch.softmax(model(x,edges,attr)[mapping],dim=1).numpy()
        p5=torch.softmax(model(graph.x[subset5],edges5,graph.edge_attr[emask5])[mapping5],dim=1).numpy()
    parity=float(np.max(abs(p4-p5)));assert parity<2e-6
    result=dict(sampled_targetids=vac.TARGETID.to_numpy()[selected].tolist(),four_hop_nodes=len(subset),
                four_hop_edges=edges.shape[1],five_hop_max_difference=parity,checkpoints=results,
                source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                identity='Cached catalogue IDs and redshifts exactly match saved VAC. Model inference only for eight diagnostic galaxies; no model comparison or retraining.')
    OUT.mkdir(exist_ok=True);(OUT/'checkpoint_replay.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
