"""Audit passive spray feedback and summarize independently solved ladder arms."""
import argparse
import json
from pathlib import Path
import numpy as np
from inspect_phase_repair import analyse, final
from spray_events_2d import load_frame


def passive(root):
    a,b=root/'escape_only',root/'full_r2';job=json.loads((a/'job.json').read_text())
    rows=[]
    for arm in ('python','bifrost'):
        maxerr={k:0. for k in ('pos','vel','rad')};checked=0;first=None
        for j in list(range(0,job['steps'],job['dump_every']))+['final']:
            pair=[final(p,arm) if j=='final' else load_frame(p,arm,j) for p in (a,b)]
            primary=[{k:z[~q['spr']] for k,z in q.items()} for q in pair]
            for k in ('pos','vel','typ','rad'):
                x,y=[q[k] for q in primary]
                if x.shape!=y.shape:
                    if first is None:first=dict(step=j,channel=k,shape=[x.shape,y.shape])
                    continue
                if k in maxerr:maxerr[k]=max(maxerr[k],float(abs(x.astype(float)-y.astype(float)).max()))
                if not np.array_equal(x,y) and first is None:first=dict(step=j,channel=k)
                checked+=1
        rows.append(dict(implementation=arm,bit_identical=first is None,first_difference=first,
                         primary_arrays_checked=checked,max_difference=maxerr))
    return rows


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args()
    arms={}
    for name in ('base','escape_only','full_r2','phase_escape_full'):
        r=analyse(a.root/name);arms[name]={k:v for k,v in r.items() if k not in ('rows','job')}
        if (a.root/name/'runtime.json').exists():arms[name]['runtime']=json.loads((a.root/name/'runtime.json').read_text())
        for impl in ('python','bifrost'):
            ev=json.loads((a.root/name/impl/'ledger.json').read_text())
            arms[name]['physics'][impl]['emitted_total']=sum(row.get('accepted',0) for row in ev)
    r=dict(arms=arms,passive_spray=passive(a.root))
    (a.root/'ladder_summary.json').write_text(json.dumps(r,indent=2))
    print(json.dumps(r,indent=2))
    if not all(v['bit_identical'] for v in r['passive_spray']):raise SystemExit('Passive-spray feedback check failed')
    if not all(v['similarity_pass'] for v in arms.values()):raise SystemExit('Independent sequence checks remain failed')


if __name__=='__main__':main()
