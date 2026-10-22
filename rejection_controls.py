import sys,json,copy,time
from pathlib import Path
root=Path(__file__).resolve().parent
sys.path.insert(0,str(root/'sources/roehrig/scripts'))
sys.path.insert(0,str(root/'independent'))
from certify_turan import certify,iv
from check_bounds import check_c42
from flint import ctx
iv.dps=100
ctx.prec=768
data=json.loads((root/'certificates/c42_256.json').read_text())
cases=[]
x=copy.deepcopy(data);x['C']='0.6888';cases.append(('too-small radius',x))
x=copy.deepcopy(data);x['eta']=[['0','0'] for _ in x['eta']];cases.append(('zero profile: analytic failure despite valid radii',x))
x=copy.deepcopy(data);x['breaks'][2]=x['breaks'][1];cases.append(('repeated threshold',x))
x=copy.deepcopy(data);x['C']='0.6888863709';cases.append(('near-threshold false bound',x))
results=[]
for name,x in cases:
    for engine,fn,terms in [('interval',certify,240),('Arb',check_c42,400)]:
        start=time.monotonic()
        try:fn(x,terms)
        except (AssertionError,ArithmeticError,ValueError) as exc:
            r={'case':name,'checker':engine,'result':'REJECTED','reason':str(exc),'seconds':round(time.monotonic()-start,3)}
            print(json.dumps(r),flush=True);results.append(r)
        else:raise RuntimeError('ACCEPTED INVALID CERTIFICATE: '+name+' '+engine)
print(json.dumps({'status':'PASS','controls':results},indent=2))
