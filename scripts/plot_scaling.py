"""W2 scaling plot: median steps-to-solution vs n (Gaussian constraint bits), S3 annealed-lambda variants.
Censored medians (fewer than half solved at the cap) drawn hollow at the cap."""
import json, collections, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
nb={}
for l in open('logs/nbits_ladder.jsonl'): x=json.loads(l); nb[(x['v'],x['d'])]=x['n_bits']
runs=collections.defaultdict(list)
for f in ['logs/grid1.jsonl','logs/grid3.jsonl']:
    try:
        for l in open(f):
            r=json.loads(l)
            if r['kind']=='s3sched': runs[(r['lam_min'],r['v'],r['d'])].append(r)
    except FileNotFoundError: pass
colors={0.0:'#2a78d6',0.3:'#eb6834',1.0:'#1baf7a'}
fig,ax=plt.subplots(figsize=(7,4.5),dpi=150)
for lm,col in colors.items():
    xs=[];ys=[];cx=[];cy=[]
    for (l,v,d),rs in sorted(runs.items(),key=lambda t:nb[(t[0][1],t[0][2])]):
        if l!=lm: continue
        st=sorted(r['steps'] if r['solved'] else np.inf for r in rs)
        med=st[(len(st)-1)//2]
        if np.isfinite(med): xs.append(nb[(v,d)]); ys.append(med)
        else: cx.append(nb[(v,d)]); cy.append(max(r['steps'] for r in rs))
    ax.plot(xs,ys,'-o',color=col,lw=2,ms=7,label=f'S3, λ floor {lm}')
    if cx: ax.plot(cx,cy,'^',color=col,ms=9,mfc='white',mew=2)
ax.set_yscale('log',base=2); ax.set_xlabel('n (constraint bits, Gaussian estimate)'); ax.set_ylabel('median flips to solution')
ax.set_title('S3 tabu: cost vs n on the b=0 ladder (hollow ▲ = unsolved at cap)',fontsize=10)
ax.grid(alpha=0.25); ax.legend(frameon=False,fontsize=9)
for s in ['top','right']: ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig('report/w2_scaling.png'); print('ok')
