"""Generate the Section 2 workflow schematic (Figure 1 of the manuscript):
the hierarchy from the input parameters theta = (d, r, b) to one network
realization (geometry -> intersection graph -> flow -> transport) and then to
the across-realization distributions the two surrogate models predict.

Run standalone:  `.venv/bin/python "Updated Code/figures/generate_workflow_figure.py"`
Writes fig_workflow.png and fig_workflow.pdf next to the other figures.

The DFN pipeline below is a Python port of the MATLAB routines in
"Updated Code/" -- generate_fractures.m, fracture_intersections.m,
build_network.m, backbone.m, solve_flow.m, flow_velocity.m and
transport_one.m -- reimplemented here so the whole figure builds in one
process. It is a faithful port of the algorithms and formulas, not a bit-exact
replica of any particular MATLAB seed: numpy and MATLAB use different PRNGs,
so this draws an independent realization from the same distributions, which is
what a schematic needs. Any change to the .m files should be mirrored here.
"""
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from matplotlib.gridspec import GridSpec
from matplotlib.collections import LineCollection

OUT = Path(__file__).resolve().parent

# ---------------------------------------------------------------------------
# 1. DFN pipeline (port of the MATLAB routines)
# ---------------------------------------------------------------------------
def generate_fractures(Lx, Ly, Nx, Ny, lx, ly, b, rng):
    xfmid = np.column_stack([rng.uniform(0, Lx, Nx), np.sort(rng.uniform(0, Ly, Nx))])
    yfmid = np.column_stack([np.sort(rng.uniform(0, Lx, Ny)), rng.uniform(0, Ly, Ny)])
    lfx = lx * rng.uniform(0, 1, Nx)
    lfy = ly * rng.uniform(0, 1, Ny)
    bfx = np.full(Nx, b); bfy = np.full(Ny, b)
    xfini = np.column_stack([xfmid[:,0]-0.5*lfx, xfmid[:,1]])
    xfend = np.column_stack([xfmid[:,0]+0.5*lfx, xfmid[:,1]])
    yfini = np.column_stack([yfmid[:,0], yfmid[:,1]-0.5*lfy])
    yfend = np.column_stack([yfmid[:,0], yfmid[:,1]+0.5*lfy])
    return xfmid, xfini, xfend, bfx, yfmid, yfini, yfend, bfy


def fracture_intersections(Lx, Ly, Nx, Ny, xfmid, xfini, xfend, yfmid, yfini, yfend):
    posnode, nfracx, nfracy = [], [], []
    for ix in range(Nx):
        xini, xend = xfini[ix,0], xfend[ix,0]
        ynode = xfmid[ix,1]
        if xini < 0:
            for iy in range(Ny):
                if (yfmid[iy,0] < xend) and (ynode > yfini[iy,1]) and (ynode < yfend[iy,1]):
                    posnode.append([0.0, ynode]); nfracx.append(ix+1); nfracy.append(0); break
        for iy in range(Ny):
            xn = yfmid[iy,0]
            if (xn > xini) and (xn < xend) and (ynode > yfini[iy,1]) and (ynode < yfend[iy,1]):
                posnode.append([xn, ynode]); nfracx.append(ix+1); nfracy.append(iy+1)
        if xend > Lx:
            for iy in range(Ny):
                if (yfmid[iy,0] > xini) and (ynode > yfini[iy,1]) and (ynode < yfend[iy,1]):
                    posnode.append([Lx, ynode]); nfracx.append(ix+1); nfracy.append(Ny+1); break
    return np.array(posnode).reshape(-1,2), np.array(nfracx,int), np.array(nfracy,int), len(posnode)


def build_network(nnode, posnode, nfracx, nfracy, Ny, bfx, bfy, visc):
    ln, la, lt = [], [], []
    for i in range(nnode-1):
        if nfracx[i+1] == nfracx[i]:
            dx = posnode[i+1,0]-posnode[i,0]
            ln.append([i,i+1]); la.append(bfx[nfracx[i]-1])
            lt.append(bfx[nfracx[i]-1]**3/(12*visc)/dx)
    ids = np.argsort(nfracy, kind='stable')
    for i in range(nnode-1):
        fa, fb = nfracy[ids[i]], nfracy[ids[i+1]]
        if fa == fb and 1 <= fa <= Ny:
            a, b_ = ids[i], ids[i+1]
            dy = posnode[b_,1]-posnode[a,1]
            ln.append([a,b_]); la.append(bfy[fa-1])
            lt.append(bfy[fa-1]**3/(12*visc)/dy)
    return np.array(ln).reshape(-1,2), np.array(la), np.array(lt), len(ln)


def backbone(linknods, linktrans, posnode, nnode, Lx):
    if len(linknods) == 0:
        return np.array([],int), np.arange(nnode), np.array([],int), np.array([],int)
    A = csr_matrix((np.ones(len(linknods)), (linknods[:,0], linknods[:,1])), shape=(nnode,nnode))
    A = A + A.T
    ncomp, lab = connected_components(A, directed=False)
    inl = np.where(posnode[:,0] == 0)[0]
    out = np.where(posnode[:,0] == Lx)[0]
    good = set(lab[inl]) & set(lab[out])
    bb = np.where(np.isin(lab, list(good)))[0]
    iso = np.setdiff1d(np.arange(nnode), bb)
    return bb, iso, inl, out


def solve_flow(linknods, linktrans, posnode, nnode, bb, inl, out, pin, pout):
    """Nodal pressures on the backbone with Dirichlet inlet/outlet."""
    pos = {n:k for k,n in enumerate(bb)}
    m = len(bb)
    L = np.zeros((m,m))
    for (u,v), T in zip(linknods, linktrans):
        if u in pos and v in pos:
            iu, iv = pos[u], pos[v]
            L[iu,iu] += T; L[iv,iv] += T; L[iu,iv] -= T; L[iv,iu] -= T
    fixed_n = np.intersect1d(np.union1d(inl,out), bb)
    fixed = np.array([pos[n] for n in fixed_n])
    pfix = np.array([pin if posnode[n,0]==0 else pout for n in fixed_n])
    free = np.setdiff1d(np.arange(m), fixed)
    p = np.zeros(m); p[fixed] = pfix
    if len(free):
        p[free] = np.linalg.solve(L[np.ix_(free,free)], -L[np.ix_(free,fixed)] @ pfix)
    Q = np.zeros(len(linknods))
    for k,((u,v),T) in enumerate(zip(linknods, linktrans)):
        if u in pos and v in pos:
            Q[k] = T*(p[pos[u]] - p[pos[v]])
    return p, Q


def transport_one(linknods, Q, linkapert, posnode, inl_bb, out_set, rng, maxstep=10000):
    inlinks = [k for k,(u,v) in enumerate(linknods) if u in inl_bb and Q[k] > 0]
    if not inlinks: return None, None
    w = np.array([Q[k] for k in inlinks]); w = w/w.sum()
    k0 = inlinks[rng.choice(len(inlinks), p=w)]
    node = linknods[k0,0]
    path=[node]; t=0.0; times=[0.0]
    for _ in range(maxstep):
        cand=[]
        for k,(u,v) in enumerate(linknods):
            if u==node and Q[k]>0: cand.append((k,v,Q[k]))
            elif v==node and Q[k]<0: cand.append((k,u,-Q[k]))
        if not cand: return None, None
        ws=np.array([c[2] for c in cand]); ws=ws/ws.sum()
        k,nxt,q = cand[rng.choice(len(cand), p=ws)]
        ln = abs(posnode[nxt,0]-posnode[node,0])+abs(posnode[nxt,1]-posnode[node,1])
        vel = q/linkapert[k]
        t += ln/abs(vel)
        node = nxt; path.append(node); times.append(t)
        if node in out_set: return np.array(path), np.array(times)
    return None, None


def realize(Lx, Ly, d, r, b, visc, pin, pout, rng):
    Nx, Ny = int(round(Lx/d)), int(round(Ly/d))
    lx, ly = r*Lx, r*Ly
    xfmid,xfini,xfend,bfx,yfmid,yfini,yfend,bfy = generate_fractures(Lx,Ly,Nx,Ny,lx,ly,b,rng)
    posnode,nfx,nfy,nn = fracture_intersections(Lx,Ly,Nx,Ny,xfmid,xfini,xfend,yfmid,yfini,yfend)
    D = dict(xfini=xfini,xfend=xfend,yfini=yfini,yfend=yfend,posnode=posnode,nnode=nn,connected=False)
    if nn == 0: return D
    ln,la,lt,nl = build_network(nn,posnode,nfx,nfy,Ny,bfx,bfy,visc)
    D.update(linknods=ln,linkapert=la,linktrans=lt,nlink=nl)
    if nl <= 1: return D
    bb,iso,inl,out = backbone(ln,lt,posnode,nn,Lx)
    D.update(bb=bb,iso=iso,inl=inl,out=out)
    if len(bb)==0: return D
    p,Q = solve_flow(ln,lt,posnode,nn,bb,inl,out,pin,pout)
    inl_bb=set(np.intersect1d(inl,bb).tolist()); out_bb=set(np.intersect1d(out,bb).tolist())
    Qtot = sum(Q[k] for k,(u,v) in enumerate(ln) if u in inl_bb)
    D.update(p=p,Q=Q,Qtot=Qtot,inl_bb=inl_bb,out_bb=out_bb,connected=Qtot>0)
    return D

# ---------------------------------------------------------------------------
# 2. Ensemble at a fixed theta
# ---------------------------------------------------------------------------
Lx, Ly, visc, pin, pout = 150., 50., 1e-3, 1e5, 0.0
D_, R_, B_ = 2.0, 0.50, 1.5e-4     # theta used for the schematic
NREAL, NPART = 1000, 200
QGRID = np.array([0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.15,0.2,0.3,
                  0.4,0.5,0.6,0.7,0.8,0.9,0.95])


def fast_paths(D, rng, npart):
    ln,Q,la,pos = D['linknods'],D['Q'],D['linkapert'],D['posnode']
    out_bb = D['out_bb']
    adj = {}
    for k,(u,v) in enumerate(ln):
        if Q[k] > 0:  adj.setdefault(u,[]).append((k,v, Q[k]))
        elif Q[k] < 0: adj.setdefault(v,[]).append((k,u,-Q[k]))
    inlinks = [(k,ln[k,0]) for k,(u,v) in enumerate(ln) if u in D['inl_bb'] and Q[k] > 0]
    if not inlinks: return [],[]
    w = np.array([Q[k] for k,_ in inlinks]); w = w/w.sum()
    paths, times = [], []
    for _ in range(npart):
        node = inlinks[rng.choice(len(inlinks), p=w)][1]
        p=[node]; t=0.0
        for _ in range(5000):
            c = adj.get(node)
            if not c: t=None; break
            ws = np.array([x[2] for x in c]); ws=ws/ws.sum()
            k,nxt,q = c[rng.choice(len(c), p=ws)]
            t += (abs(pos[nxt,0]-pos[node,0])+abs(pos[nxt,1]-pos[node,1]))/abs(q/la[k])
            node = nxt; p.append(node)
            if node in out_bb: break
        else: t=None
        if t is not None: paths.append(np.array(p)); times.append(t)
    return paths, times


def build_ensemble():
    """Run NREAL realizations at the schematic's theta; return the ensemble
    statistics, one legible realization with its particle paths, and nine
    thumbnails for the connected/disconnected panel."""
    t0 = time.time()
    Qtots, tq_rows, nconn = [], [], 0
    for i in range(NREAL):
        rng = np.random.default_rng(10000 + i)
        D = realize(Lx, Ly, D_, R_, B_, visc, pin, pout, rng)
        if not D['connected']:
            Qtots.append(0.0)
            continue
        nconn += 1
        Qtots.append(D['Qtot'])
        if len(tq_rows) < 250:
            _, ts = fast_paths(D, rng, NPART)
            if len(ts) >= 20:
                ts = np.sort(ts)
                tq_rows.append([ts[min(int(np.floor(q * len(ts))), len(ts) - 1)] for q in QGRID])
    Qtots = np.array(Qtots)
    tq = np.array(tq_rows)
    p0 = 1 - nconn / NREAL
    print(f'{time.time()-t0:.0f}s  connected {nconn}/{NREAL}  p0={p0:.3f}  tq rows={len(tq)}')

    # a legible single realization for panels (a)-(d)
    best = None
    for seed in range(400):
        Dr = realize(Lx, Ly, D_, R_, B_, visc, pin, pout, np.random.default_rng(seed))
        if not Dr['connected']:
            continue
        nb = len(Dr['bb'])
        if 65 <= nb <= 95 and len(Dr['inl_bb']) >= 2 and len(Dr['out_bb']) >= 2:
            score = abs(nb - 80)
            if best is None or score < best[0]:
                best = (score, seed, Dr)
    _, seed, Dr = best
    paths, ts = fast_paths(Dr, np.random.default_rng(seed + 1), 400)
    print(f'realization seed={seed} nodes={Dr["nnode"]} links={Dr["nlink"]} '
          f'bb={len(Dr["bb"])} Q={Dr["Qtot"]:.3e} paths={len(paths)}')

    # nine thumbnails for panel (e)
    thumbs = []
    for t in range(9):
        Dt = realize(Lx, Ly, D_, R_, B_, visc, pin, pout, np.random.default_rng(500 + t))
        thumbs.append(dict(xfini=Dt['xfini'], xfend=Dt['xfend'], yfini=Dt['yfini'],
                           yfend=Dt['yfend'], connected=Dt['connected']))

    E = dict(Qtots=Qtots, tq=tq, QGRID=QGRID, p0=p0,
             theta=np.array([D_, R_, B_]), nconn=nconn, nreal=NREAL)
    return E, (Dr, paths, ts, seed), thumbs


# ---------------------------------------------------------------------------
# 3. Figure
# ---------------------------------------------------------------------------
E, (R_real, paths_, ts_, seed_), TH = build_ensemble()
D = dict(D=R_real, paths=paths_[:250], ts=ts_[:250], seed=seed_)
BLUE, ORANGE, INK, MUTE = '#2E6FA7', '#C8622D', '#22303C', '#B9C2C9'
plt.rcParams.update({'font.size':12,'axes.linewidth':0.9,'xtick.major.width':0.7,
                     'ytick.major.width':0.7,'axes.edgecolor':'#5A6570',
                     'xtick.color':'#5A6570','ytick.color':'#5A6570',
                     'axes.labelcolor':INK,'text.color':INK})

R, paths, ts = D['D'], D['paths'], D['ts']
d,r,b = E['theta']; p0 = float(E['p0']); Lx,Ly = 150.,50.
pos, ln, Q, bb = R['posnode'], R['linknods'], R['Q'], R['bb']
bbs, iso = set(bb.tolist()), set(R['iso'].tolist())

fig = plt.figure(figsize=(16.5, 8.6))
gs = GridSpec(2, 4, figure=fig, hspace=0.42, wspace=0.30, height_ratios=[0.80, 1.0],
              left=0.052, right=0.985, top=0.885, bottom=0.085)

def dom(ax, ticks=False):
    ax.set_xlim(-3, Lx+3); ax.set_ylim(-5, 94); ax.set_aspect('equal')
    ax.add_patch(plt.Rectangle((0,0),Lx,Ly,fill=False,ec='#5A6570',lw=0.8,zorder=1))
    if not ticks: ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)

def panel(ax, tag, title, dx=-0.015):
    ax.set_title(f'  {title}', fontsize=14, loc='left', pad=8, color=INK)
    ax.text(dx, 1.055, tag, transform=ax.transAxes, fontsize=16,
            fontweight='bold', va='top', ha='right', color=INK)

# (a) fractures ------------------------------------------------------------
ROW1=[]
ax = fig.add_subplot(gs[0,0]); dom(ax); ROW1.append(ax)
ax.add_collection(LineCollection(
    [[(R['xfini'][i,0],R['xfini'][i,1]),(R['xfend'][i,0],R['xfend'][i,1])] for i in range(len(R['xfini']))],
    colors=BLUE, lw=1.0, zorder=2))
ax.add_collection(LineCollection(
    [[(R['yfini'][i,0],R['yfini'][i,1]),(R['yfend'][i,0],R['yfend'][i,1])] for i in range(len(R['yfini']))],
    colors=ORANGE, lw=1.0, zorder=2))
panel(ax,'a','Fracture realization')
ax.text(2,88,r'$x$-family, $N_x$ fractures',color=BLUE,fontsize=11.5,va='center')
ax.text(2,73,r'$y$-family, $N_y$ fractures',color=ORANGE,fontsize=11.5,va='center')
ax.text(2,58,r'$L_x \times L_y = 150 \times 50$ m',fontsize=11.5,color='#5A6570',va='center')

# (b) intersection graph ---------------------------------------------------
ax = fig.add_subplot(gs[0,1]); dom(ax); ROW1.append(ax)
seg_bb = [[(pos[u,0],pos[u,1]),(pos[v,0],pos[v,1])] for u,v in ln if u in bbs and v in bbs]
seg_is = [[(pos[u,0],pos[u,1]),(pos[v,0],pos[v,1])] for u,v in ln if not(u in bbs and v in bbs)]
ax.add_collection(LineCollection(seg_is, colors=MUTE, lw=0.8, zorder=2))
ax.add_collection(LineCollection(seg_bb, colors=INK, lw=1.0, zorder=3))
ii = np.array(sorted(iso), int)
if len(ii): ax.scatter(pos[ii,0],pos[ii,1],s=7,c=MUTE,zorder=4,lw=0)
ax.scatter(pos[bb,0],pos[bb,1],s=9,c=INK,zorder=5,lw=0)
inl = np.array(sorted(R['inl_bb']),int); out = np.array(sorted(R['out_bb']),int)
ax.scatter(pos[inl,0],pos[inl,1],s=42,facecolors='none',edgecolors=BLUE,lw=1.4,zorder=6)
ax.scatter(pos[out,0],pos[out,1],s=42,facecolors='none',edgecolors=ORANGE,lw=1.4,zorder=6)
panel(ax,'b','Intersection graph')
ax.text(2,88,r'%d nodes,  %d links,  %d in backbone'%(R['nnode'],R['nlink'],len(bb)),
        color='#5A6570',fontsize=11.5,va='center')
ax.text(2,73,'backbone',color=INK,fontsize=11.5,va='center')
ax.text(62,73,'isolated',color='#8A949C',fontsize=11.5,va='center')
ax.text(2,58,'inlet nodes',color=BLUE,fontsize=11.5,va='center')
ax.text(62,58,'outlet nodes',color=ORANGE,fontsize=11.5,va='center')

# (c) flow ------------------------------------------------------------------
ax = fig.add_subplot(gs[0,2]); dom(ax); ROW1.append(ax)
aq = np.abs(Q); mx = aq[aq>0].max()
segs, lws = [], []
for k,(u,v) in enumerate(ln):
    if u in bbs and v in bbs and aq[k] > 0:
        segs.append([(pos[u,0],pos[u,1]),(pos[v,0],pos[v,1])])
        lws.append(0.35 + 3.1*(aq[k]/mx)**0.45)
ax.add_collection(LineCollection(segs, colors='#7E8992', lw=lws, zorder=2))
sc = ax.scatter(pos[bb,0],pos[bb,1],s=17,c=R['p'][:len(bb)]/1e5,cmap='Blues',
                vmin=0,vmax=1,zorder=4,lw=0.25,edgecolors='#5A6570')
cb = fig.colorbar(sc, ax=ax, fraction=0.028, pad=0.015); CB_C = cb.ax
cb.set_label(r'$p/p_{\rm in}$', fontsize=12); cb.ax.tick_params(labelsize=10.5); cb.outline.set_visible(False)
panel(ax,'c','Pressure and flow')
ax.text(2,88,r'line width $\propto |Q_e|$',color='#5A6570',fontsize=11.5,va='center')
ax.text(2,73,r'$Q=%.1f\times10^{-7}$ m$^2$/s'%(R['Qtot']*1e7),
        color=INK,fontsize=11.5,va='center')

# (d) transport -------------------------------------------------------------
ax = fig.add_subplot(gs[0,3]); dom(ax); ROW1.append(ax)
ax.add_collection(LineCollection(seg_bb, colors='#DCE1E5', lw=0.8, zorder=2))
for p_ in paths[:250]:
    ax.plot(pos[p_,0], pos[p_,1], color=ORANGE, lw=0.8, alpha=0.10, zorder=3, solid_capstyle='round')
ax.scatter(pos[inl,0],pos[inl,1],s=26,facecolors='none',edgecolors=BLUE,lw=1.1,zorder=6)
ax.scatter(pos[out,0],pos[out,1],s=26,facecolors='none',edgecolors=ORANGE,lw=1.1,zorder=6)
panel(ax,'d','Advective particle tracking')
axi = ax.inset_axes([0.62,0.58,0.36,0.38])
axi.hist(np.log10(ts), bins=22, color=ORANGE, alpha=0.85, lw=0)
axi.set_yticks([]); axi.tick_params(labelsize=9.5, length=3)
axi.set_xlabel(r'$\log_{10} t_{\rm bt}$  (one realization)', fontsize=10, labelpad=2)
for s in ('top','right','left'): axi.spines[s].set_visible(False)
axi.patch.set_alpha(0.85)
ax.text(2,88,r'$N_{\rm particle}$ particles,',color='#5A6570',fontsize=11.5,va='center')
ax.text(2,73,r'flux-weighted',color='#5A6570',fontsize=11.5,va='center')

# (e) ensemble of realizations ---------------------------------------------
ax = fig.add_subplot(gs[1,0]); ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values(): s.set_visible(False)
ax.set_xlim(0,3); ax.set_ylim(0,3)
for i,t in enumerate(TH):
    cx, cy = i%3, 2-i//3
    conn = t['connected']; col = INK if conn else MUTE
    segs=[[(t['xfini'][j,0],t['xfini'][j,1]),(t['xfend'][j,0],t['xfend'][j,1])] for j in range(len(t['xfini']))]
    segs+=[[(t['yfini'][j,0],t['yfini'][j,1]),(t['yfend'][j,0],t['yfend'][j,1])] for j in range(len(t['yfini']))]
    segs=[[(cx+0.06+0.88*x/Lx, cy+0.30+0.40*y/Ly) for x,y in s_] for s_ in segs]
    ax.add_collection(LineCollection(segs, colors=col, lw=0.35, zorder=2))
    ax.add_patch(plt.Rectangle((cx+0.06,cy+0.30),0.88,0.40,fill=False,
                 ec=BLUE if conn else MUTE, lw=1.0 if conn else 0.7, zorder=3))
    ax.text(cx+0.5, cy+0.16, 'connected' if conn else 'disconnected', ha='center',
            fontsize=9.5, color=BLUE if conn else '#8A949C')
panel(ax,'e',r'Ensemble at fixed $\theta$')
ax.text(0.015,-0.075,r'$N_{\rm network}$ realizations $\rightarrow$ $p_0 = %.2f$'%p0,
        transform=ax.transAxes,fontsize=11.5,color='#5A6570')

# (f) flow-rate distribution -------------------------------------------------
ax = fig.add_subplot(gs[1,1])
lq = np.log10(E['Qtots'][E['Qtots']>0])
ax.hist(lq, bins=34, color=BLUE, alpha=0.85, lw=0)
ax.set_xlabel(r'$\log_{10} Q$   (m$^2$/s)'); ax.set_ylabel('realizations')
for s in ('top','right'): ax.spines[s].set_visible(False)
panel(ax,'f','Total-flow-rate distribution',dx=-0.16)
ax.text(0.015,-0.235,r'target of the flow surrogate',transform=ax.transAxes,
        fontsize=11.5,color='#5A6570')

# (g) percentile-time curves --------------------------------------------------
ax = fig.add_subplot(gs[1,2])
qg, tq = E['QGRID'], E['tq']
for row in tq[:160]:
    ax.plot(np.log10(row), qg, color=ORANGE, lw=0.5, alpha=0.16)
ax.plot(np.log10(np.median(tq,axis=0)), qg, color=INK, lw=1.8, label='median realization')
ax.set_xlabel(r'$\log_{10} t$   (s)'); ax.set_ylabel(r'percentile $q$')
ax.set_ylim(0,1); ax.legend(frameon=False, fontsize=11.5, loc='lower right')
for s in ('top','right'): ax.spines[s].set_visible(False)
panel(ax,'g',r'Breakthrough times $t_q$',dx=-0.15)

# (h) empirical F(q, log t) ---------------------------------------------------
ax = fig.add_subplot(gs[1,3])
lt = np.log10(tq); grid = np.linspace(lt.min(), lt.max(), 170)
Fs = np.array([[np.mean(lt[:,j] <= g) for g in grid] for j in range(len(qg))])
im = ax.pcolormesh(grid, qg, Fs, cmap='viridis', vmin=0, vmax=1, shading='auto', rasterized=True)
cs = ax.contour(grid, qg, Fs, levels=[0.1,0.5,0.9], colors='white', linewidths=0.8, alpha=0.85)
ax.clabel(cs, fmt='%.1f', fontsize=9.5, colors='white')
cb = fig.colorbar(im, ax=ax, fraction=0.036, pad=0.015)
cb.set_label(r'$F(\theta,q,t)$', fontsize=12); cb.ax.tick_params(labelsize=10.5); cb.outline.set_visible(False)
ax.set_xlabel(r'$\log_{10} t$   (s)'); ax.set_ylabel(r'percentile $q$')
for s in ax.spines.values(): s.set_visible(False)
panel(ax,'h','Breakthrough-time surface',dx=-0.15)
ax.text(0.015,-0.235,r'target of the transport surrogate',transform=ax.transAxes,
        fontsize=11.5,color='#5A6570')

# headers / flow arrows -------------------------------------------------------
fig.text(0.052,0.975,r'Input parameters  $\theta = (d,\,r,\,b) = (%.1f\ \mathrm{m},\ %.2f,\ %.2f\times10^{-4}\ \mathrm{m})$'%(d,r,b*1e4),
         fontsize=15.5, color=INK, va='top')
fig.text(0.052,0.545,'Across $N_{\\rm network}$ realizations at the same $\\theta$:  the distributions the surrogates are trained to predict',
         fontsize=12.5, color='#5A6570', style='italic', va='top')
fig.canvas.draw()
bb1 = [a_.get_position() for a_ in ROW1]
ymid = bb1[0].y0 + 0.34*bb1[0].height
for i in range(3):
    x0 = CB_C.get_position().x1 if i == 2 else bb1[i].x1
    x1 = bb1[i+1].x0
    fig.text(0.5*(x0+x1), ymid, r'$\rightarrow$', fontsize=21, color='#8A949C',
             ha='center', va='center')
fig.text(0.016, 0.50, r'repeat $\times\, N_{\rm network}$', fontsize=11.5, color='#8A949C',
         rotation=90, va='center', ha='center')
fig.patches.append(matplotlib.patches.FancyArrowPatch(
    (0.036, 0.62), (0.036, 0.40), transform=fig.transFigure, figure=fig,
    arrowstyle='-|>', mutation_scale=11, color='#8A949C', lw=0.9))
fig.savefig(OUT / 'fig_workflow.png', dpi=220, bbox_inches='tight', facecolor='white')
fig.savefig(OUT / 'fig_workflow.pdf', bbox_inches='tight', facecolor='white')
print(f'saved {OUT / "fig_workflow.png"} and .pdf')
