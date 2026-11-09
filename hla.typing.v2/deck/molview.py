"""Pseudo-3D backbone drawings of crystal structures for the deck (matplotlib only).

The C-alpha trace of each chain is smoothed with a cubic spline and drawn as a tube: helices wide, strands
medium, loops thin; segments are painted back to front with a dark rim, and far segments are faded towards
white (depth cueing). Secondary structure comes from the HELIX/SHEET records of the PDB file.
"""
from pathlib import Path

import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgb
from scipy.interpolate import CubicSpline


def read_pdb(path, chains):
    """C-alpha coordinates and secondary structure of the named chains (first model, altloc A or blank)."""
    ca = {c: {} for c in chains}
    ss = {c: {} for c in chains}
    for ln in Path(path).read_text().splitlines():
        rec = ln[:6].strip()
        if rec == 'ENDMDL':
            break
        if rec == 'ATOM' and ln[12:16].strip() == 'CA' and ln[16] in ' A' and ln[21] in chains:
            ca[ln[21]].setdefault(int(ln[22:26]), np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]))
        elif rec == 'HELIX' and ln[19] in chains:
            for r in range(int(ln[21:25]), int(ln[33:37]) + 1):
                ss[ln[19]][r] = 'H'
        elif rec == 'SHEET' and ln[21] in chains:
            for r in range(int(ln[22:26]), int(ln[33:37]) + 1):
                ss[ln[21]].setdefault(r, 'E')
    out = {}
    for c in chains:
        res = np.array(sorted(ca[c]))
        out[c] = {'res': res, 'xyz': np.array([ca[c][r] for r in res]),
                  'ss': np.array([ss[c].get(r, 'C') for r in res])}
    return out


def centroid(mol, sel):
    """Mean C-alpha position over [(chain, first, last), ...]."""
    pts = [mol[c]['xyz'][(mol[c]['res'] >= a) & (mol[c]['res'] <= b)] for c, a, b in sel]
    return np.vstack(pts).mean(axis=0)


def direction(mol, chain, a, b):
    m = mol[chain]
    return m['xyz'][m['res'] == b][0] - m['xyz'][m['res'] == a][0]


def frame(up, along):
    """Orthonormal frame: e_up, e_along (orthogonalised), e_depth = along x up."""
    u = up / np.linalg.norm(up)
    x = along - along.dot(u) * u
    x /= np.linalg.norm(x)
    z = np.cross(x, u)
    return u, x, z


def trace(mol, chain, n=8):
    """Smoothed trace: points, residue number and secondary structure of each point; split at chain breaks."""
    m = mol[chain]
    gaps = np.where(np.linalg.norm(np.diff(m['xyz'], axis=0), axis=1) > 4.3)[0]
    pieces = np.split(np.arange(len(m['res'])), gaps + 1)
    out = []
    for idx in pieces:
        if len(idx) < 2:
            continue
        t = np.arange(len(idx))
        tt = np.linspace(0, len(idx) - 1, (len(idx) - 1) * n + 1)
        sp = CubicSpline(t, m['xyz'][idx], axis=0)
        near = np.clip(np.round(tt).astype(int), 0, len(idx) - 1)
        out.append((sp(tt), m['res'][idx][near], m['ss'][idx][near]))
    return out


WIDTH = {'H': 1.0, 'E': 0.72, 'C': 0.36}


def draw(ax, mol, view, colour, width=6.0, chains=None, keep=None, fade=0.35, rim=0.5, layers=24, z_order=2,
         shift=(0, 0)):
    """Draw the chains on ax in the view (e_x, e_y, e_depth); colour(chain, res) -> colour; keep(chain, res)
    -> bool selects residues. Segments are grouped into depth layers painted back to front; within a layer
    every rim is drawn before every face, so a tube shows no seams. Returns the projected C-alpha positions
    {chain: {res: (x, y)}} for labelling."""
    ex, ey, ez = view
    runs = []                                                   # (depth, layer key, xy, colour, width)
    proj = {}
    for c in chains or mol:
        m = mol[c]
        proj[c] = {r: (p.dot(ex) + shift[0], p.dot(ey) + shift[1]) for r, p in zip(m['res'], m['xyz'])}
        for pts, res, ss in trace(mol, c):
            x, y, z = pts @ ex + shift[0], pts @ ey + shift[1], pts @ ez
            for k in range(len(x) - 1):
                if keep is not None and not keep(c, res[k]):
                    continue
                runs.append(((z[k] + z[k + 1]) / 2, (x[k], y[k]), (x[k + 1], y[k + 1]),
                             to_rgb(colour(c, res[k])), width * WIDTH[ss[k]]))
    if not runs:
        return proj
    depth = np.array([r[0] for r in runs])
    lo, hi = depth.min(), depth.max()
    lay = np.minimum(((depth - lo) / ((hi - lo) or 1) * layers).astype(int), layers - 1)
    for L in range(layers):
        idx = np.where(lay == L)[0]
        if not len(idx):
            continue
        dn = (L + 0.5) / layers                                 # 0 far, 1 near
        segs = [[runs[i][1], runs[i][2]] for i in idx]
        base = np.array([runs[i][3] for i in idx])
        w = np.array([runs[i][4] for i in idx])
        face = base * (1 - fade * (1 - dn)) + fade * (1 - dn)
        edge = face * (1 - rim)
        ax.add_collection(LineCollection(segs, colors=edge, linewidths=w + 1.5, capstyle='round',
                                         joinstyle='round', zorder=z_order + 2 * L / layers))
        ax.add_collection(LineCollection(segs, colors=face, linewidths=w, capstyle='round',
                                         joinstyle='round', zorder=z_order + (2 * L + 1) / layers))
    return proj


def tilt(view, deg):
    """Rotate the view about its screen x axis: positive tips the top of the molecule towards the viewer."""
    ex, ey, ez = view
    a = np.radians(deg)
    return ex, ey * np.cos(a) + ez * np.sin(a), -ey * np.sin(a) + ez * np.cos(a)


def spin(view, deg):
    """Rotate the view about its screen y axis."""
    ex, ey, ez = view
    a = np.radians(deg)
    return ex * np.cos(a) + ez * np.sin(a), ey, -ex * np.sin(a) + ez * np.cos(a)


def peptide(ax, mol, chain, view, colour, size=60, z_order=10, shift=(0, 0)):
    """Peptide as balls (C-alpha) and sticks, on top of the drawing."""
    ex, ey, _ = view
    p = mol[chain]['xyz']
    x, y = p @ ex + shift[0], p @ ey + shift[1]
    ax.plot(x, y, color=colour, lw=2.4, zorder=z_order, solid_capstyle='round')
    ax.scatter(x, y, s=size, color=colour, edgecolor='white', lw=0.8, zorder=z_order + 0.1)
    return x, y


def extent(mol, view, chains=None, shift=(0, 0)):
    ex, ey, _ = view
    pts = np.vstack([mol[c]['xyz'] for c in (chains or mol)])
    return (pts @ ex + shift[0]).min(), (pts @ ex + shift[0]).max(), (pts @ ey + shift[1]).min(), \
        (pts @ ey + shift[1]).max()
