#!/usr/bin/env python3
"""Scale-aware Joukowsky-Topological Pore Descriptor (JTPD) analysis.

The analysis is restricted to segmented 2-D rock-section fields.  It computes
physically calibrated geometric/topological features, a contour-preserving
Joukowsky shape response, synthetic topology checks, and resolution sensitivity.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import contourpy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage

from rock_joukowsky import anisotropy, joukowsky_metrics_v2, polygon_area, polygon_perimeter

ROOT = Path(__file__).resolve().parent
DATA_ROOT = ROOT / "data"
MASK_ROOT = DATA_ROOT / "structures"
METADATA = DATA_ROOT / "Lam_065_metadata.json"
OUT = ROOT / "results"
FIG = ROOT / "figures"
OUT.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

MASKS = [
    MASK_ROOT / "patch_y3800_x3800_c0_mask.png",
    MASK_ROOT / "patch_y7600_x19000_c0_mask.png",
    MASK_ROOT / "patch_y7600_x53200_c0_mask.png",
]
NAMES = ["R1", "R2", "R3"]
A_GRID = np.array([0.00, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35], dtype=float)
MAIN_TARGET = 256
SCALE_TARGETS = [128, 256, 512, 1024]
MIN_SHAPE_AREA_UM2 = 800.0
N_CONTOUR = 128

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10.5,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "legend.fontsize": 9,
    "figure.dpi": 150,
})

def metadata():
    return json.loads(METADATA.read_text(encoding="utf-8"))

def source_pixel_um():
    m = metadata()
    return float(m["physical_scaling"]["pixel_size_x_meters"]) * 1e6

def load_solid_mask(path: Path, target: int) -> np.ndarray:
    rgb = np.asarray(Image.open(path).convert("RGB"))
    solid = (rgb[:, :, 1].astype(int) > rgb[:, :, 0].astype(int) + 20) & (
        rgb[:, :, 1].astype(int) > rgb[:, :, 2].astype(int) + 20
    )
    h, w = solid.shape
    if target == h and target == w:
        return solid
    if h % target == 0 and w % target == 0:
        fy, fx = h // target, w // target
        reduced = solid.reshape(target, fy, target, fx).mean(axis=(1, 3))
        return reduced >= 0.5
    ys = np.linspace(0, h, target + 1, dtype=int)
    xs = np.linspace(0, w, target + 1, dtype=int)
    out = np.zeros((target, target), bool)
    for iy in range(target):
        for ix in range(target):
            out[iy, ix] = solid[ys[iy]:ys[iy+1], xs[ix]:xs[ix+1]].mean() >= 0.5
    return out

def effective_spacing_um(path: Path, target: int) -> float:
    with Image.open(path) as im:
        width = im.size[0]
    return source_pixel_um() * width / target

def label(mask: np.ndarray, connectivity: int):
    st = ndimage.generate_binary_structure(2, 2 if connectivity == 8 else 1)
    return ndimage.label(mask, structure=st)

def topology(mask: np.ndarray) -> dict:
    fg, beta0 = label(mask, 8)
    bg, n_bg = label(~mask, 4)
    border_ids = set(np.unique(np.r_[bg[0, :], bg[-1, :], bg[:, 0], bg[:, -1]]))
    border_ids.discard(0)
    beta1 = sum(1 for k in range(1, n_bg + 1) if k not in border_ids)
    px = py = 0
    for k in range(1, beta0 + 1):
        c = fg == k
        px += int(c[:, 0].any() and c[:, -1].any())
        py += int(c[0, :].any() and c[-1, :].any())
    return {
        "beta0": int(beta0),
        "beta1": int(beta1),
        "euler": int(beta0 - beta1),
        "percolates_x": int(px),
        "percolates_y": int(py),
    }

def phase_stats(mask: np.ndarray) -> dict:
    lab, n = label(mask, 8)
    sizes = np.bincount(lab.ravel())[1:] if n else np.array([], dtype=int)
    return {
        "fraction": float(mask.mean()),
        "n_components": int(n),
        "largest_component_fraction": float(sizes.max() / mask.sum()) if sizes.size else 0.0,
    }

def interface_length(mask: np.ndarray, spacing_um: float) -> float:
    transitions = sum(np.abs(np.diff(mask.astype(np.int8), axis=ax)).sum() for ax in (0, 1))
    return float(transitions * spacing_um)

def _clean_closed_line(points: np.ndarray) -> np.ndarray:
    p = np.asarray(points, float)
    if len(p) >= 2 and np.linalg.norm(p[0] - p[-1]) < 1e-9:
        p = p[:-1]
    keep = np.r_[True, np.linalg.norm(np.diff(p, axis=0), axis=1) > 1e-9]
    return p[keep]

def resample_closed(points: np.ndarray, n: int = N_CONTOUR) -> np.ndarray:
    p = _clean_closed_line(points)
    if len(p) < 3:
        raise ValueError("contour has fewer than 3 distinct points")
    q = np.vstack([p, p[0]])
    seg = np.linalg.norm(np.diff(q, axis=0), axis=1)
    if seg.sum() <= 0:
        raise ValueError("degenerate contour")
    cum = np.r_[0.0, np.cumsum(seg)]
    t = np.linspace(0, cum[-1], n, endpoint=False)
    out = np.empty((n, 2), float)
    j = 0
    for i, ti in enumerate(t):
        while j + 1 < len(cum) and cum[j + 1] <= ti:
            j += 1
        f = (ti - cum[j]) / max(seg[j], 1e-15)
        out[i] = q[j] + f * (q[j + 1] - q[j])
    return out

def outer_contour(component: np.ndarray, spacing_um: float) -> np.ndarray:
    padded = np.pad(component.astype(float), 1)
    cg = contourpy.contour_generator(z=padded)
    lines = cg.lines(0.5)
    candidates = []
    for line in lines:
        p = _clean_closed_line(np.asarray(line, float) - 1.0)
        if len(p) >= 3:
            candidates.append(p)
    if not candidates:
        raise ValueError("no contour")
    p = max(candidates, key=lambda x: polygon_area(x))
    return resample_closed(p * spacing_um)

def shape_features(pore: np.ndarray, spacing_um: float) -> list[dict]:
    lab, n = label(pore, 8)
    sizes = np.bincount(lab.ravel())[1:] if n else np.array([], dtype=int)
    min_pix = max(9, int(math.ceil(MIN_SHAPE_AREA_UM2 / (spacing_um ** 2))))
    rows = []
    for idx in np.argsort(sizes)[::-1]:
        npix = int(sizes[idx])
        if npix < min_pix:
            continue
        comp = lab == int(idx + 1)
        try:
            poly = outer_contour(comp, spacing_um)
            a0 = polygon_area(poly)
            l0 = polygon_perimeter(poly)
            if a0 <= 0 or l0 <= 0:
                continue
            row = {
                "component_pixels": npix,
                "component_area_um2": npix * spacing_um ** 2,
                "contour_area_um2": a0,
                "contour_perimeter_um": l0,
                "anisotropy": anisotropy(poly),
            }
            for a in A_GRID:
                jm = joukowsky_metrics_v2(poly, a=float(a))
                row[f"ja_{a:.2f}"] = jm["j_area"] / a0
                row[f"jl_{a:.2f}"] = jm["j_perimeter"] / l0
            rows.append(row)
        except (ValueError, FloatingPointError, np.linalg.LinAlgError):
            continue
    return rows

def sample_record(name: str, path: Path, target: int, include_shape: bool = True) -> dict:
    solid = load_solid_mask(path, target)
    pore = ~solid
    spacing = effective_spacing_um(path, target)
    side_mm = target * spacing / 1000.0
    area_mm2 = side_mm ** 2
    ts = topology(pore)
    ps = phase_stats(pore)
    ss = phase_stats(solid)
    rec = {
        "sample": name,
        "target_px": target,
        "effective_pixel_um": spacing,
        "field_side_mm": side_mm,
        "field_area_mm2": area_mm2,
        "pore_area_fraction": ps["fraction"],
        "solid_area_fraction": ss["fraction"],
        "pore_components": ps["n_components"],
        "largest_pore_fraction": ps["largest_component_fraction"],
        "solid_components": ss["n_components"],
        "largest_solid_fraction": ss["largest_component_fraction"],
        **ts,
        "beta0_density_mm2": ts["beta0"] / area_mm2,
        "beta1_density_mm2": ts["beta1"] / area_mm2,
        "euler_density_mm2": ts["euler"] / area_mm2,
        "interface_length_mm_per_mm2": (interface_length(solid, spacing) / 1000.0) / area_mm2,
    }
    if include_shape:
        shapes = shape_features(pore, spacing)
        rec["shape_components"] = shapes
        for a in A_GRID:
            for prefix in ("ja", "jl"):
                vals = np.array([r[f"{prefix}_{a:.2f}"] for r in shapes], float)
                rec[f"{prefix}_{a:.2f}_median"] = float(np.median(vals)) if vals.size else np.nan
                rec[f"{prefix}_{a:.2f}_q25"] = float(np.quantile(vals, .25)) if vals.size else np.nan
                rec[f"{prefix}_{a:.2f}_q75"] = float(np.quantile(vals, .75)) if vals.size else np.nan
    return rec

def validate_topology() -> dict:
    n = 96
    y, x = np.ogrid[:n, :n]
    disk = (x-48)**2 + (y-48)**2 <= 24**2
    annulus = ((x-48)**2 + (y-48)**2 <= 28**2) & ((x-48)**2 + (y-48)**2 >= 12**2)
    two = ((x-30)**2 + (y-48)**2 <= 13**2) | ((x-66)**2 + (y-48)**2 <= 13**2)
    channel = np.zeros((n, n), bool); channel[40:56, :] = True
    cases = {
        "disk": (disk, (1, 0, 0, 0)),
        "annulus": (annulus, (1, 1, 0, 0)),
        "two_disks": (two, (2, 0, 0, 0)),
        "x_channel": (channel, (1, 0, 1, 0)),
    }
    out = {}
    for key, (mask, expected) in cases.items():
        t = topology(mask)
        got = (t["beta0"], t["beta1"], int(t["percolates_x"] > 0), int(t["percolates_y"] > 0))
        out[key] = {"expected": expected, "got": got, "pass": got == expected}
    return out

def validate_joukowsky_invariance() -> dict:
    p = np.array([[0,0],[2.1,.2],[2.7,1.4],[1.3,2.4],[-.4,1.7],[-.8,.5]], float)
    def ratios(q):
        basea, basel = polygon_area(q), polygon_perimeter(q)
        j = joukowsky_metrics_v2(q, .25)
        return np.array([j["j_area"]/basea, j["j_perimeter"]/basel, j["j_anisotropy"]])
    ref = ratios(p)
    th = .73
    R = np.array([[np.cos(th), -np.sin(th)],[np.sin(th), np.cos(th)]])
    variants = {
        "translation": p + [10.0, -7.0],
        "scale": p * 4.3,
        "rotation": p @ R.T,
        "reflection": p * [-1, 1],
        "reversal": p[::-1],
    }
    errors = {}
    for k, q in variants.items():
        val = ratios(q)
        errors[k] = float(np.max(np.abs(val-ref) / np.maximum(np.abs(ref), 1e-15)))
    return {"reference": ref.tolist(), "max_relative_error": errors, "overall_max": max(errors.values())}

def figure_workflow(main):
    r = main[1]
    solid = load_solid_mask(MASKS[1], MAIN_TARGET)
    pore = ~solid
    lab, _ = label(pore, 8)
    fig, ax = plt.subplots(1, 4, figsize=(14.5, 3.6), constrained_layout=True)
    ax[0].imshow(solid, cmap="gray_r", interpolation="nearest")
    ax[0].set_title("(a) Binary rock section")
    ax[0].axis("off")
    ax[1].imshow(lab, cmap="nipy_spectral", interpolation="nearest")
    ax[1].set_title("(b) Pore components")
    ax[1].axis("off")
    shapes = r["shape_components"]
    if shapes:
        # Use a synthetic irregular contour only for a clean illustration of the mapping.
        th = np.linspace(0, 2*np.pi, 128, endpoint=False)
        rr = 1 + .22*np.cos(3*th) + .12*np.sin(2*th)
        poly = np.c_[rr*np.cos(th), rr*np.sin(th)]
        jm = []
        centered = poly - poly.mean(axis=0)
        z = centered[:,0] + 1j*centered[:,1]
        rad = np.mean(np.abs(z)); u = z/rad
        edge = z[1]-z[0]; oriented = u*np.conj(edge/abs(edge))
        w = rad*(oriented + .25**2/oriented)
        tr = np.c_[w.real,w.imag]
        ax[2].plot(poly[:,0], poly[:,1], lw=2, label="original")
        ax[2].plot(tr[:,0], tr[:,1], lw=2, label=r"$J_{0.25}$")
        ax[2].axis("equal"); ax[2].legend(frameon=False)
    ax[2].set_title("(c) Nonlinear contour probe")
    ax[2].set_xticks([]); ax[2].set_yticks([])
    ax[3].axis("off")
    ax[3].text(.02,.92, r"$\mathbf{D}_{JT}(\ell)=$", fontsize=14, va="top")
    ax[3].text(.02,.76, "pore-area fraction\ncomponent density\nBetti densities / Euler density\ndirectional spanning\ninterface density\nanisotropy\nJoukowsky area & perimeter response",
               va="top", linespacing=1.45)
    ax[3].text(.02,.08, f"Observation scale: {r['effective_pixel_um']:.2f} µm px⁻¹", fontsize=9)
    ax[3].set_title("(d) Interpretable feature vector")
    for ext in ("pdf","png"):
        fig.savefig(FIG/f"figure1_workflow.{ext}", dpi=500 if ext=="png" else None, bbox_inches="tight")
    plt.close(fig)

def figure_fields(main):
    fig = plt.figure(figsize=(13.5, 6.8), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.15, .85])
    for j,(name,path,r) in enumerate(zip(NAMES,MASKS,main)):
        solid = load_solid_mask(path, MAIN_TARGET)
        a = fig.add_subplot(gs[0,j])
        a.imshow(solid, cmap="gray_r", interpolation="nearest")
        a.set_title(f"{name}: $\\phi_p$={r['pore_area_fraction']:.3f}")
        a.axis("off")
        # 500 µm scale bar
        bar_px = 500.0/r["effective_pixel_um"]
        y = MAIN_TARGET-13; x0=14
        a.plot([x0,x0+bar_px],[y,y],color="white",lw=4,solid_capstyle="butt")
        a.text(x0+bar_px/2,y-7,"500 µm",color="white",ha="center",va="bottom",fontsize=8)
    x = np.arange(3)
    a = fig.add_subplot(gs[1,0])
    a.bar(x,[r["pore_area_fraction"] for r in main])
    a.set_xticks(x,NAMES); a.set_ylim(0,1)
    a.set_ylabel("2-D pore-area fraction")
    a.set_title("(d) Pore abundance")
    a = fig.add_subplot(gs[1,1])
    width=.36
    a.bar(x-width/2,[r["beta0_density_mm2"] for r in main],width,label=r"$\beta_0/A$")
    a.bar(x+width/2,[r["beta1_density_mm2"] for r in main],width,label=r"$\beta_1/A$")
    a.set_xticks(x,NAMES); a.set_ylabel("mm$^{-2}$"); a.legend(frameon=False)
    a.set_title("(e) Topological densities")
    a = fig.add_subplot(gs[1,2])
    M=np.array([[int(r["percolates_x"]>0),int(r["percolates_y"]>0)] for r in main])
    im=a.imshow(M,vmin=0,vmax=1,cmap="Greys",aspect="auto")
    a.set_xticks([0,1],["X spanning","Y spanning"]); a.set_yticks(x,NAMES)
    for i in range(3):
        for j in range(2):
            a.text(j,i,"yes" if M[i,j] else "no",ha="center",va="center",
                   color="white" if M[i,j] else "black",fontweight="bold")
    a.set_title("(f) Directional connectivity")
    for ext in ("pdf","png"):
        fig.savefig(FIG/f"figure2_fields_topology.{ext}", dpi=500 if ext=="png" else None, bbox_inches="tight")
    plt.close(fig)

def figure_joukowsky(main):
    fig, ax = plt.subplots(1,3,figsize=(13.8,4.1),constrained_layout=True)
    for r in main:
        med=[r[f"ja_{a:.2f}_median"] for a in A_GRID]
        lo=[r[f"ja_{a:.2f}_q25"] for a in A_GRID]
        hi=[r[f"ja_{a:.2f}_q75"] for a in A_GRID]
        ax[0].plot(A_GRID,med,marker="o",label=r["sample"])
        ax[0].fill_between(A_GRID,lo,hi,alpha=.15)
        medl=[r[f"jl_{a:.2f}_median"] for a in A_GRID]
        lol=[r[f"jl_{a:.2f}_q25"] for a in A_GRID]
        hil=[r[f"jl_{a:.2f}_q75"] for a in A_GRID]
        ax[1].plot(A_GRID,medl,marker="o",label=r["sample"])
        ax[1].fill_between(A_GRID,lol,hil,alpha=.15)
    ax[0].axhline(1,lw=.8,ls="--"); ax[1].axhline(1,lw=.8,ls="--")
    ax[0].set(xlabel=r"Joukowsky strength $\alpha$",ylabel=r"$A_J/A$",title="(a) Area response")
    ax[1].set(xlabel=r"Joukowsky strength $\alpha$",ylabel=r"$L_J/L$",title="(b) Perimeter response")
    ax[0].legend(frameon=False); ax[1].legend(frameon=False)
    vals=[]
    for r in main:
        vals.append([q["ja_0.25"] for q in r["shape_components"]])
    ax[2].boxplot(vals,tick_labels=NAMES,showfliers=False)
    ax[2].axhline(1,lw=.8,ls="--")
    ax[2].set(ylabel=r"$A_J/A$ at $\alpha=0.25$",title="(c) Component-level response")
    for a in ax: a.grid(alpha=.2)
    for ext in ("pdf","png"):
        fig.savefig(FIG/f"figure3_joukowsky_response.{ext}", dpi=500 if ext=="png" else None, bbox_inches="tight")
    plt.close(fig)

def figure_scale(scale):
    df=pd.DataFrame(scale)
    fig,ax=plt.subplots(2,2,figsize=(10.5,7.4),constrained_layout=True)
    for name in NAMES:
        d=df[df["sample"]==name].sort_values("effective_pixel_um")
        ax[0,0].plot(d["effective_pixel_um"],d["pore_area_fraction"],marker="o",label=name)
        ax[0,1].plot(d["effective_pixel_um"],d["beta0_density_mm2"],marker="o",label=name)
        ax[1,0].plot(d["effective_pixel_um"],d["beta1_density_mm2"],marker="o",label=name)
        ax[1,1].plot(d["effective_pixel_um"],d["interface_length_mm_per_mm2"],marker="o",label=name)
    ax[0,0].set_ylabel("2-D pore-area fraction"); ax[0,0].set_title("(a) Pore abundance")
    ax[0,1].set_ylabel(r"$\beta_0/A$ [mm$^{-2}$]"); ax[0,1].set_title("(b) Component density")
    ax[1,0].set_ylabel(r"$\beta_1/A$ [mm$^{-2}$]"); ax[1,0].set_title("(c) Loop density")
    ax[1,1].set_ylabel(r"Interface [mm mm$^{-2}$]"); ax[1,1].set_title("(d) Interface density")
    for a in ax.flat:
        a.set_xlabel("Effective pixel size [µm]")
        a.invert_xaxis()
        a.grid(alpha=.22)
    ax[0,0].legend(frameon=False,ncol=3)
    for ext in ("pdf","png"):
        fig.savefig(FIG/f"figure4_scale_sensitivity.{ext}", dpi=500 if ext=="png" else None, bbox_inches="tight")
    plt.close(fig)

def main():
    main_records=[sample_record(n,p,MAIN_TARGET,True) for n,p in zip(NAMES,MASKS)]
    scale_records=[sample_record(n,p,t,False) for n,p in zip(NAMES,MASKS) for t in SCALE_TARGETS]
    topo_check=validate_topology()
    inv=validate_joukowsky_invariance()
    if not all(v["pass"] for v in topo_check.values()):
        raise RuntimeError(f"Topology validation failed: {topo_check}")
    # Compact tables.
    compact_keys=[
        "sample","effective_pixel_um","field_side_mm","pore_area_fraction",
        "pore_components","largest_pore_fraction","beta0","beta1","euler",
        "beta0_density_mm2","beta1_density_mm2","euler_density_mm2",
        "percolates_x","percolates_y","interface_length_mm_per_mm2",
        "ja_0.25_median","jl_0.25_median",
    ]
    pd.DataFrame([{k:r.get(k) for k in compact_keys} for r in main_records]).to_csv(OUT/"jtpd_summary.csv",index=False)
    pd.DataFrame(scale_records).drop(columns=["shape_components"],errors="ignore").to_csv(OUT/"scale_sensitivity.csv",index=False)
    payload={
        "metadata":metadata(),
        "method":{
            "main_target_px":MAIN_TARGET,
            "source_pixel_um":source_pixel_um(),
            "effective_main_pixel_um":effective_spacing_um(MASKS[0],MAIN_TARGET),
            "connectivity":"8-connected pore foreground; 4-connected complement for holes",
            "min_shape_area_um2":MIN_SHAPE_AREA_UM2,
            "contour_samples":N_CONTOUR,
            "joukowsky_alpha_grid":A_GRID.tolist(),
        },
        "samples":main_records,
        "topology_validation":topo_check,
        "joukowsky_invariance":inv,
    }
    (OUT/"jtpd_results.json").write_text(json.dumps(payload,indent=2,ensure_ascii=False),encoding="utf-8")
    figure_workflow(main_records)
    figure_fields(main_records)
    figure_joukowsky(main_records)
    figure_scale(scale_records)
    print(pd.read_csv(OUT/"jtpd_summary.csv").to_string(index=False))
    print("\\nTopology checks:",topo_check)
    print("Joukowsky invariance max relative error:",inv["overall_max"])

if __name__=="__main__":
    main()
