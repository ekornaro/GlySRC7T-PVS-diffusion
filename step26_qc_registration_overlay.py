#!/usr/bin/env python3
"""
Visual registration check: the DKI MD map warped to T1w space
(dkifit_MD_toT1w_Warped.nii.gz, written by step15) overlaid on the T1w
(<sub>_ses-01_T1wbrain_restore.nii.gz) in 12 axial slices per participant.
One PNG per slice in <pvs_root>/QC_MD_overlays_png/<sub>/. The
co-registration supplementary figure was assembled by hand from these PNGs.

Usage:  python step26_qc_registration_overlay.py
"""
import argparse
from pathlib import Path
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from pipeline_config import cfg


# -----------------------------
# USER SETTINGS
# ----------------------------
BASE_DIR = Path(cfg("pvs_root"))
OUT_DIR = BASE_DIR / "QC_MD_overlays_png"
N_SLICES = 12                 # number of axial slices to export per subject
ALPHA = 0.45                  # overlay transparency
T1_CLIP = (1, 99)             # percentiles for T1 display
MD_CLIP = (5, 99)             # percentiles for MD display
FIG_DPI = 150
MAX_DEPTH = 8                 # how deep to search within each subject folder
# -----------------------------


def rglob_limited(root: Path, pattern: str, max_depth: int):
    """Depth-limited recursive glob."""
    root = root.resolve()
    for p in root.rglob(pattern):
        try:
            rel = p.resolve().relative_to(root)
        except Exception:
            continue
        if len(rel.parts) <= max_depth:
            yield p


def load_canonical_3d(path: Path) -> np.ndarray:
    img = nib.load(str(path))
    img = nib.as_closest_canonical(img)
    data = img.get_fdata(dtype=np.float32)
    if data.ndim != 3:
        raise ValueError(f"Expected 3D NIfTI for {path.name}, got shape {data.shape}")
    return data


def choose_axial_slices_from_brain(t1: np.ndarray, n_slices: int):
    """Choose axial slice indices where there is brain signal."""
    # crude brain mask: non-zero + robust threshold
    nz = t1[t1 > 0]
    if nz.size == 0:
        # fallback: mid slices
        z = t1.shape[2]
        return np.linspace(int(z*0.25), int(z*0.75), n_slices).astype(int)

    thr = np.percentile(nz, 20)  # loose threshold
    mask = t1 > thr
    counts = mask.sum(axis=(0, 1))  # per-slice brain area

    valid = np.where(counts > np.percentile(counts[counts > 0], 20))[0]
    if valid.size < n_slices:
        z = t1.shape[2]
        return np.linspace(int(z*0.25), int(z*0.75), n_slices).astype(int)

    # spread evenly across valid range
    return np.linspace(valid.min(), valid.max(), n_slices).astype(int)


def robust_limits(vol: np.ndarray, pct_lo: float, pct_hi: float):
    v = vol[np.isfinite(vol)]
    v = v[v != 0]
    if v.size == 0:
        return 0.0, 1.0
    lo = np.percentile(v, pct_lo)
    hi = np.percentile(v, pct_hi)
    if hi <= lo:
        hi = lo + 1e-6
    return float(lo), float(hi)


def save_overlay_png(t1: np.ndarray, md: np.ndarray, z: int, out_png: Path,
                     alpha: float, t1_lim, md_lim):
    # Use transpose so "imshow" displays in radiological-ish view consistently
    t1_slice = np.rot90(t1[:, :, z])
    md_slice = np.rot90(md[:, :, z])

    fig = plt.figure(figsize=(6, 6), dpi=FIG_DPI)
    ax = plt.gca()
    ax.imshow(t1_slice, cmap="gray", vmin=t1_lim[0], vmax=t1_lim[1], interpolation="nearest")
    ax.imshow(md_slice, cmap="jet", vmin=md_lim[0], vmax=md_lim[1],
              alpha=alpha, interpolation="nearest")
    ax.set_axis_off()
    fig.tight_layout(pad=0)
    fig.savefig(out_png, bbox_inches="tight", pad_inches=0)
    plt.close(fig)


def main():
    if not BASE_DIR.exists():
        raise FileNotFoundError(f"BASE_DIR does not exist: {BASE_DIR}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    subjects = sorted([p for p in BASE_DIR.iterdir() if p.is_dir() and p.name.startswith("sub-")])
    print(f"Found {len(subjects)} subject folders under: {BASE_DIR}")

    for sub_dir in subjects:
        sub_id = sub_dir.name

        # Find files (depth-limited)
        md_candidates = list(rglob_limited(sub_dir, "dkifit_MD_toT1w_Warped.nii.gz", MAX_DEPTH))
        t1_candidates = list(rglob_limited(sub_dir, f"{sub_id}_ses-01_T1wbrain_restore.nii.gz", MAX_DEPTH))

        if len(md_candidates) == 0 or len(t1_candidates) == 0:
            print(f"[SKIP] {sub_id}: missing files "
                  f"(MD found={len(md_candidates)}, T1 found={len(t1_candidates)})")
            continue

        md_path = md_candidates[0]
        t1_path = t1_candidates[0]

        try:
            md = load_canonical_3d(md_path)
            t1 = load_canonical_3d(t1_path)
        except Exception as e:
            print(f"[FAIL] {sub_id}: could not load/convert to 3D canonical. {e}")
            continue

        if md.shape != t1.shape:
            print(f"[WARN] {sub_id}: shape mismatch T1 {t1.shape} vs MD {md.shape}. "
                  f"This suggests resampling/registration output mismatch. Skipping.")
            continue

        # Robust display limits
        t1_vmin, t1_vmax = robust_limits(t1, *T1_CLIP)
        md_vmin, md_vmax = robust_limits(md, *MD_CLIP)

        z_slices = choose_axial_slices_from_brain(t1, N_SLICES)

        sub_out = OUT_DIR / sub_id
        sub_out.mkdir(parents=True, exist_ok=True)

        for z in z_slices:
            out_png = sub_out / f"{sub_id}_MDonT1w_z{int(z):03d}.png"
            save_overlay_png(
                t1=t1, md=md, z=int(z),
                out_png=out_png,
                alpha=ALPHA,
                t1_lim=(t1_vmin, t1_vmax),
                md_lim=(md_vmin, md_vmax),
            )

        print(f"[OK] {sub_id}: saved {len(z_slices)} PNGs -> {sub_out}")

    print(f"\nDone. Output folder: {OUT_DIR}")


if __name__ == "__main__":
    argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter).parse_args()
    main()
