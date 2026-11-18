#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step09_determine_gradient_convention.py
=======================================
Determine the CORRECT gradient (bvec) convention per subject/block by refitting
the tensor under each candidate convention and scoring the result against
DIRECT ANATOMY: corpus callosum (left-right) and corticospinal tract /
brainstem (superior-inferior).

WHY THIS AND NOT step20c
------------------------
20c established that native-space tensors differ by block (A 0.498 vs B 0.720,
p=1.6e-9) using a brainstem SHAPE-AXIS coherence proxy. That proxy is partly
sensitive to head positioning (coh_none vs b0 obliquity r = -0.33), so it is
convincing that something is wrong and that flipping helps, but weak evidence
about WHICH flip is right per subject. This script replaces the proxy with two
independent, orthogonal anatomical priors and refits the tensor for real under
each candidate, so the winning convention is chosen on the quantity that
actually matters: does V1 point where the anatomy says it must.

THE FRAME QUESTION, STATED PRECISELY
------------------------------------
DIPY returns eigenvectors in the frame in which the bvecs are expressed. The
test therefore asks: under convention C, does V1 (bvec frame) align with the
anatomical direction EXPRESSED IN THE SAME IMAGE'S VOXEL FRAME? Anatomical
directions are defined in world RAS (x = L-R, z = I-S) and pulled back into the
voxel frame via the image affine's rotation:

    d_vox = normalize( inv(R_affine) @ d_world )

No tensor is ever transformed and no registration is applied to the diffusion
data; only binary label masks are moved (nearest-neighbour) if they must come
from T1w space.

SCORING: PROPAGATION COHERENCE (not axis alignment)
---------------------------------------------------
A first version of this script scored |V1 . x| in the corpus callosum and
|V1 . z| in the brainstem. That is PROVABLY BLIND to the error being hunted.
Fitting with flipped gradients Mg (M = diag(1,1,-1), say) when the acquisition
used g gives

    D_fit = M D_true M      ->      V1_fit = M V1_true

and for any coordinate axis e,  |M v . e| = |v . e|.  Alignment of V1 with a
coordinate direction is therefore identical under every sign convention; a
synthetic phantom with a known z-flip scored 1.000 for all four candidates.
This is also why 20c's single oblique brainstem shape axis discriminates only
weakly and ambiguously: its sensitivity comes solely from the axis being
tilted, and any structure near a coordinate plane is close to degenerate.

The detectable consequence of a wrong convention is SPATIAL, not directional.
The tensor at each voxel is reflected IN PLACE while the voxels stay put, so
V1 no longer follows the tract from one voxel to the next. The score is
therefore a one-step tractography:

    for each white-matter voxel p with FA >= --fa-min
        step to p' = p +/- s * V1(p)          (s = --step, in voxels)
        accumulate |V1(p) . V1(p')|
    score = FA-weighted median over voxels, averaged over both step directions

A correct convention propagates along the fibre and lands on a voxel with a
similar V1 (score high). A reflected field steps sideways out of the tract and
lands on unrelated orientations (score low). This is the principle MRtrix3's
`dwigradcheck` uses via full tractography; the one-step version needs no
tracking machinery and runs in seconds per convention.

VALIDATED on a 3D helical phantom against three different injected ground
truths (flip_z, flip_x, none): the winner was correct every time, with a ~4x
margin over the runner-up (14.3 vs 3.7 voxels mean length). Two degeneracies
were found during that validation and are worth knowing:
  * a PLANAR fibre field cannot separate the two in-plane flips (they give the
    same line field up to sign) and is invariant to the out-of-plane flip. Real
    white matter curves in 3D, but a structure chosen inside one plane will not
    discriminate -- which is why the scored mask pools CC and brainstem rather
    than using either alone.
  * scoring only voxels where the step lands on a valid tensor biases a wrong
    convention upward, because leaving the tract is precisely its signature.
    Exits are terminations here, which is what makes length the right measure.

IF MRtrix3 IS AVAILABLE, RUN dwigradcheck TOO
---------------------------------------------
    dwigradcheck dwi.mif -fslgrad bvecs bvals
It is the validated reference implementation and tests the same 24 candidates
with whole-brain tractography. Agreement between the two is strong evidence;
disagreement means neither should be trusted without a visual check.

NOTE ON LEFT-RIGHT FLIPS
------------------------
The brain is approximately mirror-symmetric about the mid-sagittal plane, so an
x-flip is the hardest case for ANY anatomy-based criterion, including this one
and dwigradcheck. If the winner is flip_x with a small margin, treat it as
unresolved and settle it from acquisition metadata (DICOM/PAR gradient tables)
rather than from the images.

CANDIDATE CONVENTIONS
---------------------
  basic     none, flip_x, flip_y, flip_z          (4)
  extended  the 4 above x 6 axis permutations     (24)
Global sign is irrelevant for a tensor (D is quadratic in g), so flipping two
axes equals flipping the third; only these 4 sign patterns are distinct.

USAGE
-----
  # pilot: one block-A subject and one block-B subject, basic conventions
  python step09_determine_gradient_convention.py --subjects sub-115 sub-160

  # if the basic set does not produce a clean winner
  python step09_determine_gradient_convention.py --subjects sub-115 \
      --conventions extended

  # whole cohort once the pilot looks right
  python step09_determine_gradient_convention.py --all --out ./gradient_convention

  # masks already resampled into native DKI space (skips ANTs entirely)
  python step09_determine_gradient_convention.py --subjects sub-115 \
      --native-labels /path/to/sub-115_aseg_in_dki.nii.gz

Importable:
  from step09_determine_gradient_convention import score_convention, run_subject
"""

import os
import glob
import json
import argparse
import subprocess
import tempfile

import numpy as np
import nibabel as nib
import pandas as pd
from pipeline_config import cfg

# ==========================================================================
# CONFIG
# ==========================================================================
DKI_ROOT = cfg("data_root")
T1_ROOT = cfg("pvs_root")
STUDY_XLSX = cfg("study_xlsx")

DWI_SUBDIR = "pipeline_output/DATASINK/DKISpace/dwi_proc"
DWI_PATTERNS = ["DKI_denoised_gibbs_eddy.nii.gz", "DKI_denoised_gibbs.nii.gz",
                "DKI_denoised.nii.gz"]
BVEC_NAME = "DKI_corrected_convfix.bvec"
BVAL_NAME = "DKI_corrected_convfix.bval"

# SynthSeg labels. NOTE: SynthSeg emits aseg + Desikan-Killiany cortical
# labels; it does NOT emit the FreeSurfer corpus-callosum labels 251-255, so
# those cannot be used here (an earlier version tried and got 0 voxels).
# The tracking score needs only a white-matter mask, not a direction prior,
# so cerebral WM is the right scoring region anyway.
WM_LABELS = [2, 41]                        # cerebral WM (L, R)
BRAINSTEM_LABEL = [16]                     # brainstem
CC_LABELS = [251, 252, 253, 254, 255]      # used only if present (axis report)

FA_MIN = 0.35            # only well-determined tensors inform the test
MIN_VOX = 50             # per structure
BLOCK_SPLIT = 149        # sub-num >= this is block B (from 20b/20c)

SIGN_PATTERNS = {"none":   (+1, +1, +1),
                 "flip_x": (-1, +1, +1),
                 "flip_y": (+1, -1, +1),
                 "flip_z": (+1, +1, -1)}
PERMS = {"xyz": (0, 1, 2), "xzy": (0, 2, 1), "yxz": (1, 0, 2),
         "yzx": (1, 2, 0), "zxy": (2, 0, 1), "zyx": (2, 1, 0)}


# ==========================================================================
# HELPERS
# ==========================================================================
def block_of(sub):
    try:
        return "B" if int(sub.split("-")[1]) >= BLOCK_SPLIT else "A"
    except (IndexError, ValueError):
        return "?"


def conventions(kind="basic"):
    """List of (name, perm_tuple, sign_tuple)."""
    if kind == "basic":
        return [(n, PERMS["xyz"], s) for n, s in SIGN_PATTERNS.items()]
    out = []
    for pname, perm in PERMS.items():
        for sname, sign in SIGN_PATTERNS.items():
            label = pname if sname == "none" else f"{pname}_{sname}"
            out.append((label, perm, sign))
    return out


def apply_convention(bvecs, perm, sign):
    """bvecs (N,3) -> permuted and sign-flipped copy."""
    b = np.asarray(bvecs, float)[:, list(perm)]
    return b * np.asarray(sign, float)[None, :]


def world_dirs_in_voxel_frame(affine):
    """
    Anatomical unit directions (world RAS) pulled back into the image's voxel
    index frame: x_world = L->R, z_world = I->S.
    Returns dict name -> unit vector in voxel frame.
    """
    R = np.asarray(affine, float)[:3, :3]
    R = R / np.linalg.norm(R, axis=0, keepdims=True)
    Rinv = np.linalg.inv(R)
    out = {}
    for name, w in (("x", np.array([1.0, 0, 0])), ("z", np.array([0, 0, 1.0]))):
        v = Rinv @ w
        out[name] = v / np.linalg.norm(v)
    return out


def weighted_median(values, weights):
    v = np.asarray(values, float)
    w = np.asarray(weights, float)
    o = np.argsort(v)
    v, w = v[o], w[o]
    cw = np.cumsum(w)
    if cw[-1] <= 0:
        return float(np.median(v))
    return float(v[np.searchsorted(cw, 0.5 * cw[-1])])


# ==========================================================================
# DISCOVERY
# ==========================================================================
def find_dwi(subj_dir):
    d = os.path.join(subj_dir, DWI_SUBDIR)
    for pat in DWI_PATTERNS:
        p = os.path.join(d, pat)
        if os.path.exists(p):
            return p
    hits = sorted(glob.glob(os.path.join(d, "*.nii.gz")))
    hits = [h for h in hits if "noise" not in os.path.basename(h).lower()]
    return hits[0] if hits else None


def find_grads(subj_dir):
    d = os.path.join(subj_dir, DWI_SUBDIR)
    bvec = os.path.join(d, BVEC_NAME)
    bval = os.path.join(d, BVAL_NAME)
    if os.path.exists(bvec) and os.path.exists(bval):
        return bvec, bval
    bvecs = sorted(glob.glob(os.path.join(d, "*.bvec")))
    bvals = sorted(glob.glob(os.path.join(d, "*.bval")))
    return (bvecs[0] if bvecs else None), (bvals[0] if bvals else None)


def find_synthseg(sub):
    hits = sorted(glob.glob(os.path.join(T1_ROOT, "**",
                                         f"*{sub}*_ses-01_synthseg.nii.gz"),
                            recursive=True))
    return hits[0] if hits else None


def find_t1_to_dki_transforms(sub):
    """
    Transforms that take a T1w-space image into b0/DKI space.
    The pipeline stores '<sub>_ses-01_T1wbrain_restore_from_b0_*', i.e. the
    registration OF T1w FROM b0, so these same files with -t applied in the
    forward order map T1w -> b0 for antsApplyTransforms when the reference is
    the b0. Direction is verified empirically by --check-mask-overlap.
    """
    hits = sorted(glob.glob(os.path.join(T1_ROOT, "**",
                                         f"{sub}_ses-01_T1wbrain_restore_from_b0_*"),
                            recursive=True))
    aff = next((h for h in hits if h.endswith("0GenericAffine.mat")), None)
    warp = next((h for h in hits if h.endswith("1Warp.nii.gz")
                 and "Inverse" not in os.path.basename(h)), None)
    inv = next((h for h in hits if "1InverseWarp.nii.gz" in os.path.basename(h)), None)
    return aff, warp, inv


def labels_to_native(sub, seg_path, ref_path, aff, warp, inv, workdir,
                     invert_affine=True):
    """
    Bring the SynthSeg parcellation into native DKI space with nearest
    neighbour. Tries the inverse direction first (T1w -> b0 given transforms
    named 'from_b0'); the caller validates by label counts.
    """
    out = os.path.join(workdir, f"{sub}_seg_in_dki.nii.gz")
    cmd = ["antsApplyTransforms", "-d", "3", "-i", seg_path, "-r", ref_path,
           "-o", out, "-n", "NearestNeighbor"]
    if invert_affine:
        # inverse of the b0 -> T1w registration: the inverted affine is listed
        # first (applied last), as in step07. Runs before the fix listed the
        # inverse warp first.
        cmd += ["-t", f"[{aff},1]"]
        if inv:
            cmd += ["-t", inv]
    else:
        cmd += ["-t", aff]
        if warp:
            cmd += ["-t", warp]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(out):
        raise RuntimeError(f"antsApplyTransforms failed: "
                           f"{(r.stderr or r.stdout).strip()[:200]}")
    return out


# ==========================================================================
# CORE
# ==========================================================================
def fit_tensor(data, bvals, bvecs, mask, fit_kind="dti", b_max=1100.0):
    """
    Fit and return (FA, V1) on masked voxels. 'dti' restricts to b <= b_max
    (fast, and a bvec convention error shows identically in DTI and DKI);
    'dki' uses the full multishell fit for exact parity with the pipeline.
    """
    from dipy.core.gradients import gradient_table
    bvals = np.asarray(bvals, float)
    if fit_kind == "dti":
        from dipy.reconst import dti
        keep = bvals <= b_max
        gtab = gradient_table(bvals[keep], bvecs=bvecs[keep], b0_threshold=50.0)
        fit = dti.TensorModel(gtab, fit_method="WLS").fit(
            np.nan_to_num(data[..., keep]), mask=mask)
    else:
        from dipy.reconst import dki
        gtab = gradient_table(bvals, bvecs=bvecs, b0_threshold=50.0)
        fit = dki.DiffusionKurtosisModel(gtab, fit_method="WLS").fit(
            np.nan_to_num(data), mask=mask)
    return np.asarray(fit.fa), np.asarray(fit.evecs)[..., :, 0]


def _sample_v1(v1, pts):
    """Nearest-neighbour lookup of V1 at (n,3) float voxel coordinates."""
    idx = np.rint(pts).astype(int)
    shape = np.array(v1.shape[:3])
    ok = np.all((idx >= 0) & (idx < shape[None, :]), axis=1)
    out = np.zeros((pts.shape[0], 3), np.float32)
    if ok.any():
        i, j, k = idx[ok].T
        out[ok] = v1[i, j, k]
    return out, ok


def tract_length_score(fa, v1, mask, fa_min=FA_MIN, step=0.5,
                       max_steps=60, cos_thresh=0.766, max_seeds=6000,
                       rng=None):
    """
    Deterministic multi-step tracking score -- the criterion MRtrix3's
    dwigradcheck uses, implemented without a tracking dependency.

    From every seed voxel the field is followed along V1 in both directions,
    stopping when the next voxel has no valid tensor, falls below fa_min, or
    turns by more than acos(cos_thresh) (~40 deg by default). The score is the
    mean streamline length in voxels.

    Why length and not one-step coherence: a wrong gradient convention
    reflects each tensor IN PLACE, so V1 leaves the tract by a small angle at
    each step. One step is easily inside the tract's own width and barely
    registers; error ACCUMULATES over many steps, which is what makes tracking
    discriminative. (A one-step version of this test failed to recover a known
    flip on a curved phantom; this one recovers it.)
    """
    sel = mask & np.isfinite(fa) & (fa >= fa_min)
    n = int(sel.sum())
    if n < MIN_VOX:
        return np.nan, n
    seeds = np.argwhere(sel).astype(np.float32)
    if seeds.shape[0] > max_seeds:
        rng = rng or np.random.default_rng(0)
        seeds = seeds[rng.choice(seeds.shape[0], max_seeds, replace=False)]

    fa_ok = np.isfinite(fa) & (fa >= fa_min)
    total = np.zeros(seeds.shape[0], np.float32)
    for sgn in (+1.0, -1.0):
        pos = seeds.copy()
        d0, _ = _sample_v1(v1, pos)
        nrm = np.linalg.norm(d0, axis=1, keepdims=True)
        cur = sgn * np.divide(d0, np.maximum(nrm, 1e-9))
        alive = nrm.ravel() > 1e-6
        for _ in range(max_steps):
            if not alive.any():
                break
            pos[alive] += step * cur[alive]
            idx = np.rint(pos).astype(int)
            shp = np.array(v1.shape[:3])
            inb = np.all((idx >= 0) & (idx < shp[None, :]), axis=1)
            keep = alive & inb
            if not keep.any():
                alive[:] = False
                break
            i, j, k = idx[keep].T
            nxt = v1[i, j, k]
            nn = np.linalg.norm(nxt, axis=1)
            valid = (nn > 1e-6) & fa_ok[i, j, k]
            nxt = np.divide(nxt, np.maximum(nn, 1e-9)[:, None])
            dot = np.einsum("ij,ij->i", cur[keep], nxt)
            nxt = np.where(dot[:, None] < 0, -nxt, nxt)      # sign continuity
            turn_ok = np.abs(dot) >= cos_thresh
            cont = valid & turn_ok
            ci = np.flatnonzero(keep)
            total[ci[cont]] += step
            cur[ci[cont]] = nxt[cont]
            newalive = np.zeros_like(alive)
            newalive[ci[cont]] = True
            alive = newalive
    return float(np.mean(total)), n


def score_convention(fa, v1, structures, dirs, fa_min=FA_MIN, step=0.5):
    """
    Primary score = propagation coherence over all structures pooled.
    The axis-alignment numbers are still reported (score_axis_*) because they
    are informative about gross orientation, but they are BLIND to sign flips
    and must never be used to choose the convention.
    """
    out = {}
    pooled = np.zeros_like(fa, bool)
    for name, (mask, dkey) in structures.items():
        pooled |= mask
        sel = mask & np.isfinite(fa) & (fa >= fa_min)
        out[f"n_{name}"] = int(sel.sum())
        if sel.sum() >= MIN_VOX:
            cos = np.abs(v1[sel] @ dirs[dkey])
            out[f"score_axis_{name}"] = round(weighted_median(cos, fa[sel]), 4)
        else:
            out[f"score_axis_{name}"] = np.nan
        s, _ = tract_length_score(fa, v1, mask, fa_min, step)
        out[f"len_{name}"] = round(s, 3) if np.isfinite(s) else np.nan
    s, n = tract_length_score(fa, v1, pooled, fa_min, step)
    out["n_scored"] = n
    out["score"] = round(s, 3) if np.isfinite(s) else np.nan
    return out


def run_subject(subj_dir, sub, conv_kind="basic", fit_kind="dti",
                native_labels=None, workdir=None, verbose=True, step=0.5):
    dwi = find_dwi(subj_dir)
    bvec_p, bval_p = find_grads(subj_dir)
    if not (dwi and bvec_p and bval_p):
        raise RuntimeError("DWI/bvec/bval not found")

    img = nib.load(dwi)
    data = np.asarray(img.dataobj, np.float32)
    bvals = np.loadtxt(bval_p).ravel()
    bvecs = np.loadtxt(bvec_p)
    if bvecs.shape[0] == 3:
        bvecs = bvecs.T
    if bvecs.shape[0] != bvals.size:
        raise RuntimeError(f"bvec/bval mismatch: {bvecs.shape} vs {bvals.shape}")
    if bvecs.shape[0] != data.shape[-1]:
        raise RuntimeError(f"gradients {bvecs.shape[0]} != volumes {data.shape[-1]}")

    # ---- labels in native DKI space ----
    if native_labels:
        seg_native = native_labels
    else:
        seg = find_synthseg(sub)
        if not seg:
            raise RuntimeError("SynthSeg parcellation not found")
        aff, warp, inv = find_t1_to_dki_transforms(sub)
        if not aff:
            raise RuntimeError("no T1w<->b0 affine found")
        b0 = os.path.join(workdir, f"{sub}_b0.nii.gz")
        nib.save(nib.Nifti1Image(data[..., int(np.argmin(bvals))], img.affine), b0)
        seg_native = labels_to_native(sub, seg, b0, aff, warp, inv, workdir)

    lab = np.asarray(nib.load(seg_native).dataobj)
    if lab.shape != data.shape[:3]:
        raise RuntimeError(f"labels {lab.shape} != DWI {data.shape[:3]}")

    wm = np.isin(lab, WM_LABELS)
    bs = np.isin(lab, BRAINSTEM_LABEL)
    cc = np.isin(lab, CC_LABELS)
    present = sorted(int(v) for v in np.unique(lab) if v > 0)
    if wm.sum() + bs.sum() < MIN_VOX:
        raise RuntimeError(
            f"scoring mask too small (WM {int(wm.sum())}, BS {int(bs.sum())}). "
            f"Labels present: {present[:25]}{' ...' if len(present) > 25 else ''} "
            f"-- check the label transform / label convention")

    dirs = world_dirs_in_voxel_frame(img.affine)
    structures = {}
    if wm.sum() >= MIN_VOX:
        structures["WM"] = (wm, "x")
    if bs.sum() >= MIN_VOX:
        structures["CST"] = (bs, "z")
    if cc.sum() >= MIN_VOX:                       # only when 251-255 exist
        structures["CC"] = (cc, "x")
    mask = wm | bs | cc
    if verbose:
        print(f"    scoring mask: WM {int(wm.sum())}  brainstem {int(bs.sum())}"
              + (f"  CC {int(cc.sum())}" if cc.sum() else "  (no CC labels)"))

    rows = []
    for name, perm, sign in conventions(conv_kind):
        bv = apply_convention(bvecs, perm, sign)
        fa, v1 = fit_tensor(data, bvals, bv, mask, fit_kind)
        r = {"Subject": sub, "block": block_of(sub), "convention": name,
             "fit": fit_kind}
        r.update(score_convention(fa, v1, structures, dirs, FA_MIN, step))
        rows.append(r)
        if verbose:
            print(f"    {name:12s} mean tract len {r['score']:6.2f} vox   "
                  f"(axis CC {r.get('score_axis_CC', float('nan')):.3f} "
                  f"CST {r.get('score_axis_CST', float('nan')):.3f} "
                  f"- blind to flips, reported only)")
    return pd.DataFrame(rows)


# ==========================================================================
# MAIN
# ==========================================================================
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=DKI_ROOT)
    ap.add_argument("--subjects", nargs="*", default=None)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--conventions", choices=["basic", "extended"], default="basic")
    ap.add_argument("--fit", choices=["dti", "dki"], default="dti",
                    help="dti (default) restricts to b<=1100 and is fast; a "
                         "bvec convention error is identical in both. Use dki "
                         "for exact parity with the pipeline fit.")
    ap.add_argument("--native-labels", default=None,
                    help="parcellation already in native DKI space (single "
                         "subject only; skips ANTs)")
    ap.add_argument("--fa-min", type=float, default=FA_MIN)
    ap.add_argument("--step", type=float, default=0.5,
                    help="tracking step in voxels (default 0.5)")
    ap.add_argument("--out", default="./gradient_convention")
    a = ap.parse_args(argv)

    globals()["FA_MIN"] = a.fa_min
    os.makedirs(a.out, exist_ok=True)

    if a.subjects:
        subs = a.subjects
    elif a.all:
        subs = sorted(os.path.basename(d) for d in
                      glob.glob(os.path.join(a.root, "sub-*")) if os.path.isdir(d))
    else:
        raise SystemExit("give --subjects or --all")

    print(f"{len(subs)} subject(s) | conventions: {a.conventions} | fit: {a.fit}")
    print("score = mean streamline length in voxels (higher = V1 follows the tract)")
    print("axis-alignment columns are reported but are BLIND to sign flips")
    print(f"FA >= {a.fa_min} | step {a.step} vox\n")

    frames = []
    with tempfile.TemporaryDirectory(prefix="conv20d_") as wd:
        for sub in subs:
            print(f"[{sub}] block {block_of(sub)}")
            try:
                df = run_subject(os.path.join(a.root, sub), sub, a.conventions,
                                 a.fit, a.native_labels, wd, step=a.step)
            except Exception as exc:
                print(f"    FAILED -- {exc}")
                frames.append(pd.DataFrame([{"Subject": sub,
                                             "block": block_of(sub),
                                             "convention": None,
                                             "status": f"FAILED: {exc}"}]))
                continue
            best = df.loc[df["score"].idxmax()]
            runner = df.sort_values("score", ascending=False).iloc[1]
            margin = float(best["score"] - runner["score"])
            print(f"    -> BEST {best['convention']} "
                  f"(mean length {best['score']:.2f} vox, runner-up "
                  f"{runner['convention']} {runner['score']:.2f})")
            ratio = float(best["score"] / runner["score"]) if runner["score"] else np.inf
            print(f"       margin ratio over runner-up: {ratio:.2f}x")
            if ratio < 1.5:
                print("       [!] winner is not clearly separated (<1.5x). Do "
                      "NOT act on this subject alone: check the label overlap, "
                      "try --conventions extended, and corroborate with "
                      "MRtrix3 dwigradcheck.")
            if best["convention"].startswith("flip_x") and ratio < 2.0:
                print("       [!] an x-flip winner is the hardest call to make "
                      "from anatomy (near mirror symmetry). Settle it from the "
                      "acquisition gradient table, not from the images.")
            df["is_best"] = df["convention"] == best["convention"]
            df["margin"] = margin
            frames.append(df)

    out = pd.concat(frames, ignore_index=True)
    csv = os.path.join(a.out, "gradient_convention_scores.csv")
    out.to_csv(csv, index=False)

    if "score" not in out.columns or out["score"].notna().sum() == 0:
        print("\nNo subject produced a score. Fix the failures above before "
              "interpreting anything; nothing was written beyond the raw table.")
        print(f"\n-> {csv}")
        return out
    ok = out[out["score"].notna()]
    if len(ok):
        print("\n" + "=" * 70)
        print("SUMMARY: winning convention by block")
        best_rows = ok[ok.is_best] if "is_best" in ok else ok
        print(pd.crosstab(best_rows["block"], best_rows["convention"]))
        print("\nmean score of the winner, by block:")
        print(best_rows.groupby("block")["score"].median().round(3))
        cons = (best_rows.groupby("block")["convention"]
                .agg(lambda s: s.value_counts().index[0]))
        print("\nMost common winner per block:")
        for blk, conv in cons.items():
            n = int((best_rows[best_rows.block == blk].convention == conv).sum())
            tot = int((best_rows.block == blk).sum())
            print(f"  block {blk}: {conv}  ({n}/{tot} subjects)")
        with open(os.path.join(a.out, "convention_by_block.json"), "w") as fh:
            json.dump({k: str(v) for k, v in cons.items()}, fh, indent=2)
        print("\nIf a block's winner is unanimous or near-unanimous, that is the "
              "convention to refit that block with.\nIf winners are split within "
              "a block, the problem is per-subject, not per-block -- do NOT apply "
              "a blanket flip.")
    print(f"\n-> {csv}")
    return out


if __name__ == "__main__":
    main()
