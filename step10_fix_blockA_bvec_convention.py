#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step10_fix_blockA_bvec_convention.py
====================================
Write a z-flipped copy of the stored bvecs for the block-A subjects, based on
the step20d determination:

    block A (sub-102 .. sub-146):  stored bvecs need the z row NEGATED
    block B (sub-149 .. sub-172):  stored bvecs are correct as-is

EVIDENCE (step20d, 65/65 subjects, zero exceptions)
---------------------------------------------------
Tract-length scoring of a refit under all four sign conventions:
    block A: 41/41 subjects win with flip_z  (margin >= 1.23x, median 1.5x)
    block B: 24/24 subjects win with none    (margin >= 1.35x, median 1.5x)
    winner tract length A 9.94 vs B 9.81 voxels -> the flip restores block A
    to block B's quality exactly.
Independently concordant with step20c's block-level verdict (coherence 0.498
-> 0.701 under flip_z in block A). NOTE this gradient-convention boundary
(146|149) is NOT the b0-vs-T1w reflection boundary from step20b (140|141);
the two block-associated quirks are separate phenomena.

WHAT THIS SCRIPT DOES (and deliberately does not do)
----------------------------------------------------
For every block-A subject it finds <dwi_proc>/DKI_corrected.bvec and writes

    DKI_corrected_convfix.bvec        (z row negated)
    DKI_corrected_convfix.bval        (unchanged copy, so the pair travels
                                       together and downstream globbing can
                                       key on the _convfix stem)
    bvec_convention_provenance.json   (what was done, from what, when, why)

For block-B subjects it ALSO writes _convfix files -- as UNCHANGED copies --
so that after this script runs, EVERY subject has a DKI_corrected_convfix.*
pair and the refit (step14) can be pointed at one stem for the whole cohort.
The provenance json records flip_z vs identity per subject.

THE ORIGINALS ARE NEVER MODIFIED. The stored DKI_corrected.bvec is the
provenance of what the scanner/pipeline produced; the _convfix copy is the
analysis input. Anything that edits originals in place destroys the ability
to audit this correction later.

SAFETY
------
* dry-run by default; --apply to write
* refuses to overwrite an existing _convfix file unless --overwrite
* verifies after writing: reloads the pair, checks shapes match the original,
  checks z row == -original (A) or == original (B), checks unit norms
  preserved, and checks the bval copy is byte-identical
* block membership is an explicit table in this file (SUBJECT ASSIGNMENT
  section), not an inequality on the subject number, so sub-147/148 style
  gaps and future subjects cannot be silently misclassified

USAGE
-----
  python step10_fix_blockA_bvec_convention.py                # dry run
  python step10_fix_blockA_bvec_convention.py --apply
  python step10_fix_blockA_bvec_convention.py --subjects sub-115 --apply

AFTER THIS: re-run step14 with the tensor fit pointed at the _convfix pair
for ALL subjects (one code path, no per-block branching in the fit), then
re-run step20d on the refit tensors: every subject should now win 'none'.
"""

import os
import glob
import json
import argparse
import datetime
import shutil

import numpy as np
from pipeline_config import cfg

# ==========================================================================
# CONFIG
# ==========================================================================
DKI_ROOT = cfg("data_root")
DWI_SUBDIR = "pipeline_output/DATASINK/DKISpace/dwi_proc"
BVEC_NAME = "DKI_corrected.bvec"
BVAL_NAME = "DKI_corrected.bval"
FIX_STEM = "DKI_corrected_convfix"

# ==========================================================================
# SUBJECT ASSIGNMENT (explicit, from step20d/20c -- not an inequality)
# ==========================================================================
BLOCK_A = [102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114,
           115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127,
           129, 130, 131, 132, 133, 134, 135, 136, 138, 139, 140, 141, 142,
           143, 144, 145, 146]
BLOCK_B = [149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161,
           162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172]
ACTION = {**{n: "flip_z" for n in BLOCK_A}, **{n: "identity" for n in BLOCK_B}}


def subject_action(sub):
    """'flip_z' | 'identity' | None (unknown subject -> refuse, don't guess)."""
    try:
        return ACTION.get(int(sub.split("-")[1]))
    except (IndexError, ValueError):
        return None


# ==========================================================================
# CORE
# ==========================================================================
def load_bvec(path):
    b = np.loadtxt(path)
    if b.ndim != 2 or 3 not in b.shape:
        raise ValueError(f"{path}: unexpected shape {b.shape}")
    transposed = b.shape[1] == 3 and b.shape[0] != 3   # (N,3) storage
    return (b.T if transposed else b), transposed      # always return (3,N)


def plan_subject(subj_dir, sub, overwrite=False):
    """Return an action dict (no writes)."""
    action = subject_action(sub)
    if action is None:
        return {"Subject": sub, "status": "UNKNOWN SUBJECT -- not in the "
                "block table; add it to BLOCK_A/BLOCK_B explicitly"}
    d = os.path.join(subj_dir, DWI_SUBDIR)
    bvec, bval = os.path.join(d, BVEC_NAME), os.path.join(d, BVAL_NAME)
    if not os.path.exists(bvec) or not os.path.exists(bval):
        return {"Subject": sub, "status": f"MISSING {BVEC_NAME}/{BVAL_NAME}"}
    out_bvec = os.path.join(d, FIX_STEM + ".bvec")
    out_bval = os.path.join(d, FIX_STEM + ".bval")
    exists = os.path.exists(out_bvec) or os.path.exists(out_bval)
    if exists and not overwrite:
        return {"Subject": sub, "status": "SKIP: _convfix exists "
                "(use --overwrite to redo)", "action": action}
    return {"Subject": sub, "status": "ok", "action": action,
            "bvec": bvec, "bval": bval,
            "out_bvec": out_bvec, "out_bval": out_bval}


def apply_subject(plan):
    """Write the _convfix pair + provenance, then verify. Returns status str."""
    B, transposed = load_bvec(plan["bvec"])            # (3, N)
    fixed = B.copy()
    if plan["action"] == "flip_z":
        fixed[2, :] *= -1.0

    out = fixed.T if transposed else fixed             # preserve storage layout
    np.savetxt(plan["out_bvec"], out, fmt="%.10g")
    shutil.copyfile(plan["bval"], plan["out_bval"])

    prov = {"script": "step10_fix_blockA_bvec_convention.py",
            "action": plan["action"],
            "source_bvec": os.path.basename(plan["bvec"]),
            "evidence": "step20d gradient_convention_scores.csv "
                        "(41/41 A=flip_z, 24/24 B=none)",
            "written": datetime.datetime.now().isoformat(timespec="seconds")}
    with open(os.path.join(os.path.dirname(plan["out_bvec"]),
                           "bvec_convention_provenance.json"), "w") as fh:
        json.dump(prov, fh, indent=2)

    # ---- verification ----
    C, _ = load_bvec(plan["out_bvec"])
    if C.shape != B.shape:
        return f"VERIFY FAIL: shape {C.shape} != {B.shape}"
    want = B.copy()
    if plan["action"] == "flip_z":
        want[2, :] *= -1.0
    if not np.allclose(C, want, atol=1e-8):
        return "VERIFY FAIL: written bvec does not match intended transform"
    n_old = np.linalg.norm(B, axis=0)
    n_new = np.linalg.norm(C, axis=0)
    if not np.allclose(n_old, n_new, atol=1e-8):
        return "VERIFY FAIL: gradient norms changed"
    with open(plan["bval"], "rb") as f1, open(plan["out_bval"], "rb") as f2:
        if f1.read() != f2.read():
            return "VERIFY FAIL: bval copy not byte-identical"
    return "written+verified"


# ==========================================================================
# MAIN
# ==========================================================================
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=DKI_ROOT)
    ap.add_argument("--subjects", nargs="*", default=None)
    ap.add_argument("--apply", action="store_true",
                    help="actually write; default is dry run")
    ap.add_argument("--overwrite", action="store_true",
                    help="redo existing _convfix files")
    a = ap.parse_args(argv)

    dirs = ([os.path.join(a.root, s) for s in a.subjects] if a.subjects
            else sorted(glob.glob(os.path.join(a.root, "sub-*"))))
    mode = "APPLY" in ("APPLY" if a.apply else "") or a.apply
    print(f"Mode: {'APPLY' if a.apply else 'DRY RUN (nothing will be written)'}")
    print(f"{len(dirs)} subjects | flip_z for block A "
          f"({len(BLOCK_A)} listed), identity copy for block B "
          f"({len(BLOCK_B)} listed)\n")

    n_flip = n_copy = n_skip = n_fail = 0
    for d in dirs:
        sub = os.path.basename(d.rstrip("/"))
        plan = plan_subject(d, sub, a.overwrite)
        if plan["status"] != "ok":
            print(f"  {sub}: {plan['status']}")
            n_skip += plan["status"].startswith("SKIP")
            n_fail += not plan["status"].startswith("SKIP")
            continue
        if not a.apply:
            print(f"  {sub}: would write {FIX_STEM}.bvec/.bval  "
                  f"[{plan['action']}]")
            n_flip += plan["action"] == "flip_z"
            n_copy += plan["action"] == "identity"
            continue
        res = apply_subject(plan)
        print(f"  {sub}: {res}  [{plan['action']}]")
        if res == "written+verified":
            n_flip += plan["action"] == "flip_z"
            n_copy += plan["action"] == "identity"
        else:
            n_fail += 1

    verb = "written" if a.apply else "to write (DRY RUN)"
    print(f"\n{n_flip} flip_z + {n_copy} identity {verb} | "
          f"{n_skip} skipped | {n_fail} problems")
    if not a.apply and (n_flip or n_copy):
        print("re-run with --apply to commit")
    if a.apply and not n_fail:
        print(f"\nNEXT: point step14's gradient discovery at '{FIX_STEM}' for "
              f"ALL subjects and refit with --overwrite; then re-run step20d "
              f"on the refit tensors -- every subject should now win 'none'.")


if __name__ == "__main__":
    main()
