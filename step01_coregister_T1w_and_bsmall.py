#!/usr/bin/env python3
"""
Register the first DKI volume (b-small) of every participant to the T1w with
ANTs (T1w fixed, b0 moving). Writes <sub>_ses-01_T1wbrain_restore_from_b0_*
(0GenericAffine.mat, 1Warp, 1InverseWarp, Warped) in <sub>/DKI_in_T1w_space/.
No arguments: every participant under data_root is processed.

Usage:  python step01_coregister_T1w_and_bsmall.py
"""
import argparse
import os
import glob
import subprocess
import nibabel as nib
from pipeline_config import cfg

# ---- CONFIG ----
DKI_ROOT = cfg("data_root")
T1_ROOT  = cfg("pvs_root")

# Maximum depth (0 = root only, 1 = one level down, etc.)
MAX_DEPTH = 5

# ANTs registration parameters (adjust if you like)
ANTS_CMD = "antsRegistration"
FSLMATHS_CMD = "fslmaths"

def find_dki_files(root):
    """
    Walk DKI_ROOT up to MAX_DEPTH and yield paths to 4D NIfTI files
    whose name contains '_eddy.nii.gz'.
    """
    for dirpath, dirnames, filenames in os.walk(root):
        # Limit depth
        rel = os.path.relpath(dirpath, root)
        depth = 0 if rel == "." else rel.count(os.sep) + 1
        if depth > MAX_DEPTH:
            # Don't descend further
            dirnames[:] = []
            continue

        for f in filenames:
          if "DKI" in f:
            if f.endswith("_eddy.nii.gz"):
               if "DATASINK" in dirpath:
                yield os.path.join(dirpath, f)

def extract_subject_id_from_dki(dki_path):
    """
    Heuristic: subject ID = basename before '_eddy'.
    e.g. Sub-105_eddy.nii.gz -> 'Sub-105'
    """
    #base = os.path.basename(dki_path)
    #subj_part = base.split("_eddy")[0]
    subj_part = dki_path.split("/")[-6]
    return subj_part

def find_t1_for_subject(subj_id, t1_root):
    """
    Search recursively under T1_ROOT for a file that:
    - contains subj_id
    - contains '_T1wbrain_restore.nii.gz'
    Return first match or None.
    """
    pattern = os.path.join(t1_root, f"**/{subj_id}*T1wbrain_restore.nii.gz")
    matches = glob.glob(pattern, recursive=True)
    if not matches:
        return None
    # If more than one, just pick the first
    return matches[0]

def extract_first_volume(dki_path):
    """
    Load 4D DKI NIfTI, extract first volume, save as 3D NIfTI
    in the same folder with suffix '_b0.nii.gz'.
    Returns path to new 3D file.
    """
    img = nib.load(dki_path)
    data = img.get_fdata()
    if data.ndim != 4:
        raise ValueError(f"{dki_path} is not 4D (shape: {data.shape})")

    first_vol = data[..., 0]
    new_img = nib.Nifti1Image(first_vol, img.affine, img.header)
    # Make sure header knows it's 3D
    new_img.header.set_data_shape(first_vol.shape)

    base = os.path.basename(dki_path)
    base_noext = base.replace(".nii.gz", "")
    out_name = base_noext + "_b0.nii.gz"
    out_path = os.path.join(os.path.dirname(dki_path), out_name)
    nib.save(new_img, out_path)
    return out_path

def run_ants_registration(moving_b0, fixed_t1):
    """
    Run antsRegistration:
    - fixed  = T1w
    - moving = b0
    Output warped image in the T1 folder.
    """
    t1_dir = os.path.dirname(fixed_t1)
    t1_base = os.path.basename(fixed_t1).replace(".nii.gz", "")
    warped_dir = t1_dir.replace("Segmentation", "DKI_in_T1w_space")
    os.makedirs(warped_dir, exist_ok=True)
    prefix = os.path.join(warped_dir, f"{t1_base}_from_b0_")
    reg_mask = fixed_t1.replace("T1wbrain_restore", "PVSROI")
    reg_mask_bin = reg_mask.replace(".nii", "_binarized.nii")

    #cmd = [FSLMATHS_CMD, reg_mask, "-thr", "0.5", "-bin", reg_mask_bin]
    #subprocess.run(cmd, check=True)

    cmd = [
    ANTS_CMD,
    "--dimensionality", "3",
    "--float", "0",
    "--output", f"[{prefix},{prefix}Warped.nii.gz]",
    "--interpolation", "Linear",
    "--use-histogram-matching", "1",
    "--winsorize-image-intensities", "[0.005,0.995]",
    # Initial transform (center-of-mass)
    "--initial-moving-transform", f"[{fixed_t1},{moving_b0},1]",

    # --- Rigid ---
    "--transform", "Rigid[0.1]",
    "--metric", f"Mattes[{fixed_t1},{moving_b0},1,64,Regular,0.2]",
    "--convergence", "[1000x500x250,1e-6,10]",
    "--shrink-factors", "4x2x1",
    "--smoothing-sigmas", "2x1x0vox",

    # --- Affine ---
    "--transform", "Affine[0.1]",
    "--metric", f"Mattes[{fixed_t1},{moving_b0},1,64,Regular,0.2]",
    "--convergence", "[1000x500x250,1e-6,10]",
    "--shrink-factors", "4x2x1",
    "--smoothing-sigmas", "2x1x0vox",

    # --- Mild SyN non-linear refinement ---
    "--transform", "SyN[0.1,3,0]",
    "--metric", f"CC[{fixed_t1},{moving_b0},1,4]",
    "--convergence", "[70x50x20,1e-6,10]",
    "--shrink-factors", "4x2x1",
    "--smoothing-sigmas", "2x1x0vox",
]


    print("\nRunning ANTs registration:")
    print("  Fixed (T1): ", fixed_t1)
    print("  Moving (b0):", moving_b0)
    print("  Cmd:", " ".join(cmd))

    subprocess.run(cmd, check=True)

    warped_image = f"{prefix}Warped.nii.gz"
    print("  -> Warped image saved as:", warped_image)
    return warped_image

def main():
    dki_files = list(find_dki_files(DKI_ROOT))
    if not dki_files:
        print("No _eddy.nii.gz files found.")
        return

    print(f"Found {len(dki_files)} DKI files.")

    for dki_path in dki_files:
        print("\n=== Processing ===")
        print("DKI:", dki_path)

        subj_id = extract_subject_id_from_dki(dki_path)
        print("Subject ID:", subj_id)

        t1_path = find_t1_for_subject(subj_id, T1_ROOT)
        if t1_path is None:
            print("  !! No T1w file found for this subject. Skipping.")
            continue

        print("T1w:", t1_path)

        # 1) Extract first b-volume
        try:
            b0_path = extract_first_volume(dki_path)
        except Exception as e:
            print(f"  !! Failed to extract first volume: {e}")
            continue

        print("b0 extracted:", b0_path)

        # 2) Run ANTs registration (T1 as fixed, b0 as moving)
        try:
            run_ants_registration(b0_path, t1_path)
        except subprocess.CalledProcessError as e:
            print(f"  !! ANTs registration failed: {e}")
            continue

if __name__ == "__main__":
    argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter).parse_args()
    main()
