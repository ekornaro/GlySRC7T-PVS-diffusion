"""
step22_create_csv_based_on_niftis_and_SPM12ICV.py
(from create_csv_based_on_niftis_and_SPM12ICV_2026-09-29_median.py)

CHANGELOG 2026-10-01
--------------------
 * --subjects sub-NNN [sub-MMM ...]: process ONLY these subjects and MERGE
   their rows into the existing spreadsheets: each file is read, any old rows
   of these subjects are dropped, the new rows are added, and every other
   subject's row is kept unchanged. Without --subjects the script behaves as
   before (all subjects, files rewritten).
   Safeguards: a requested subject must be marked 'Yes' in the study xlsx and
   have a folder on disk; an existing file must contain all three group sheets
   (a 'ppcs'-only file from the pre-fix script stops the run) and the same
   columns; every file is backed up as <file>.xlsx.bak_<timestamp> before it
   is rewritten (the .bak suffix keeps it out of the stats scripts' scans); a
   missing file is never created from the selected subjects alone.
   The summary of a merge run is creation_summary_<mode>_merge-<subjects>.csv
   and reports the per-group counts AFTER the merge against the full cohort.
 * Command line: python <script> [vol|diff|both] [--subjects ...]

CHANGELOG 2026-09-29
--------------------

CHANGELOG 2026-09-29 (details)
------------------------------
 * SPM12 TIV WAS READ FROM ONE SHEET ONLY. TIV_spm12_results.xlsx has one
   sheet per group (PPCS, HC, Boxers, ...). pd.read_excel() without
   sheet_name reads only the first (PPCS), so every HC and RHI subject hit
   "No TIV found ... skipping" and was dropped from EVERY spreadsheet, in both
   volumetry and diffusion mode (the lookup precedes the mode branch). Files
   written with the old script contain a 'ppcs' sheet only.
   Now: HC, PPCS and Boxers sheets are all read, subjects are matched on the
   exact ID (sub-NNN) instead of a substring search, and a missing TIV is a
   hard error at START-UP, before any image is loaded.
 * Subject lists use the same inclusion rule as the other rebuttal scripts
   ('Will we used it?' == Yes); blank and 'Total:' rows no longer become
   'sub-nan' / 'sub-Total:'.
 * Every spreadsheet is checked after it is built: the number of subjects per
   group is printed when it differs from the group's subject count, a file with
   no rows at all is not written (openpyxl cannot save a workbook without a
   sheet), and creation_summary_<mode>.csv lists n per group for every file,
   so a missing group can no longer pass unnoticed.
"""
import os
import shutil
import argparse
import datetime
import pandas as pd
import nibabel as nib
import numpy as np
import warnings
from pipeline_config import cfg
warnings.filterwarnings("ignore")

# Define spatial resolution (voxel volume in mm^3) in advance as a global variable
voxel_volume_mm3 = 0.512 # for 0.8x0.8x0.8 mm3 voxels

def get_subjects(sheet_name, col="Will we used it?"):
    """Subjects marked 'Yes' in the study xlsx (same rule as the other rebuttal
    scripts). Blank and 'Total:' rows are skipped."""
    df = pd.read_excel(info_lund_7t, sheet_name=sheet_name)
    if col not in df.columns:
        raise SystemExit(f"'{col}' missing in sheet {sheet_name} of {info_lund_7t}")
    df = df[df[col].astype(str).str.strip().str.lower() == "yes"]
    out = set()
    for v in df['Subj'].dropna():
        s = str(v).strip()
        if s.replace(".0", "").isdigit():
            out.add(f"sub-{int(float(s))}")
    return out

def get_region_voxel_counts_from_data(parcel_data):
    """
    Given a flattened parcellation map, returns a dictionary mapping each region label (as int)
    to its voxel count.
    """
    labels, counts = np.unique(parcel_data, return_counts=True)
    return {int(round(label)): count for label, count in zip(labels, counts)}

def volumetry_process_subject(sub, base_dir, label_name, tissue_size_in_voxels,region_mask, total_brain_vol,voxel_volume=voxel_volume_mm3):
    ...
    # Compute "global" (WHOLE) averages:
    
    subID = sub.split("-")[-1]
    rows = []
    # Construct file paths
    pvs_path    = os.path.join(base_dir, sub, "PVS", f"{sub}_ses-01_PVS_postprocessed.nii.gz")
    # Check existence of required files
    if not (os.path.exists(pvs_path)):
        return rows
    try:
        # Load images once per subject and flatten
        pvs_data = nib.load(pvs_path).get_fdata().ravel()
    except Exception as e:
        print(f"Error loading images for {sub}: {e}")
        return rows
    combined_mask = (pvs_data > 0) & region_mask
    if not isinstance(tissue_size_in_voxels,np.int64):
        print(tissue_size_in_voxels)
        print(type(tissue_size_in_voxels))
        print(sub)
        print(label_name)
    tissue_size_in_mm = tissue_size_in_voxels * voxel_volume
    relative_tissue_size_in_mm = tissue_size_in_mm/total_brain_vol
    num_pvs_voxels = np.sum(pvs_data[combined_mask] > 0)
    volume_pvs_voxels = np.sum(pvs_data[combined_mask] > 0) * voxel_volume
    relative_volume_pvs_voxels = volume_pvs_voxels/tissue_size_in_mm
    rows.append([sub, label_name, tissue_size_in_voxels, tissue_size_in_mm,relative_tissue_size_in_mm, num_pvs_voxels,volume_pvs_voxels,relative_volume_pvs_voxels])
    return rows

def volumetry_process_subject_eTIV(sub, base_dir, label_name, total_brain_vol,voxel_volume=voxel_volume_mm3):
    ...
    # Compute "global" (WHOLE) averages:
    
    subID = sub.split("-")[-1]
    rows = []
    # Construct file paths
    pvs_path    = os.path.join(base_dir, sub, "PVS", f"{sub}_ses-01_PVS_postprocessed.nii.gz")
    # Check existence of required files
    if not (os.path.exists(pvs_path)):
        return rows
    try:
        # Load images once per subject and flatten
        pvs_data = nib.load(pvs_path).get_fdata().ravel()
    except Exception as e:
        print(f"Error loading images for {sub}: {e}")
        return rows
    tissue_size_in_voxels = total_brain_vol / voxel_volume
    tissue_size_in_mm = total_brain_vol
    relative_tissue_size_in_mm = tissue_size_in_mm/total_brain_vol
    num_pvs_voxels = np.sum(pvs_data > 0)
    volume_pvs_voxels = np.sum(pvs_data > 0) * voxel_volume
    relative_volume_pvs_voxels = volume_pvs_voxels/total_brain_vol
    rows.append([sub, label_name, tissue_size_in_voxels, tissue_size_in_mm,relative_tissue_size_in_mm, num_pvs_voxels,volume_pvs_voxels,relative_volume_pvs_voxels])
    return rows

def volumetry_process_subject_csf_fast(sub, base_dir, label_name, parcel_data, total_brain_vol,voxel_volume=voxel_volume_mm3):
    ...
    # Compute "global" (WHOLE) averages: 
    subID = sub.split("-")[-1]
    rows = []
    # Construct file paths
    pvs_path    = os.path.join(base_dir, sub, "PVS", f"{sub}_ses-01_PVS_postprocessed.nii.gz")
    csv_path    = os.path.join(base_dir, sub, "FreeSurfer", "synthseg", "volumetry_scores.csv")
    # Check existence of required files
    if not (os.path.exists(pvs_path)):
        return rows
    try:
        # Load images once per subject and flatten
        pvs_data = nib.load(pvs_path).get_fdata().ravel()
    except Exception as e:
        print(f"Error loading images for {sub}: {e}")
        return rows
    tissue_size_voxels = np.sum(parcel_data>0.5)
    tissue_size_mm = np.sum(parcel_data>0.5) * voxel_volume
    tissue_relative_size_mm = tissue_size_mm/total_brain_vol
    combined_mask = (pvs_data > 0) & (parcel_data > 0.5)
    pvs_in_csf_voxels = np.sum(pvs_data[combined_mask]>0)
    pvs_in_csf_mm = np.sum(pvs_data[combined_mask]>0)*voxel_volume
    pvs_in_csf_relative_mm = pvs_in_csf_mm/tissue_size_mm
    rows.append([sub, "CSF_fast", tissue_size_voxels, tissue_size_mm,tissue_relative_size_mm, pvs_in_csf_voxels,pvs_in_csf_mm,pvs_in_csf_relative_mm])
    return rows

def diffusion_process_subject(sub, base_dir, diffusion_map, label_name, region_mask):
    ...
    # Compute "global" (WHOLE) averages:
    
    subID = sub.split("-")[-1]
    rows = []
    # Construct file paths
    pvs_path    = os.path.join(base_dir, sub, "PVS", f"{sub}_ses-01_PVS_postprocessed.nii.gz")
    if diffusion_map == "MK" or diffusion_map == "AK" or diffusion_map == "RK":
        map_path    = os.path.join(base_dir, sub, "DKI_in_T1w_space", f"{diffusion_map}_medfilt_toT1w_Warped.nii.gz")
    else:
        map_path    = os.path.join(base_dir, sub, "DKI_in_T1w_space", f"{diffusion_map}_toT1w_Warped.nii.gz")
    csv_path    = os.path.join(base_dir, sub, "FreeSurfer", "synthseg", "volumetry_scores.csv")
    # Check existence of required files
    if not (os.path.exists(pvs_path) and os.path.exists(map_path)):
        return rows
    try:
        # Load images once per subject and flatten
        pvs_data = nib.load(pvs_path).get_fdata().ravel()
        map_data = nib.load(map_path).get_fdata().ravel()
    except Exception as e:
        print(f"Error loading images for {sub}: {e}")
        return rows
    tissue_size = np.count_nonzero(region_mask)
    pvs_size = np.count_nonzero(pvs_data)
    if label_name == "eTIV":
        map_in_region = _med(map_data, diffusion_map)
        pvs_vals      = map_data[pvs_data > 0]
        map_in_pvs    = _med(pvs_vals, diffusion_map)
        dex_tiss      = (np.sum(pvs_vals[_valid(pvs_vals, diffusion_map)])/tissue_size
                         if np.any(pvs_data > 0) else np.nan)
        pvs_weighted_map = map_in_pvs*(pvs_size/tissue_size) if np.any(pvs_data > 0) else np.nan
        reg_med       = _med(map_data[region_mask], diffusion_map)
        dex_ratio     = ((map_in_pvs/reg_med)*(pvs_size/tissue_size)
                         if (np.any(pvs_data > 0) and np.isfinite(reg_med) and reg_med != 0) else np.nan)
    else:
        map_region = map_data[region_mask]
        map_in_region = _med(map_region, diffusion_map)
        combined_mask = (pvs_data > 0) & region_mask
        pvs_size_within_ROI = np.count_nonzero(combined_mask)
        map_pvs = map_data[combined_mask]
        map_in_pvs = _med(map_pvs, diffusion_map)
        dex_tiss = map_in_pvs/tissue_size if np.isfinite(map_in_pvs) else np.nan
        pvs_weighted_map = map_in_pvs*(pvs_size_within_ROI/tissue_size) if np.isfinite(map_in_pvs) else np.nan
        dex_ratio = ((map_in_pvs/map_in_region)*(pvs_size_within_ROI/tissue_size)
                     if (np.isfinite(map_in_pvs) and np.isfinite(map_in_region) and map_in_region != 0) else np.nan)
        
    rows.append([sub, diffusion_map, label_name, map_in_region, map_in_pvs, dex_tiss,pvs_weighted_map,dex_ratio])
    return rows

def diffusion_process_subject_csf_fast(sub, base_dir, diffusion_map, label_name, parcel_data):
    ...
    # Compute "global" (WHOLE) averages:
    
    subID = sub.split("-")[-1]
    rows = []
    # Construct file paths
    pvs_path    = os.path.join(base_dir, sub, "PVS", f"{sub}_ses-01_PVS_postprocessed.nii.gz")
    if diffusion_map == "MK" or diffusion_map == "AK" or diffusion_map == "RK":
        map_path    = os.path.join(base_dir, sub, "DKI_in_T1w_space", f"{diffusion_map}_medfilt_toT1w_Warped.nii.gz")
    else:
        map_path    = os.path.join(base_dir, sub, "DKI_in_T1w_space", f"{diffusion_map}_toT1w_Warped.nii.gz")
    csv_path    = os.path.join(base_dir, sub, "FreeSurfer", "synthseg", "volumetry_scores.csv")
    # Check existence of required files
    if not (os.path.exists(pvs_path) and os.path.exists(map_path)):
        return rows
    try:
        # Load images once per subject and flatten
        pvs_data = nib.load(pvs_path).get_fdata().ravel()
        map_data = nib.load(map_path).get_fdata().ravel()
    except Exception as e:
        print(f"Error loading images for {sub}: {e}")
        return rows
    map_region = map_data[parcel_data>0.5]
    map_in_region = _med(map_region, diffusion_map)
    combined_mask = (pvs_data > 0) & (parcel_data>0.5)
    map_in_pvs = _med(map_data[combined_mask], diffusion_map)
        
    rows.append([sub, diffusion_map, label_name, map_in_region, map_in_pvs,0,0,0])
    return rows

################################################################################################################################################################
# Load synthseg labels
synthseg_sprd = cfg("synthseg_labels_xlsx")
df_labels = pd.read_excel(synthseg_sprd, sheet_name="Labels")
labels_names_synthseg  = df_labels.iloc[:, 0].tolist()
labels_values_synthseg = df_labels.iloc[:, 1].tolist()
eTIV_label = 10000 # this is an arbitary value I gave to the eTIV segmentation - it is not included in the official synthseg LUT
CSF_fast_label = 240 # this is an arbitary value I gave to the CSF segmentation derived by the CSF-fast algorithm - it is not included in the official synthseg LUT
ppcs_hc_to_boxers_lim = 148 # this is the ID of the last "non-boxer" subject we scanned - boxers' ID start from 149 and further
pvs_cases = ["conservative"]#conservative, ambitious
# =====================================================================
# WHAT TO EXTRACT  --  choose ONE mode, on the command line or here.
# =====================================================================
#   python <script>            -> diffusion metrics (default)
#   python <script> vol        -> PVS volumetry
#   python <script> diff       -> diffusion metrics (explicit)
#   python <script> both       -> volumetry first, then diffusion
# The two lists are kept side by side so neither has to be commented
# out (the previous edit-one-comment-the-other pattern is exactly how a
# run ends up using the wrong list without anyone noticing).

# Every scalar *_toT1w_Warped map present in <sub>/DKI_in_T1w_space/ after the
# step15_v2 uniform re-warp (2026-09-01) + step16 v3 orientation maps.
# EXCLUDED on purpose: Dxx..Dyz (tensor components, not metrics), V1 (4-D
# eigenvector field), and the FWE family (FAt/MDt/ADt/RDt/MKt/AKt/RKt/f) --
# never warped, and excluded by the step23 reliability gate.
# theta / clusterLinearity / distToPVS are extracted only if you add them
# here; they are geometry/QC descriptors, NOT hypothesis-test outcomes, and
# were removed from the tested family on purpose.
DIFFUSION_METRICS = [
    # standard DKI (all shells)
    "FA", "MD", "AD", "RD", "MK_mf", "AK_mf", "RK_mf",
    # low-b (b <= 500) maps, step13 v2
    "lowb_ADC", "lowb_MD", "lowb_AD", "lowb_RD", "lowb_FA",
    # orientation-referenced maps, step16 v3 (header-corrected frame)
    "Dpar", "Dperp", "Dpar_minus_Dperp",
]

# "VOL" is a MODE, not a map: it routes to volumetry_process_subject* and
# writes PVS counts/volumes instead of median map values.
VOLUMETRY = ["VOL"]

_MODES = {"vol": VOLUMETRY, "volumetry": VOLUMETRY,
          "diff": DIFFUSION_METRICS, "diffusion": DIFFUSION_METRICS,
          "both": VOLUMETRY + DIFFUSION_METRICS}
_ap = argparse.ArgumentParser(description="Per-subject volumetry / diffusion spreadsheets")
_ap.add_argument("mode", nargs="?", default="diff",
                 help="vol | diff | both (default: diff)")
_ap.add_argument("--subjects", nargs="+", default=None, metavar="sub-NNN",
                 help="process only these subjects and merge them into the existing "
                      "spreadsheets (all other rows are kept)")
_args = _ap.parse_args()
_arg = _args.mode.lower()
if _arg not in _MODES:
    raise SystemExit(f"unknown mode '{_arg}'. Use one of: {sorted(_MODES)}")
SUBSET = None
if _args.subjects:
    SUBSET = sorted({s if s.startswith("sub-") else f"sub-{s}" for s in _args.subjects})
_STAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
ccases = _MODES[_arg]
EXTRACTION_MODE = _arg
print(f"[mode] {_arg}: {len(ccases)} case(s) -> "
      f"{ccases if _arg.startswith(('vol', 'both')) else 'diffusion metrics'}")
if SUBSET:
    print(f"[mode] MERGE: only {SUBSET}; existing spreadsheets are updated, other rows kept")

SIGNED_METRICS = {"Dpar_minus_Dperp", "theta"}
SUMMARY = []          # one dict per spreadsheet: n per group (written at the end)


def _valid(vals, metric):
    """Voxels that carry a real value for this metric.
    Non-negative metrics (diffusivities, kurtosis, FA, distances): value > 0,
    which also drops the 0 background fill and NaN.
    Signed metrics (Dpar_minus_Dperp, theta): any finite value -- negative
    Dpar_minus_Dperp is physically meaningful and must NOT be discarded."""
    vals = np.asarray(vals)
    if metric in SIGNED_METRICS:
        return np.isfinite(vals)
    return vals > 0


def _med(vals, metric):
    v = np.asarray(vals)[_valid(vals, metric)]
    return float(np.median(v)) if v.size else np.nan
#segmentations = ["CSF_fast"]
#segmentations = ["synthseg"]
segmentations = ["synthseg","scPVSROI","PVSROI"]
# Load subject info for each group
info_lund_7t = cfg("study_xlsx")

# --- self-check: does the metric list match what is on disk? -----------------
def _verify_ccases_against_disk(base_dir, subjects, ccases):
    import glob as _glob
    for sub in subjects:
        d = os.path.join(base_dir, sub, "DKI_in_T1w_space")
        if not os.path.isdir(d):
            continue
        on_disk = {os.path.basename(p).replace("_toT1w_Warped.nii.gz", "")
                   for p in _glob.glob(os.path.join(d, "*_toT1w_Warped.nii.gz"))}
        tensor = {"Dxx", "Dyy", "Dzz", "Dxy", "Dxz", "Dyz", "V1", "V2", "V3"}
        missing = [c for c in ccases if c not in on_disk]
        unlisted = sorted(on_disk - set(ccases) - tensor)
        print(f"[ccases check on {sub}] {len(ccases)} listed; "
              f"missing on disk: {missing or 'none'}; "
              f"on disk but not listed: {unlisted or 'none'}")
        return
    print("[ccases check] no subject with a DKI_in_T1w_space folder found")

hc_subjects     = get_subjects('HC')
ppcs_subjects   = get_subjects('PPCS')
boxers_subjects = get_subjects('Boxers')
#severe_subjects = get_subjects('severeTBI')

# Directories for imaging data
SMP12_results_spreadsheet = cfg("tiv_spm12_xlsx")
TIV_SHEETS = ("HC", "PPCS", "Boxers")


def load_tiv_mm3(path, sheets=TIV_SHEETS):
    """{'sub-NNN': TIV in mm^3} from ALL group sheets (the file has one sheet
    per group; reading only the first silently dropped HC and RHI)."""
    try:
        xls = pd.ExcelFile(path)
    except Exception as e:
        raise SystemExit(f"cannot read TIV spreadsheet {path}: {e}")
    missing_sheets = [sh for sh in sheets if sh not in xls.sheet_names]
    if missing_sheets:
        raise SystemExit(f"{path}: sheets {missing_sheets} not found "
                         f"(present: {xls.sheet_names})")
    tiv = {}
    for sh in sheets:
        df = xls.parse(sh)
        if "subject" not in df.columns or "TIV_ml" not in df.columns:
            raise SystemExit(f"{path} [{sh}]: need columns 'subject' and 'TIV_ml', "
                             f"found {list(df.columns)[:6]}")
        df = df.dropna(subset=["subject", "TIV_ml"])
        ids = df["subject"].astype(str).str.extract(r"^(sub-\d+)")[0]
        for sid, ml in zip(ids, df["TIV_ml"]):
            if pd.isna(sid):
                continue
            if sid in tiv and not np.isclose(tiv[sid], ml * 1000.0):
                raise SystemExit(f"{path}: {sid} has two different TIV values")
            tiv[sid] = float(ml) * 1000.0            # mL -> mm^3
    print(f"[TIV] {len(tiv)} subjects read from sheets {list(sheets)} of {path}")
    return tiv


TIV_MM3 = load_tiv_mm3(SMP12_results_spreadsheet)
ppcs_boxers_hc_dir_init = cfg("pvs_segmentations_dir")
for pvscc,pvsccase in enumerate(pvs_cases):   
 # Directories and file paths
 # Volumetry and diffusion write to SEPARATE subtrees so a mode switch can
 # never overwrite the other mode's spreadsheets. Tree and study-xlsx label
 # follow the 2026-09 rebuttal outputs (the old path pointed at .../results/
 # with a 2025-08-28 label, i.e. last year's tree).
 output_xlsx_dir1_tmpl = os.path.join(cfg("results_dir"),
                     "spreadsheets/tissues-PVS_VOL-DIFF_nifti-based_and_SPM12-eTIV-based/"
                     "%s/{subtree}/based_on_basic_info_created_by_Evgenios_updated_on_"
                     "2026-08-28_PVS_segmentation_xlsx" % pvsccase)
 # (created per case below, once the subtree is known)
 
 ppcs_boxers_hc_dir = os.path.join(ppcs_boxers_hc_dir_init,"Lund_7T_%s"%pvsccase,"Rob_bfc","using_FLAIR")

 # Get available subject folders for each group.
 all_ppcs_boxers_hc = sorted(os.listdir(ppcs_boxers_hc_dir))

 # Filter subjects for each group based on available processed folders.
 subjects_hc     = sorted(sub for sub in all_ppcs_boxers_hc if sub in hc_subjects)
 subjects_ppcs   = sorted(sub for sub in all_ppcs_boxers_hc if sub in ppcs_subjects)
 #subjects_severe = sorted(sub for sub in all_ppcs_boxers_hc if sub in severe_subjects)
 subjects_boxers = sorted(sub for sub in all_ppcs_boxers_hc if sub in boxers_subjects)
 # For each group, process subjects, filter for the current tissue label, and write a sheet.
 groups = {
             "hc": (subjects_hc, ppcs_boxers_hc_dir),
             "ppcs": (subjects_ppcs, ppcs_boxers_hc_dir),
             "boxers": (subjects_boxers, ppcs_boxers_hc_dir),
             #"severe-tbi": (subjects_severe, ppcs_boxers_hc_dir, ppcs_hc_dir_proc)
        }
 # full cohort sizes (the expected counts, also in a merge run)
 FULL_N = {"hc": len(subjects_hc), "ppcs": len(subjects_ppcs), "boxers": len(subjects_boxers)}
 if SUBSET:
     _known = set(subjects_hc) | set(subjects_ppcs) | set(subjects_boxers)
     _bad = [s for s in SUBSET if s not in _known]
     if _bad:
         raise SystemExit(f"--subjects {_bad}: not marked 'Yes' in {info_lund_7t} "
                          f"or no folder in {ppcs_boxers_hc_dir}")
     subjects_hc = [s for s in subjects_hc if s in SUBSET]
     subjects_ppcs = [s for s in subjects_ppcs if s in SUBSET]
     subjects_boxers = [s for s in subjects_boxers if s in SUBSET]
     groups = {"hc": (subjects_hc, ppcs_boxers_hc_dir),
               "ppcs": (subjects_ppcs, ppcs_boxers_hc_dir),
               "boxers": (subjects_boxers, ppcs_boxers_hc_dir)}
     print("[subjects] merge run: " + ", ".join(
         f"{g}: {groups[g][0]}" for g in groups if groups[g][0]))
 # TIV must exist for every subject that will be processed -- fail NOW, not
 # after hours of NIfTI loading (the old script skipped them silently)
 _no_tiv = sorted(sub for grp in (subjects_hc, subjects_ppcs, subjects_boxers)
                  for sub in grp if sub not in TIV_MM3)
 if _no_tiv:
     raise SystemExit(f"No SPM12 TIV for {_no_tiv} in {SMP12_results_spreadsheet}")
 print(f"[subjects] hc={len(subjects_hc)} ppcs={len(subjects_ppcs)} "
       f"boxers={len(subjects_boxers)} (study xlsx 'Yes' AND folder on disk); "
       f"TIV found for all")
 # one-time check that ccases matches the maps actually on disk (first subject)
 if "_ccases_checked" not in globals():
     _verify_ccases_against_disk(ppcs_boxers_hc_dir,
                                 list(subjects_hc) + list(subjects_ppcs) + list(subjects_boxers),
                                 [c for c in ccases if c != "VOL"])   # VOL is a mode, not a map
     _ccases_checked = True
 # We'll use only one diffusion metric here (e.g., "FA")
 for cc,ccase in enumerate(ccases):
  # For each tissue label (from synthseg labels), create a single Excel file with four sheets.
  for se,seg in enumerate(segmentations):
   if seg != "CSF_fast":
      labels_values = labels_values_synthseg[:-1]
      labels_names = labels_names_synthseg[:-1]
   else:
      labels_values = [labels_values_synthseg[-1]]
      labels_names = [labels_names_synthseg[-1]]
   for label_value, label_name in zip(labels_values, labels_names): 
     if ccase != "VOL" and label_name == "eTIV":
       continue 
    #if label_name == "right cerebellum cortex":   
     toProceed = False
     for group_name, (subject_list, base_dir) in groups.items():  
      for sub in subject_list:
                segnow_path = os.path.join(base_dir, sub, "Segmentation", f"{sub}_ses-01_%s.nii.gz"%seg) 
                parcel_data = nib.load(segnow_path).get_fdata().ravel()
                synthseg_regions = get_region_voxel_counts_from_data(parcel_data)
                try:
                    if "background" not in label_name and "3rd" not in label_name and "4th" not in label_name and "brain-stem" not in label_name and "CSF" not in label_name and "midbrain" not in label_name and "eTIV" not in label_name:
                        label_vals = label_value.split(" ")
                        label_val_int = int(label_vals[0]) 
                    else:
                        label_val_int = int(label_value)
                except Exception:
                    continue 
                if int(label_val_int) in [int(k) for k in synthseg_regions.keys()] or label_val_int == eTIV_label or label_val_int == CSF_fast_label:  
                    toProceed = True
                    continue
      if toProceed:
         continue
     if toProceed:             
      label_namef = label_name.replace(" ","_") 
      if label_namef == "CSF_fast" and seg != "CSF_fast":
         continue
      output_xlsx_dir1 = output_xlsx_dir1_tmpl.format(
          subtree="volumetry" if ccase == "VOL" else "diffusion")
      os.makedirs(output_xlsx_dir1, exist_ok=True)
      output_xlsx_dir = os.path.join(output_xlsx_dir1,label_namef,seg)
      os.makedirs(output_xlsx_dir, exist_ok=True)
      output_file = os.path.join(output_xlsx_dir, f"{ccase}.xlsx")
      print("Processing case %s on %s by %s... "%(ccase,label_name,seg))   
      group_frames = {}                                   # written below
      if True:
         for group_name, (subject_list, base_dir) in groups.items():
            rows = [] 
            for sub in subject_list:
                #print("Processing case %s on %s by %s for subject %s... "%(ccase,label_name,seg,sub))   
                subID = sub.split("-")[-1]
                # exact-ID lookup; presence was verified at start-up
                total_brain_vol = TIV_MM3[sub]                # mm^3
                segnow_path = os.path.join(base_dir, sub, "Segmentation", f"{sub}_ses-01_%s.nii.gz"%seg) 
                parcel_data = nib.load(segnow_path).get_fdata().ravel()
                # Get region counts from the parcellation
                synthseg_regions = get_region_voxel_counts_from_data(parcel_data)
                try: 
                    if "background" not in label_name and "3rd" not in label_name and "4th" not in label_name and "brain-stem" not in label_name and "CSF" not in label_name and "midbrain" not in label_name and "eTIV" not in label_name:
                        label_vals = label_value.split(" ")
                        label_val_int = [int(label_vals[0]), int(label_vals[1])]
                    else:
                        label_val_int = [int(label_value)]
                except Exception:
                    continue
                if label_val_int[0] not in synthseg_regions and label_val_int[0] != eTIV_label and label_val_int[0] != CSF_fast_label: 
                    continue
                elif label_val_int[0] == eTIV_label or label_val_int[0] == CSF_fast_label:
                    size_of_label = []
                else: 
                    if "background" not in label_name and "3rd" not in label_name and "4th" not in label_name and "brain-stem" not in label_name and "CSF" not in label_name and "midbrain" not in label_name and "eTIV" not in label_name:
                        size_of_label = synthseg_regions[label_val_int[0]] + synthseg_regions[label_val_int[1]]
                    else:
                        size_of_label = synthseg_regions[label_val_int[0]]
                # Create a boolean mask for the current tissue
                if "background" not in label_name and "3rd" not in label_name and "4th" not in label_name and "brain-stem" not in label_name and "CSF" not in label_name and "midbrain" not in label_name and "eTIV" not in label_name:
                    region_mask = np.isin(parcel_data, [label_val_int[0],label_val_int[1]])
                else:
                    region_mask = (parcel_data == label_val_int[0])
                if not np.any(region_mask) and label_val_int[0] != eTIV_label and label_val_int[0] != CSF_fast_label:
                    continue
                # Process the subject and then filter rows for the current tissue label.
                #print(f"Creating spreadsheet for tissue: {label_name} and case: {ccase} and subject: {sub}")
                if ccase == "VOL":
                    if seg == "CSF_fast": 
                        subject_rows = volumetry_process_subject_csf_fast(sub, base_dir, label_name,parcel_data,total_brain_vol)
                    else:
                        if label_namef == "eTIV":
                            subject_rows = volumetry_process_subject_eTIV(sub, base_dir,label_name, total_brain_vol)
                        else:
                            subject_rows = volumetry_process_subject(sub, base_dir,label_name, size_of_label, region_mask,total_brain_vol)
                else:
                    if seg != "CSF_fast":
                        subject_rows = diffusion_process_subject(sub, base_dir, ccase, label_name, region_mask)
                    else:
                        subject_rows = diffusion_process_subject_csf_fast(sub, base_dir, ccase, label_name, parcel_data)
                # Only take rows where the tissue-label matches the current label (skip the "WHOLE" row).
                for r in subject_rows:
                    if ccase != "VOL":
                        to_check = r[2]
                    else:
                        if seg != "CSF_fast":
                            to_check = r[1]  
                        else:
                            to_check = r[2]
                    if to_check == label_name:
                        rows.append(r)            
                if ccase == "VOL":
                    df_group = pd.DataFrame(rows, columns=["subject", "tissue-label", "tissue-size-in-voxels","tissue-size-in-mm", "relative-tissue-size-in-mm","pvs-size-within-tissue-in-voxels","pvs-size-within-tissue-in-mm","relative-pvs-size-within-tissue-in-mm"])   
                else: 
                    df_group = pd.DataFrame(rows, columns=["subject", "diffusion-metric", "tissue-label", "diffusion_within_tissue", "diffusion_within_pvs","Dex-tissue","PVS-Weighted-diffusion","Dex-ratio"])
            if len(rows)>0:
                group_frames[group_name] = df_group
      # ---- merge run: fold the selected subjects into the existing file ----
      if SUBSET:
          if not group_frames:
              print(f"  [i] nothing to add for {SUBSET} -> {output_file} unchanged")
              continue
          if not os.path.exists(output_file):
              print(f"  [!] {output_file} does not exist -> NOT created from {SUBSET} "
                    f"alone (run the full script first)")
              SUMMARY.append({"mode": EXTRACTION_MODE, "case": ccase, "segmentation": seg,
                              "label": label_name, "file": output_file,
                              "merge": "skipped: no existing file"})
              continue
          old = pd.read_excel(output_file, sheet_name=None)
          lacking = [g for g in groups if FULL_N[g] > 0 and g not in old]
          if lacking:
              raise SystemExit(f"{output_file} has sheets {list(old)} but lacks {lacking}: "
                               f"it looks like output of the pre-fix script. Rerun the "
                               f"full script instead of merging.")
          merged = {}
          for g in list(groups) + [k for k in old if k not in groups]:
              o, n = old.get(g), group_frames.get(g)
              if o is not None and n is not None and list(o.columns) != list(n.columns):
                  raise SystemExit(f"{output_file} [{g}]: columns differ from this run "
                                   f"({list(o.columns)} vs {list(n.columns)})")
              parts = []
              if o is not None:
                  parts.append(o[~o["subject"].astype(str).isin(SUBSET)])
              if n is not None:
                  parts.append(n)
              if parts:
                  merged[g] = pd.concat(parts, ignore_index=True).sort_values(
                      "subject", kind="stable").reset_index(drop=True)
          shutil.copy2(output_file, f"{output_file}.bak_{_STAMP}")
          with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
              for g, frame in merged.items():
                  frame.to_excel(writer, index=False, sheet_name=g)
          counts = {g: (len(merged[g]) if g in merged else 0) for g in groups}
          SUMMARY.append({"mode": EXTRACTION_MODE, "case": ccase, "segmentation": seg,
                          "label": label_name, "file": output_file, "merge": "merged",
                          **{f"n_{g}": counts[g] for g in groups},
                          **{f"expected_{g}": FULL_N[g] for g in groups}})
          short = [g for g in groups if counts[g] != FULL_N[g]]
          if short:
              print("  [!] after merge, rows per group differ from the cohort: "
                    + ", ".join(f"{g} {counts[g]}/{FULL_N[g]}" for g in short))
          continue
      # write once, only groups with rows; record n per group for the summary
      counts = {g: (len(group_frames[g]) if g in group_frames else 0) for g in groups}
      expected = {g: len(groups[g][0]) for g in groups}
      SUMMARY.append({"mode": EXTRACTION_MODE, "case": ccase, "segmentation": seg,
                      "label": label_name, "file": output_file,
                      **{f"n_{g}": counts[g] for g in groups},
                      **{f"expected_{g}": expected[g] for g in groups}})
      if not group_frames:
          print(f"  [!] no rows for any group -> {output_file} NOT written")
          continue
      short = [g for g in groups if counts[g] != expected[g]]
      if short:
          print("  [!] rows per group differ from subjects per group: "
                + ", ".join(f"{g} {counts[g]}/{expected[g]}" for g in short)
                + "  (label absent in some subjects, or missing maps)")
      with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
          for group_name, frame in group_frames.items():
              frame.to_excel(writer, index=False, sheet_name=group_name)

# ---- summary: n per group for every file written (or skipped) -------------
_sum = pd.DataFrame(SUMMARY)
_sum_dir = output_xlsx_dir1_tmpl.format(subtree="").rstrip("/").rsplit("/", 1)[0]
os.makedirs(_sum_dir, exist_ok=True)
_sum_path = os.path.join(_sum_dir, f"creation_summary_{EXTRACTION_MODE}.csv" if not SUBSET
                         else f"creation_summary_{EXTRACTION_MODE}_merge-{'_'.join(SUBSET)}.csv")
_sum.to_csv(_sum_path, index=False)
if len(_sum):
    _any_empty = _sum[[c for c in _sum.columns if c.startswith("n_")]].eq(0).any(axis=1)
    print(f"\n[summary] {len(_sum)} files -> {_sum_path}")
    print(f"[summary] files with an EMPTY group: {int(_any_empty.sum())}")
print("Processing complete.")
