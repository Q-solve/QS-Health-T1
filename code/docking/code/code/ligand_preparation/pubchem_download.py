"""
pubchem_download.py

Reads SMILES from an Excel file, looks each one up on PubChem (exact match
first, falling back to a 2D similarity search), downloads structures as
individual SDF files, and merges everything into one combined SDF file.

INCREMENTAL / APPEND BEHAVIOUR:
Each run reads the existing log (_download_log.csv in OUTPUT_DIR, if present)
to see which input SMILES were already successfully resolved. Those are
skipped (not re-downloaded from PubChem) -- so you can add new rows to your
Excel file later and simply rerun the script; only the new compounds get
queried. The combined SDF file is rebuilt each run from every .sdf file
currently in OUTPUT_DIR, so it always reflects everything downloaded so far
(old + new).

COORDINATES:
PubChem SDF records can be 2D or 3D. This script requests 3D coordinates
first (RECORD_TYPE = "3d") and falls back to 2D for any CID where PubChem
doesn't have a precomputed 3D conformer (common for less-studied compounds).
The log records which type each compound actually got -- check the
'coord_type' column if you need every structure to be true 3D (you'd then
need to generate 3D conformers yourself, e.g. with RDKit, for any '2d' rows).

SMILES are POSTed (not put in the URL) since they often contain characters
(slashes, #, +, etc.) that break a URL path.

Requests are throttled to ~20 per minute to stay well within PubChem's
usage limits (max 5 requests/second).
"""

import os
import time
import pandas as pd
import requests

# ---------------------- CONFIG ---------------------------------------
INPUT_XLSX = "smiles_list.xlsx"     # path to your Excel file
SMILES_COLUMN = "Smiles"            # column name containing SMILES strings
ID_COLUMN = "compoundID"            # optional: set to None if you don't have one
SHEET_NAME = 0

OUTPUT_DIR = "pubchem_structures"   # folder where files will be saved
COMBINED_SDF_NAME = "all_compounds_combined.sdf"

SIMILARITY_THRESHOLD = 90           # Tanimoto % threshold for similarity fallback
MAX_SIMILAR_HITS = 1                # how many similar compounds per SMILES if no exact match

REQUESTS_PER_MINUTE = 20
SECONDS_BETWEEN_REQUESTS = 60.0 / REQUESTS_PER_MINUTE

PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
VERBOSE = True
# -----------------------------------------------------------------------


def throttle():
    time.sleep(SECONDS_BETWEEN_REQUESTS)


def safe_filename(s, maxlen=60):
    keep = "".join(c if c.isalnum() or c in "-_." else "_" for c in s)
    return keep[:maxlen] if keep else "compound"


def log_fail(label, resp):
    if VERBOSE:
        print(f"    [{label}] HTTP {resp.status_code}: {resp.text[:300]}")


def exact_match_cids(smiles):
    url = f"{PUG_BASE}/compound/smiles/cids/JSON"
    r = requests.post(url, data={"smiles": smiles}, timeout=30)
    throttle()
    if r.status_code == 200:
        try:
            return r.json().get("IdentifierList", {}).get("CID", [])
        except ValueError:
            log_fail("exact-match: bad JSON", r)
            return []
    log_fail("exact-match", r)
    return []


def similarity_search_cids(smiles, threshold=90, max_hits=1):
    url = f"{PUG_BASE}/compound/fastsimilarity_2d/smiles/cids/JSON"
    params = {"Threshold": threshold, "MaxRecords": max_hits}
    r = requests.post(url, params=params, data={"smiles": smiles}, timeout=30)
    throttle()

    if r.status_code == 202:
        try:
            listkey = r.json()["Waiting"]["ListKey"]
        except (ValueError, KeyError):
            log_fail("similarity: bad ListKey response", r)
            return []
        poll_url = f"{PUG_BASE}/compound/listkey/{listkey}/cids/JSON"
        for _ in range(20):
            time.sleep(2)
            pr = requests.get(poll_url, timeout=30)
            if pr.status_code == 200:
                try:
                    pdata = pr.json()
                except ValueError:
                    log_fail("similarity poll: bad JSON", pr)
                    return []
                if "IdentifierList" in pdata:
                    throttle()
                    return pdata["IdentifierList"]["CID"][:max_hits]
                if "Fault" in pdata:
                    log_fail("similarity poll: fault", pr)
                    return []
            elif pr.status_code != 202:
                log_fail("similarity poll", pr)
                break
        throttle()
        return []
    elif r.status_code == 200:
        try:
            return r.json().get("IdentifierList", {}).get("CID", [])[:max_hits]
        except ValueError:
            log_fail("similarity: bad JSON", r)
            return []
    else:
        log_fail("similarity", r)
        return []


def download_sdf(cid, out_path):
    """Try 3D first, fall back to 2D. Returns coord_type used, or None on failure."""
    for record_type in ("3d", "2d"):
        url = f"{PUG_BASE}/compound/cid/{cid}/SDF"
        params = {"record_type": record_type}
        r = requests.get(url, params=params, timeout=30)
        throttle()
        if r.status_code == 200 and r.content:
            with open(out_path, "wb") as f:
                f.write(r.content)
            return record_type
        else:
            if VERBOSE and record_type == "3d":
                print(f"    No 3D conformer for CID {cid}, falling back to 2D...")
            else:
                log_fail(f"download SDF CID {cid}", r)
    return None


def load_existing_log(log_path):
    """Return set of input_smiles that were already successfully resolved."""
    if os.path.exists(log_path):
        try:
            prev = pd.read_csv(log_path)
            done = set(prev.loc[prev["match_type"] != "none", "input_smiles"].astype(str))
            return prev.to_dict("records"), done
        except Exception:
            return [], set()
    return [], set()


def merge_all_sdfs(output_dir, combined_name):
    """Concatenate every individual .sdf file in output_dir into one combined file."""
    combined_path = os.path.join(output_dir, combined_name)
    sdf_files = sorted(
        f for f in os.listdir(output_dir)
        if f.endswith(".sdf") and f != combined_name
    )
    count = 0
    with open(combined_path, "wb") as out:
        for fname in sdf_files:
            with open(os.path.join(output_dir, fname), "rb") as f:
                content = f.read()
                out.write(content)
                if not content.rstrip().endswith(b"$$$$"):
                    out.write(b"\n$$$$\n")
                count += 1
    print(f"\nCombined SDF written: {combined_path} ({count} structures)")
    return combined_path, count


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    log_path = os.path.join(OUTPUT_DIR, "_download_log.csv")

    df = pd.read_excel(INPUT_XLSX, sheet_name=SHEET_NAME)
    if SMILES_COLUMN not in df.columns:
        raise ValueError(
            f"Column '{SMILES_COLUMN}' not found. Available columns: {list(df.columns)}"
        )

    prev_rows, already_done = load_existing_log(log_path)
    print(f"Found {len(already_done)} previously-resolved compound(s) in existing log; will skip those.")

    smiles_list = df[SMILES_COLUMN].dropna().astype(str).tolist()
    ids_list = (
        df[ID_COLUMN].astype(str).tolist() if (ID_COLUMN and ID_COLUMN in df.columns)
        else [None] * len(smiles_list)
    )
    print(f"Loaded {len(smiles_list)} SMILES from {INPUT_XLSX}")

    new_rows = []
    counts = {"exact": 0, "similar": 0, "none": 0, "skipped_cached": 0}

    for i, (smi, cmp_id) in enumerate(zip(smiles_list, ids_list), 1):
        smi = smi.strip()

        if smi in already_done:
            print(f"[{i}/{len(smiles_list)}] {smi} -- already downloaded, skipping.")
            counts["skipped_cached"] += 1
            continue

        print(f"\n[{i}/{len(smiles_list)}] {smi}")
        match_type = None
        cids = []

        try:
            cids = exact_match_cids(smi)
            if cids:
                match_type = "exact"
            else:
                cids = similarity_search_cids(smi, threshold=SIMILARITY_THRESHOLD, max_hits=MAX_SIMILAR_HITS)
                if cids:
                    match_type = f"similar(>{SIMILARITY_THRESHOLD}%)"
        except requests.RequestException as e:
            print(f"  Request error: {e}")

        if not cids:
            print("  No match found.")
            counts["none"] += 1
            new_rows.append({
                "compound_id": cmp_id, "input_smiles": smi, "match_type": "none",
                "cid": None, "file": None, "coord_type": None,
            })
            continue

        base_name = safe_filename(f"{cmp_id}_{smi}" if cmp_id else smi)

        for cid in cids:
            fname_root = f"{base_name}_CID{cid}"
            sdf_path = os.path.join(OUTPUT_DIR, f"{fname_root}.sdf")
            coord_type = download_sdf(cid, sdf_path)

            print(f"  CID {cid} ({match_type}): SDF={'ok (' + coord_type + ')' if coord_type else 'FAIL'}")

            counts["exact" if match_type == "exact" else "similar"] += 1 if coord_type else 0
            new_rows.append({
                "compound_id": cmp_id, "input_smiles": smi, "match_type": match_type,
                "cid": cid, "file": fname_root, "coord_type": coord_type,
            })

    all_rows = prev_rows + new_rows
    log_df = pd.DataFrame(all_rows)
    log_df.to_csv(log_path, index=False)

    merge_all_sdfs(OUTPUT_DIR, COMBINED_SDF_NAME)

    print("\n--- Summary (this run) ---")
    print(f"  Exact matches:      {counts['exact']}")
    print(f"  Similarity matches: {counts['similar']}")
    print(f"  No match found:     {counts['none']}")
    print(f"  Skipped (cached):   {counts['skipped_cached']}")
    print(f"  Full log:           {log_path}")


if __name__ == "__main__":
    main()
