# SqnPAN_03_check_maf_uniprot_mapping.py
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# Consistency check between the MAF protein columns and the UniProt 2026_03 reviewed sequences (the sequences AlphaFold models).
# Counts only identifier and coordinate agreement for missense rows; no mutation result of interest is computed.
# Usage: chunk <raw_dir> <records_jsonl> <start> <end>   |   summary <records_jsonl> <log_path>
import sys, os, gzip, csv, json, collections

def load_uniprot(raw_dir):
    seq = {}
    for row in csv.DictReader(gzip.open(os.path.join(raw_dir, "Uniprot", "stream.gz"), "rt"), delimiter="\t"):
        seq[row["Entry"]] = row["Sequence"]
    return seq

def chunk(raw_dir, rec_path, start, end):
    seq = load_uniprot(raw_dir)
    ids = [r["id"] for r in csv.DictReader(open(os.path.join(raw_dir, "TCGA_masked_maf_manifest_2026-10-21.txt")), delimiter="\t")]
    c = collections.Counter()
    for fid in ids[start:end]:
        d = os.path.join(raw_dir, "maf_tcga", fid)
        path = os.path.join(d, [x for x in os.listdir(d) if x.endswith(".maf.gz")][0])
        hdr = None
        with gzip.open(path, "rt") as f:
            for line in f:
                if line.startswith("#"): continue
                if hdr is None:
                    hdr = line.rstrip("\n").split("\t")
                    ix = {k: hdr.index(k) for k in ("Variant_Classification", "SWISSPROT", "Protein_position", "Amino_acids", "CANONICAL", "Transcript_ID")}
                    continue
                v = line.rstrip("\n").split("\t")
                c["rows"] += 1
                if v[ix["Variant_Classification"]] != "Missense_Mutation": continue
                c["missense"] += 1
                c["missense_canonical_YES" if v[ix["CANONICAL"]] == "YES" else "missense_canonical_not_YES"] += 1
                acc = v[ix["SWISSPROT"]].split(".")[0].split("&")[0]
                if not acc: c["no_swissprot"] += 1; continue
                if acc not in seq: c["acc_not_in_uniprot_reviewed"] += 1; continue
                c["acc_in_uniprot"] += 1
                pp = v[ix["Protein_position"]]; aa = v[ix["Amino_acids"]]
                try:
                    pos_s, tot_s = pp.split("/"); pos = int(pos_s); tot = int(tot_s); ref = aa.split("/")[0]
                except Exception:
                    c["unparsable_position_or_aa"] += 1; continue
                s = seq[acc]
                if tot == len(s): c["length_match"] += 1
                else: c["length_mismatch"] += 1
                if 1 <= pos <= len(s) and len(ref) == 1 and s[pos - 1] == ref: c["ref_aa_match"] += 1
                else: c["ref_aa_mismatch"] += 1
                if tot == len(s) and 1 <= pos <= len(s) and s[pos - 1] == ref: c["both_match"] += 1
    with open(rec_path, "a") as o: o.write(json.dumps({"start": start, "end": end, "counts": c}) + "\n")
    print("chunk", start, end, "done", dict(c))

def summary(rec_path, log_path):
    tot = collections.Counter(); seen = []
    for l in open(rec_path):
        r = json.loads(l); tot.update(r["counts"]); seen.append((r["start"], r["end"]))
    with open(log_path, "a") as log:
        for s in (f"chunks: {sorted(seen)}", *[f"{k}: {v}" for k, v in sorted(tot.items())]):
            print(s); log.write(s + "\n")

if __name__ == "__main__":
    if sys.argv[1] == "chunk": chunk(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))
    else: summary(sys.argv[2], sys.argv[3])
# ---- END OF SCRIPT ----
