# SqnPAN_01_verify_maf_download.py
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# Checks the GDC download against the manifest: file present, size, md5.
# Also checks that every file has the same header and counts data rows (shape only, no values read).
import sys, os, gzip, hashlib, csv, collections

def main(raw_dir, log_path):
    manifest = os.path.join(raw_dir, "TCGA_masked_maf_manifest_2026-10-21.txt")
    maf_dir = os.path.join(raw_dir, "maf_tcga")
    log = open(log_path, "a")
    def p(*a):
        s = " ".join(str(x) for x in a)
        print(s); log.write(s + "\n"); log.flush()
    try:
        rows = list(csv.DictReader(open(manifest), delimiter="\t"))
        p("manifest rows:", len(rows))
        missing, bad_size, bad_md5, headers, nrows = [], [], [], collections.Counter(), []
        for i, r in enumerate(rows, 1):
            path = os.path.join(maf_dir, r["id"], r["filename"])
            if not os.path.isfile(path):
                missing.append(r["id"]); continue
            if os.path.getsize(path) != int(r["size"]):
                bad_size.append(r["id"])
            h = hashlib.md5()
            with open(path, "rb") as f:
                for blk in iter(lambda: f.read(1 << 20), b""):
                    h.update(blk)
            if h.hexdigest() != r["md5"]:
                bad_md5.append(r["id"])
            with gzip.open(path, "rt") as f:
                n = 0; hdr = None
                for line in f:
                    if line.startswith("#"): continue
                    if hdr is None: hdr = line.rstrip("\n"); continue
                    n += 1
            headers[hdr] += 1; nrows.append(n)
            if i % 1000 == 0: p("progress:", i, "of", len(rows))
        p("missing files:", len(missing))
        p("size mismatches:", len(bad_size))
        p("md5 mismatches:", len(bad_md5))
        if missing[:5]: p("first missing ids:", missing[:5])
        if bad_md5[:5]: p("first md5 mismatch ids:", bad_md5[:5])
        p("distinct headers:", len(headers))
        for h, c in headers.items():
            p("  header with", len(h.split("\t")), "columns in", c, "files")
        ids = {r["id"] for r in rows}
        extra = [d for d in os.listdir(maf_dir) if d not in ids]
        p("folders not in manifest:", len(extra), extra[:5])
        nrows.sort()
        p("files with 0 data rows:", sum(1 for n in nrows if n == 0))
        p("data rows total:", sum(nrows), " min/median/max per file:", nrows[0], nrows[len(nrows)//2], nrows[-1])
    except Exception as e:
        import traceback; p("ERROR:", traceback.format_exc())
    finally:
        log.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
# ---- END OF SCRIPT ----
