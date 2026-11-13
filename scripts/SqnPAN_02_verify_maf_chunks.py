# SqnPAN_02_verify_maf_chunks.py
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# Same checks as SqnPAN_01_verify_maf_download.py (file present, size, md5, header, data row count; shape only, no values),
# but split into chunks so a stopped run can continue. Usage:
#   python SqnPAN_02_verify_maf_chunks.py chunk <raw_dir> <records_tsv> <start> <end>
#   python SqnPAN_02_verify_maf_chunks.py summary <raw_dir> <records_tsv> <log_path>
import sys, os, gzip, hashlib, csv, collections

def chunk(raw_dir, rec_path, start, end):
    rows = list(csv.DictReader(open(os.path.join(raw_dir, "TCGA_masked_maf_manifest_2026-10-21.txt")), delimiter="\t"))
    new = not os.path.exists(rec_path)
    out = open(rec_path, "a")
    try:
        if new: out.write("id\tpresent\tsize_ok\tmd5_ok\tn_cols\tn_rows\theader_md5\n")
        for r in rows[start:end]:
            path = os.path.join(raw_dir, "maf_tcga", r["id"], r["filename"])
            if not os.path.isfile(path):
                out.write(f"{r['id']}\t0\t0\t0\t0\t0\t-\n"); continue
            size_ok = int(os.path.getsize(path) == int(r["size"]))
            h = hashlib.md5()
            with open(path, "rb") as f:
                for blk in iter(lambda: f.read(1 << 20), b""):
                    h.update(blk)
            md5_ok = int(h.hexdigest() == r["md5"])
            hdr, n = None, 0
            with gzip.open(path, "rt") as f:
                for line in f:
                    if line.startswith("#"): continue
                    if hdr is None: hdr = line.rstrip("\n"); continue
                    n += 1
            hh = hashlib.md5(hdr.encode()).hexdigest()[:8]
            out.write(f"{r['id']}\t1\t{size_ok}\t{md5_ok}\t{len(hdr.split(chr(9)))}\t{n}\t{hh}\n")
            out.flush()
        print("chunk", start, end, "done")
    finally:
        out.close()

def summary(raw_dir, rec_path, log_path):
    log = open(log_path, "a")
    def p(*a):
        s = " ".join(str(x) for x in a); print(s); log.write(s + "\n")
    try:
        recs = list(csv.DictReader(open(rec_path), delimiter="\t"))
        ids = [r["id"] for r in recs]
        p("records:", len(recs), " unique ids:", len(set(ids)))
        p("missing files:", sum(r["present"] == "0" for r in recs), [r["id"] for r in recs if r["present"] == "0"])
        p("size mismatches:", sum(r["present"] == "1" and r["size_ok"] == "0" for r in recs))
        p("md5 mismatches:", sum(r["present"] == "1" and r["md5_ok"] == "0" for r in recs))
        pr = [r for r in recs if r["present"] == "1"]
        c = collections.Counter((r["header_md5"], r["n_cols"]) for r in pr)
        p("distinct headers (header md5, columns):", dict(c))
        n = sorted(int(r["n_rows"]) for r in pr)
        p("files with 0 data rows:", sum(x == 0 for x in n))
        p("data rows total:", sum(n), " min/median/max per file:", n[0], n[len(n)//2], n[-1])
    except Exception:
        import traceback; p("ERROR:", traceback.format_exc())
    finally:
        log.close()

if __name__ == "__main__":
    if sys.argv[1] == "chunk": chunk(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))
    else: summary(sys.argv[2], sys.argv[3], sys.argv[4])
# ---- END OF SCRIPT ----
