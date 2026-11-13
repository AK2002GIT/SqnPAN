# SqnPAN_10_foldx_helpers.py
# Part of the SqnPAN analysis code (A. Balasundar, 2026; MIT licence). Run order, inputs and outputs: README.md
# Reusable functions for running FoldX 5.1 on AlphaFold v6 models. SqnPAN_10_foldx_setup.py tests them; the FoldX run script (SqnPAN_11) imports them:
# reuse them, do not rewrite them.
# - model_pdb_text: the PDB text of one AlphaFold model (fragment) read from the packed tar through the SqnPAN_04 index (rule 3A; the tar is never unpacked).
# - pdb_ca_residues / check_model_sequence: residue list of a PDB (chain, number, residue) and the check against the UniProt 2026_03 sequence.
# - run_foldx: runs the FoldX program once in a working folder, stdout to a file, measures seconds and (if psutil is installed) peak memory.
# - repair / build: RepairPDB and BuildModel calls. build writes the mutation list itself and refuses a mutation whose wild-type residue
#   is not the residue in the PDB at that number (our own guard; FoldX's own behaviour on a wrong residue is tested in SqnPAN_10 section 1).
# - read_fx_table / dif_ddg / average_ddg: parsers for the FoldX output tables.
# FROM MEMORY (marked; tested in SqnPAN_10 section 1 before any run): the mutation-list argument name (--mutant-file), the list file name
#   (individual_list.txt), the line format of the list (wild type + chain + number + mutant + ";", for example KA42R;), the output file names
#   (<base>_Repair.pdb; Dif_<base>.fxout, Average_<base>.fxout, mutant PDBs <base>_<mutation>_<run>.pdb) and the name pattern of the rows of the
#   Dif table. Argument names and defaults of the physics settings are from the program's own help (FoldX help output, analysis plan Amendment 3).
import os, sys, re, time, shutil, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import SqnPAN_05_helpers as H
import SqnPAN_04_alphafold_helpers as AF

FOLDX_EXE = os.path.join(H.PAPER2, "foldx5_Windows", "foldx_1_20270131.exe")
FOLDX_MD5 = "c66ffd383d4d1d274e940bf99979bf00"            # recorded at installation
FOLDX_EXPIRY = "2026-12-31"
FOLDX_WORK = os.path.join(H.PROC_OUT, "P2_10_foldx")      # working folders, extracted models and repaired models
REPAIRED_DIR = os.path.join(FOLDX_WORK, "repaired")       # RepairPDB outputs kept for reuse by P2_11 (rule 3A: only the models FoldX needs)
# Physics settings = the FoldX 5.1 defaults printed by the program's help (analysis plan Amendment 3, decision D3 A). They are passed explicitly;
# SqnPAN_10 section 1 checks that the result is identical when they are left out (so the explicit values equal the program defaults).
PHYS = [("temperature", "298"), ("pH", "7"), ("ionStrength", "0.05"), ("water", "-IGNORE"), ("vdwDesign", "2"), ("moveNeighbours", "true")]
N_RUNS = 3                                                # brief section 8
MUTANT_ARG = "mutant-file"                                # FROM MEMORY, tested in P2_10 section 1
MUTANT_FILE = "individual_list.txt"                       # FROM MEMORY, tested in P2_10 section 1
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
         "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
ONE = {v: k for k, v in THREE.items()}

def phys_args():
    return [f"--{k}={v}" for k, v in PHYS]

def model_pdb_text(fh, idx, acc, fragment):
    return AF.read_model(fh, idx, acc, fragment, "pdb")

def pdb_ca_residues(text):
    """[(chain, residue number, insertion code, three-letter name)] for the CA atoms of ATOM records, in file order."""
    out = []
    for l in text.split("\n"):
        if l.startswith("ATOM") and l[12:16].strip() == "CA":
            out.append((l[21], int(l[22:26]), l[26].strip(), l[17:20].strip()))
    return out

def record_types(text):
    c = {}
    for l in text.split("\n"):
        k = l[:6].strip()
        if k: c[k] = c.get(k, 0) + 1
    return c

def check_model_sequence(text, seq, align_begin):
    """Compares the CA residues of a model with the UniProt sequence (model residue i = UniProt position align_begin + i - 1).
    Returns (n_residues, n_mismatch, chains, numbering_ok)."""
    res = pdb_ca_residues(text); mism = 0
    for k, (ch, num, ic, name) in enumerate(res):
        p = align_begin + num - 1
        if not (1 <= p <= len(seq)) or THREE.get(name) != seq[p - 1]: mism += 1
    numbering_ok = [r[1] for r in res] == list(range(1, len(res) + 1)) and all(r[2] == "" for r in res)
    return len(res), mism, sorted({r[0] for r in res}), numbering_ok

def run_foldx(args, cwd, tag, poll=0.5):
    """One FoldX call in folder cwd. stdout and stderr go to cwd/stdout_<tag>.txt. Returns rc, seconds, peak_mb (None without psutil), stdout path."""
    try:
        import psutil
    except ImportError:
        psutil = None
    out = os.path.join(cwd, f"stdout_{tag}.txt"); peak = None; t0 = time.time()
    with open(out, "w", encoding="utf-8", errors="replace") as fo:
        p = subprocess.Popen([FOLDX_EXE] + list(args), cwd=cwd, stdout=fo, stderr=subprocess.STDOUT)
        pp = None
        if psutil is not None:
            try: pp = psutil.Process(p.pid)
            except Exception: pp = None
        while p.poll() is None:
            if pp is not None:
                try: peak = max(peak or 0, pp.memory_info().rss)
                except Exception: pass
            time.sleep(poll)
    return {"rc": p.returncode, "seconds": time.time() - t0, "peak_mb": None if peak is None else round(peak / 2 ** 20, 1), "stdout": out,
            "args": " ".join(args)}

def stdout_tail(path, n=15):
    try:
        with open(path, encoding="utf-8", errors="replace") as f: lines = [l.rstrip() for l in f if l.strip()]
    except OSError:
        return ["(no stdout file)"]
    return lines[-n:]

def _place(src, workdir):
    os.makedirs(workdir, exist_ok=True)
    name = os.path.basename(src); dst = os.path.join(workdir, name)
    if os.path.abspath(src) != os.path.abspath(dst): shutil.copy2(src, dst)
    return name

def repair(src_pdb, workdir, phys=True, tag=None):
    """RepairPDB on src_pdb inside workdir. Adds 'repaired' (expected output path, FROM MEMORY name <base>_Repair.pdb) and 'exists'."""
    name = _place(src_pdb, workdir); base = name[:-4]
    r = run_foldx(["--command=RepairPDB", f"--pdb={name}"] + (phys_args() if phys else []), workdir, tag or f"RepairPDB_{base}")
    r["repaired"] = os.path.join(workdir, base + "_Repair.pdb"); r["exists"] = os.path.exists(r["repaired"])
    return r

def mutation_line(wt, chain, num, mt):
    return f"{wt}{chain}{num}{mt};"

def mutant_pdbs(workdir, base):
    """Mutant and wild-type PDBs written by BuildModel: <base>_<mutation>_<run>.pdb and WT_<base>_<mutation>_<run>.pdb."""
    pat = re.compile(r"(WT_)?" + re.escape(base) + r"_\d+_\d+\.pdb")
    return sorted(os.path.join(workdir, f) for f in os.listdir(workdir) if pat.fullmatch(f))

def build(pdb_path, muts, workdir, phys=True, out_pdb=None, runs=N_RUNS, tag="BuildModel", check_wt=True, keep_pdbs=False):
    """BuildModel on pdb_path (a repaired model) with the mutations muts = [(wt, chain, model number, mutant)], one mutant per line.
    out_pdb None (amendment of 2026-10-25, option B): no --out-pdb argument, so FoldX runs at its default (mutant PDBs written); True/False pass 1/0.
    keep_pdbs False: the mutant and WT PDBs are deleted after the call (the tables are read from the .fxout files); 'pdbs_written' counts them.
    check_wt: refuse (ValueError) if a wild-type letter differs from the PDB residue at that number (set False only for the FoldX behaviour test)."""
    name = _place(pdb_path, workdir); base = name[:-4]
    if check_wt:
        with open(os.path.join(workdir, name), encoding="utf-8", errors="replace") as f: res = pdb_ca_residues(f.read())
        have = {(ch, num): THREE.get(nm) for ch, num, ic, nm in res}
        for wt, ch, num, mt in muts:
            if have.get((ch, num)) != wt: raise ValueError(f"{name}: wild type {wt} at {ch}{num} but the PDB has {have.get((ch, num))}")
            if mt not in ONE or mt == wt: raise ValueError(f"bad mutant residue {mt} at {ch}{num}")
    with open(os.path.join(workdir, MUTANT_FILE), "w", newline="\n") as f:
        for m in muts: f.write(mutation_line(*m) + "\n")
    args = ["--command=BuildModel", f"--pdb={name}", f"--{MUTANT_ARG}={MUTANT_FILE}", f"--numberOfRuns={runs}"] + \
           ([] if out_pdb is None else [f"--out-pdb={1 if out_pdb else 0}"]) + (phys_args() if phys else [])
    r = run_foldx(args, workdir, tag)
    r["base"] = base; r["dif"] = os.path.join(workdir, f"Dif_{base}.fxout"); r["average"] = os.path.join(workdir, f"Average_{base}.fxout")
    pdbs = mutant_pdbs(workdir, base); r["pdbs_written"] = len(pdbs)
    if not keep_pdbs:
        for p in pdbs: os.remove(p)
    return r

def read_fx_table(path):
    """FoldX .fxout table: the header is the line whose first field is 'Pdb'; returns (header, rows as lists of fields)."""
    header = None; rows = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = [p.strip() for p in line.rstrip("\r\n").split("\t")]
            if header is None:
                if parts and parts[0] == "Pdb": header = parts
                continue
            if line.strip(): rows.append(parts)
    if header is None: raise ValueError(f"no 'Pdb' header line in {path}")
    return header, rows

DIF_RE = re.compile(r"_(\d+)_(\d+)\.pdb$")     # FROM MEMORY: <base>_<mutation number>_<run number>.pdb
AVG_RE = re.compile(r"_(\d+)(?:\.pdb)?$")      # <base>_<mutation number>, no .pdb (seen in the FoldX 5.1 output, P2_10 section 1; .pdb also accepted)

def dif_ddg(path):
    """{mutation number (1-based line of the list): {run number: ddG total energy as text}}. Values stay text so that two runs compare exactly."""
    header, rows = read_fx_table(path); ci = header.index("total energy"); out = {}
    for p in rows:
        m = DIF_RE.search(os.path.basename(p[0]))
        if m is None: raise ValueError(f"unexpected row name {p[0]} in {path}")
        out.setdefault(int(m.group(1)), {})[int(m.group(2))] = p[ci]
    return out

def average_ddg(path):
    """{mutation number: (SD text, mean total energy text)} from the Average table."""
    header, rows = read_fx_table(path); ci = header.index("total energy"); si = header.index("SD"); out = {}
    for p in rows:
        m = AVG_RE.search(os.path.basename(p[0]))
        if m is None: raise ValueError(f"unexpected row name {p[0]} in {path}")
        out[int(m.group(1))] = (p[si], p[ci])
    return out

def file_head(path, n=12):
    out = []
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            for l in f:
                if len(out) >= n: break
                out.append(l.rstrip("\r\n"))
    except OSError:
        return ["(file not found)"]
    return out
# ---- END OF HELPERS ----
