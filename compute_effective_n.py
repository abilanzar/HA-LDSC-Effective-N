#!/usr/bin/env python3
"""Replace LDSC participant N with summed cohort-specific effective N.

Inputs contain summary-level SNP identifiers and cohort membership only. Never
put participant-level records in these files or in a public repository.
"""

import argparse
import csv
import gzip
import io
from pathlib import Path


def open_text(path):
    return gzip.open(path, "rt", newline="") if str(path).endswith(".gz") else open(path, "r", newline="")


def load_counts(path):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["cohort", "cases", "controls"]:
            raise ValueError("Counts header must be: cohort, cases, controls")
        counts = {}
        for row in reader:
            cohort = row["cohort"]
            if not cohort or cohort in counts:
                raise ValueError("Cohort names must be nonempty and unique")
            cases, controls = int(row["cases"]), int(row["controls"])
            if cases <= 0 or controls <= 0:
                raise ValueError(f"Positive case/control counts required for {cohort}")
            counts[cohort] = (cases, controls)
    if not counts:
        raise ValueError("No cohort counts supplied")
    return counts


def effective_n(cases, controls):
    return 4.0 * cases * controls / (cases + controls)


def convert(sumstats, membership, counts_file, output):
    counts = load_counts(counts_file)
    output = Path(output)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite: {output}")
    if output.resolve() in (Path(sumstats).resolve(), Path(membership).resolve(), Path(counts_file).resolve()):
        raise ValueError("Output must differ from every input")

    # Validate every row before creating an output; fail closed on mismatches.
    import tempfile
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="ldsc_neff_", dir=output.parent) as tempdir:
        staged = Path(tempdir) / "validated.tsv.gz"
        with open_text(sumstats) as source, open_text(membership) as members:
            src = csv.DictReader(source, delimiter="\t")
            mem = csv.DictReader(members, delimiter="\t")
            if src.fieldnames != ["SNP", "A1", "A2", "N", "Z"]:
                raise ValueError("LDSC header must be: SNP, A1, A2, N, Z")
            if mem.fieldnames != ["SNP", *counts]:
                raise ValueError("Membership header/order must be SNP followed by count-file cohort order")
            with staged.open("wb") as raw:
                with gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as compressed:
                    with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                        writer = csv.DictWriter(text, fieldnames=src.fieldnames, delimiter="\t", lineterminator="\n")
                        writer.writeheader()
                        rows = 0
                        for original, flags in zip(src, mem):
                            rows += 1
                            if not original["SNP"] or original["SNP"] != flags["SNP"]:
                                raise ValueError(f"SNP/order mismatch at row {rows}")
                            if not all(original[col] for col in ("A1", "A2", "N", "Z")):
                                raise ValueError(f"Missing LDSC field at row {rows}")
                            float(original["N"])
                            float(original["Z"])
                            contributors = []
                            for cohort in counts:
                                flag = flags[cohort]
                                if flag not in ("0", "1"):
                                    raise ValueError(f"Membership must be 0/1 at row {rows}, cohort {cohort}")
                                if flag == "1":
                                    contributors.append(cohort)
                            if not contributors:
                                raise ValueError(f"No contributing cohort at row {rows}")
                            neff = sum(effective_n(*counts[name]) for name in contributors)
                            updated = dict(original)
                            updated["N"] = format(neff, ".12g")
                            writer.writerow(updated)
                        if next(src, None) is not None or next(mem, None) is not None:
                            raise ValueError("Input row counts differ")
        if not rows:
            raise ValueError("No SNP rows supplied")
        staged.replace(output)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sumstats", required=True, help="LDSC SNP/A1/A2/N/Z .tsv.gz")
    parser.add_argument("--membership", required=True, help="Aligned SNP-by-cohort 0/1 .tsv.gz")
    parser.add_argument("--counts", required=True, help="Cohort case/control counts .tsv")
    parser.add_argument("--out", required=True, help="New .sumstats.gz output; never overwritten")
    args = parser.parse_args()
    print(f"Validated and wrote {convert(args.sumstats, args.membership, args.counts, args.out)} SNPs")


if __name__ == "__main__":
    main()
