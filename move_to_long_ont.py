#!/usr/bin/env python3
"""
Moves everything except boundary_nodes.txt and info.txt from each
TTT_runs/<RUN_ID>/ into TTT_runs/<RUN_ID>/long_ont/, then fixes
paths in the moved run.sh scripts.
"""

import os
import shutil
import re

BASE_DIR = "/data/antipovd2/res/TTT_paper/giraffe/TTT_runs"
KEEP_FILES = {"boundary_nodes.txt", "info.txt"}

def process_run_dir(run_dir):
    """Move files from run_dir into run_dir/long_ont/, keeping boundaries and info.txt."""
    long_ont_dir = os.path.join(run_dir, "long_ont")

    entries = os.listdir(run_dir)
    to_move = [e for e in entries if e not in KEEP_FILES and e != "long_ont"]

    if not to_move:
        print(f"  Nothing to move in {run_dir}")
        return

    os.makedirs(long_ont_dir, exist_ok=True)

    for entry in to_move:
        src = os.path.join(run_dir, entry)
        dst = os.path.join(long_ont_dir, entry)
        if os.path.exists(dst):
            print(f"  WARNING: {dst} already exists, skipping {entry}")
            continue
        shutil.move(src, dst)
        print(f"  Moved: {entry}")

    # Fix paths in run.sh if it was moved
    run_sh = os.path.join(long_ont_dir, "run.sh")
    if os.path.isfile(run_sh):
        fix_run_sh(run_sh)


def fix_run_sh(run_sh_path):
    """Fix relative paths in run.sh now that it lives in long_ont/ subdirectory."""
    with open(run_sh_path, "r") as f:
        content = f.read()

    original = content

    # boundary_nodes.txt is now in parent directory
    content = content.replace(
        "--boundary-nodes boundary_nodes.txt",
        "--boundary-nodes ../boundary_nodes.txt"
    )

    # --outdir . should point to current dir (long_ont/) which is fine,
    # but if user wants results in long_ont/ we keep it as .
    # No change needed for --outdir .

    if content != original:
        with open(run_sh_path, "w") as f:
            f.write(content)
        print(f"  Fixed paths in run.sh")


def main():
    for entry in sorted(os.listdir(BASE_DIR)):
        run_dir = os.path.join(BASE_DIR, entry)
        if not os.path.isdir(run_dir):
            continue
        # Skip non-run directories (e.g. giraffe_detection)
        if not entry.startswith("tangle_") and not entry.startswith("manual_tangle_"):
            continue
        print(f"Processing {entry}:")
        process_run_dir(run_dir)


if __name__ == "__main__":
    main()
