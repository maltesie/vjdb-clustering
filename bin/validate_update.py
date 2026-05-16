#!/usr/bin/env python3
"""Compare clustering state before and after an update.
Reports any changes in cluster assignments for pre-existing sequences."""

import argparse
import sys
import pandas as pd
from collections import Counter

parser = argparse.ArgumentParser(description="Validate clustering update consistency.")
parser.add_argument("--before", required=True, help="CSV before update (e.g. vjdb1_merged_reps.csv)")
parser.add_argument("--after", required=True, help="CSV after update (e.g. vjdb1_test_merged_reps.csv)")
parser.add_argument("--new", default=None, help="Optional: CSV of just the new sequences")
args = parser.parse_args()

before_df = pd.read_csv(args.before, dtype=str)
after_df = pd.read_csv(args.after, dtype=str)

print(f"Before: {len(before_df)} sequences, {before_df['cluster_rep'].nunique()} clusters")
print(f"After:  {len(after_df)} sequences, {after_df['cluster_rep'].nunique()} clusters")
print()

# --- Check existing sequences ---
before_map = dict(zip(before_df["vjdb_id"], before_df["cluster_rep"]))
after_map = dict(zip(after_df["vjdb_id"], after_df["cluster_rep"]))

existing_ids = set(before_map.keys())
missing = existing_ids - set(after_map.keys())
if missing:
    print(f"ERROR: {len(missing)} sequences from before are missing after update")
    for sid in list(missing)[:5]:
        print(f"  {sid}")

# Check for changed cluster assignments
changed = []
for sid in existing_ids & set(after_map.keys()):
    old_rep = before_map[sid]
    new_rep = after_map[sid]
    if old_rep != new_rep:
        changed.append((sid, old_rep, new_rep))

print(f"Changed cluster assignments: {len(changed)} / {len(existing_ids)} existing sequences")

if changed:
    # Analyze the changes: are representatives getting swapped or merged?
    old_to_new_rep = Counter()
    for sid, old_rep, new_rep in changed:
        old_to_new_rep[(old_rep, new_rep)] += 1

    print(f"  Distinct (old_rep -> new_rep) transitions: {len(old_to_new_rep)}")
    print()

    # Did any old representatives get replaced?
    old_reps_before = set(before_df["cluster_rep"])
    old_reps_after = set(after_df[after_df["vjdb_id"].isin(existing_ids)]["cluster_rep"])
    replaced_reps = old_reps_before - old_reps_after
    if replaced_reps:
        print(f"  Representatives fully replaced: {len(replaced_reps)}")
        for rep in list(replaced_reps)[:10]:
            # Find what they were replaced by
            successors = set()
            for sid, old_rep, new_rep in changed:
                if old_rep == rep:
                    successors.add(new_rep)
            old_size = sum(1 for v in before_map.values() if v == rep)
            print(f"    {rep} (was {old_size} members) -> {successors}")

    # Did any clusters get merged?
    merge_targets = Counter(new_rep for _, _, new_rep in changed)
    merges = {rep: count for rep, count in merge_targets.items() if count > 1}
    if merges:
        print(f"\n  Clusters absorbing members from multiple old clusters:")
        for rep, count in sorted(merges.items(), key=lambda x: -x[1])[:10]:
            sources = set(old_rep for _, old_rep, new_rep in changed if new_rep == rep)
            print(f"    {rep} <- {len(sources)} old clusters ({count} reassigned seqs)")

    # Did any clusters get split?
    split_sources = Counter(old_rep for _, old_rep, _ in changed)
    splits = {}
    for old_rep, count in split_sources.items():
        destinations = set(new_rep for _, o, new_rep in changed if o == old_rep)
        # Also count members that stayed
        stayed = sum(1 for sid in existing_ids if before_map[sid] == old_rep and after_map.get(sid) == old_rep)
        if stayed > 0:
            destinations.add(old_rep)
        if len(destinations) > 1:
            splits[old_rep] = destinations
    if splits:
        print(f"\n  Clusters that were split across multiple new clusters:")
        for old_rep, dests in sorted(splits.items(), key=lambda x: -len(x[1]))[:10]:
            print(f"    {old_rep} -> {len(dests)} clusters: {dests}")

    print()
    print("  First 10 individual changes:")
    for sid, old_rep, new_rep in changed[:10]:
        print(f"    {sid}: {old_rep} -> {new_rep}")

# --- Check new sequences ---
if args.new:
    new_df = pd.read_csv(args.new, dtype=str)
    new_ids = set(new_df["seq_id"] if "seq_id" in new_df.columns else new_df["vjdb_id"])
    new_assigned = new_df[new_df["cluster_rep"].isin(old_reps_before if 'old_reps_before' in dir() else set())]
    print(f"\nNew sequences: {len(new_df)} total")
    print(f"  Assigned to pre-existing clusters: {len(new_assigned)}")
    print(f"  Forming new clusters: {len(new_df) - len(new_assigned)}")

# --- Summary ---
print("\n--- Summary ---")
issues = 0
if missing:
    print(f"FAIL: {len(missing)} sequences lost")
    issues += 1
if changed:
    print(f"WARN: {len(changed)} sequences changed cluster assignment")
    issues += 1
if not missing and not changed:
    print("OK: All existing sequences retained their cluster assignments")

sys.exit(1 if missing else 0)
