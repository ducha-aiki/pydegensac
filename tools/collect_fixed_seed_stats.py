#!/usr/bin/env python
import argparse
import csv
import json
import statistics
from pathlib import Path

from run_fixed_seed import run_estimation


def write_summary(rows, output_path):
    h_times = [row["h_time"] for row in rows]
    f_times = [row["f_time"] for row in rows]
    totals = [row["total_time"] for row in rows]
    h_inliers = [row["h_inliers"] for row in rows]
    f_inliers = [row["f_inliers"] for row in rows]

    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        writer.writerow(["runs", len(rows)])
        writer.writerow(["h_time_mean", statistics.mean(h_times)])
        writer.writerow(["h_time_stdev", statistics.pstdev(h_times)])
        writer.writerow(["f_time_mean", statistics.mean(f_times)])
        writer.writerow(["f_time_stdev", statistics.pstdev(f_times)])
        writer.writerow(["total_time_mean", statistics.mean(totals)])
        writer.writerow(["total_time_stdev", statistics.pstdev(totals)])
        writer.writerow(["h_inliers_mean", statistics.mean(h_inliers)])
        writer.writerow(["f_inliers_mean", statistics.mean(f_inliers)])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed-start", type=int, default=1234)
    parser.add_argument("--count", type=int, default=10)
    parser.add_argument("--stats-output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    parser.add_argument("--json-dir", type=Path, default=None)
    args = parser.parse_args()

    args.stats_output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    if args.json_dir is not None:
        args.json_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    with args.stats_output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "run",
                "seed",
                "h_time",
                "f_time",
                "total_time",
                "h_inliers",
                "f_inliers",
                "h_signature",
                "f_signature",
            ]
        )

        for offset in range(args.count):
            seed = args.seed_start + offset
            result = run_estimation(seed)
            row = {
                "run": offset + 1,
                "seed": seed,
                "h_time": result["homography"]["runtime_sec"],
                "f_time": result["fundamental"]["runtime_sec"],
                "total_time": (
                    result["homography"]["runtime_sec"]
                    + result["fundamental"]["runtime_sec"]
                ),
                "h_inliers": result["homography"]["inliers"],
                "f_inliers": result["fundamental"]["inliers"],
                "h_signature": result["homography"]["signature"],
                "f_signature": result["fundamental"]["signature"],
            }
            rows.append(row)
            writer.writerow(
                [
                    row["run"],
                    row["seed"],
                    row["h_time"],
                    row["f_time"],
                    row["total_time"],
                    row["h_inliers"],
                    row["f_inliers"],
                    row["h_signature"],
                    row["f_signature"],
                ]
            )

            if args.json_dir is not None:
                json_path = args.json_dir / f"seed_{seed}.json"
                json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

    write_summary(rows, args.summary_output)


if __name__ == "__main__":
    main()
