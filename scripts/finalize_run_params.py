#!/usr/bin/env python3
"""Merge job runtime statistics into run_params.json."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def merge_dict(base: dict[str, Any], extra: dict[str, Any]) -> dict[str, Any]:
    for key, value in extra.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            merge_dict(base[key], value)
        else:
            base[key] = value
    return base


def parse_time_file(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}

    parsed: dict[str, Any] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip()
        if not key:
            continue
        try:
            parsed[key] = float(value)
        except ValueError:
            parsed[key] = value
    return parsed


def parse_sacct(job_id: str) -> dict[str, Any]:
    if not job_id:
        return {}

    command = [
        "sacct",
        "-j",
        f"{job_id}.batch",
        "--format=Elapsed,TotalCPU,MaxRSS,MaxVMSize,State,ExitCode",
        "-P",
        "-n",
        "--units=M",
    ]
    try:
        result = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError:
        return {}

    line = result.stdout.strip().splitlines()
    if not line:
        return {}

    fields = line[0].split("|")
    if len(fields) < 6:
        return {}

    elapsed, total_cpu, max_rss, max_vms, state, exit_code = fields[:6]
    accounting: dict[str, Any] = {
        "elapsed": elapsed or None,
        "total_cpu": total_cpu or None,
        "max_rss": max_rss or None,
        "max_vms": max_vms or None,
        "state": state or None,
        "exit_code": exit_code or None,
    }
    if max_rss and max_rss.isdigit():
        accounting["max_rss_mb"] = int(max_rss)
    if max_vms and max_vms.isdigit():
        accounting["max_vms_mb"] = int(max_vms)
    return accounting


def build_runtime_summary(args: argparse.Namespace) -> dict[str, Any]:
    runtime: dict[str, Any] = {
        "started_at": args.started or None,
        "finished_at": args.finished or None,
        "hostname": args.hostname or None,
        "profile": args.profile or None,
        "exit_code": args.exit_code,
    }

    testrun: dict[str, Any] = {}
    if args.time_file is not None:
        testrun = parse_time_file(args.time_file)
    if testrun:
        runtime["testrun"] = testrun

    accounting = parse_sacct(args.sacct_job_id or "")
    if accounting:
        runtime["slurm_accounting"] = accounting

    return runtime


def print_human_summary(runtime: dict[str, Any]) -> None:
    print("Job summary:")
    if runtime.get("started_at"):
        print(f"  started:  {runtime['started_at']}")
    if runtime.get("finished_at"):
        print(f"  finished: {runtime['finished_at']}")
    if runtime.get("hostname"):
        print(f"  host:     {runtime['hostname']}")
    if runtime.get("profile"):
        print(f"  profile:  {runtime['profile']}")
    print(f"  exit:     {runtime.get('exit_code')}")

    testrun = runtime.get("testrun", {})
    if testrun.get("wall_seconds") is not None:
        print(f"  testRun wall time: {testrun['wall_seconds']:.1f} s")
    if testrun.get("max_rss_kb") is not None:
        print(
            "  testRun max RSS:   "
            f"{testrun['max_rss_kb'] / 1024:.1f} MiB"
        )

    sim = runtime.get("tick_timing", {})
    if sim.get("tick_execution_ms_avg") is not None:
        print(
            "  avg tick time:     "
            f"{sim['tick_execution_ms_avg']:.1f} ms"
        )

    accounting = runtime.get("slurm_accounting", {})
    if accounting.get("elapsed"):
        print(f"  Slurm elapsed:     {accounting['elapsed']}")
    if accounting.get("total_cpu"):
        print(f"  Slurm CPU time:    {accounting['total_cpu']}")
    if accounting.get("max_rss_mb") is not None:
        print(f"  Slurm max RSS:     {accounting['max_rss_mb']} MiB")
    elif accounting.get("max_rss"):
        print(f"  Slurm max RSS:     {accounting['max_rss']}")
    if accounting.get("state"):
        print(f"  Slurm state:       {accounting['state']}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_params", type=Path, help="Path to run_params.json")
    parser.add_argument(
        "--timing",
        type=Path,
        default=None,
        help="Optional simulation timing sidecar JSON from testRun",
    )
    parser.add_argument("--time-file", type=Path, default=None)
    parser.add_argument("--sacct-job-id", default="")
    parser.add_argument("--started", default="")
    parser.add_argument("--finished", default="")
    parser.add_argument("--hostname", default="")
    parser.add_argument("--profile", default="")
    parser.add_argument("--exit-code", type=int, default=0)
    parser.add_argument(
        "--quiet-summary",
        action="store_true",
        help="Do not print a human-readable summary to stdout",
    )
    args = parser.parse_args()

    data = load_json(args.run_params)
    runtime = build_runtime_summary(args)

    timing: dict[str, Any] = {}
    if args.timing is not None and args.timing.is_file():
        timing = load_json(args.timing)
    if timing:
        runtime["tick_timing"] = timing

    if runtime:
        data["runtime"] = runtime

    args.run_params.parent.mkdir(parents=True, exist_ok=True)
    with args.run_params.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2)
        handle.write("\n")

    if args.timing is not None and args.timing.is_file():
        args.timing.unlink()
    if args.time_file is not None and args.time_file.is_file():
        args.time_file.unlink()

    if not args.quiet_summary:
        print_human_summary(runtime)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
