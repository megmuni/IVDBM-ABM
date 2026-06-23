#!/usr/bin/env python3
"""Write test_summary.json from JUnit output after a Slurm ctest job."""

from __future__ import annotations

import argparse
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any


def parse_junit(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}

    root = ET.parse(path).getroot()
    summary: dict[str, Any] = {
        "tests": 0,
        "failures": 0,
        "errors": 0,
        "skipped": 0,
        "time_seconds": 0.0,
        "cases": [],
    }

    for suite in root.iter("testsuite"):
        summary["tests"] += int(suite.attrib.get("tests", 0))
        summary["failures"] += int(suite.attrib.get("failures", 0))
        summary["errors"] += int(suite.attrib.get("errors", 0))
        summary["skipped"] += int(suite.attrib.get("skipped", 0))
        summary["time_seconds"] += float(suite.attrib.get("time", 0.0) or 0.0)
        for case in suite.findall("testcase"):
            entry = {
                "name": case.attrib.get("name"),
                "classname": case.attrib.get("classname"),
                "time_seconds": float(case.attrib.get("time", 0.0) or 0.0),
                "status": "passed",
            }
            if case.find("failure") is not None:
                entry["status"] = "failed"
            elif case.find("error") is not None:
                entry["status"] = "error"
            elif case.find("skipped") is not None:
                entry["status"] = "skipped"
            summary["cases"].append(entry)

    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--junit", type=Path, default=None)
    parser.add_argument("--suite", default="")
    parser.add_argument("--profile", default="")
    parser.add_argument("--build-dir", default="build")
    parser.add_argument("--exit-code", type=int, default=0)
    parser.add_argument("--job-id", default="")
    args = parser.parse_args()

    junit_path = args.junit or (args.output_dir / "test_results.xml")
    junit_summary = parse_junit(junit_path)

    payload = {
        "slurm_job_id": args.job_id or None,
        "suite": args.suite or None,
        "profile": args.profile or None,
        "build_dir": args.build_dir,
        "exit_code": args.exit_code,
        "junit_file": str(junit_path),
        "results": junit_summary,
    }

    out_path = args.output_dir / "test_summary.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")

    if junit_summary:
        print(
            "Test summary: "
            f"{junit_summary.get('tests', 0)} tests, "
            f"{junit_summary.get('failures', 0)} failures, "
            f"{junit_summary.get('errors', 0)} errors"
        )
    else:
        print("Test summary: no JUnit results parsed")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
