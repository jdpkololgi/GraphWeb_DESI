#!/usr/bin/env python3
"""Read-only, bounded P12-A handoff verification; never opens catalogue payloads.

Writes an exclusive new report file. Exit 2 means handoff is incomplete, not a
scientific calibration failure. No checkpoint deserialization or inference.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path


MAX_BYTES = 32 * 1024 * 1024
FEATURES = ["base_lambda1", "base_lambda2", "base_lambda3", "redshift",
            "log_ntilde_mpc3", "cap_ngc", "log1p_random_support_boundary_distance_mpc"]
REQUIRED = {"base_checkpoint", "checkpoint", "base_summary", "completion",
            "dataset", "calibration_audit", "gaussian_baseline", "quality_thresholds"}


def digest(path: Path) -> str:
    """Bound reads even if a file grows after stat."""
    h = hashlib.sha256()
    total = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            total += len(chunk)
            if total > MAX_BYTES:
                raise ValueError("small-artifact read limit exceeded")
            h.update(chunk)
    return h.hexdigest()


def read_json(path: Path) -> dict:
    if path.stat().st_size > MAX_BYTES:
        raise ValueError("JSON exceeds small-artifact limit")
    with path.open("rb") as handle:
        raw = handle.read(MAX_BYTES + 1)
    if len(raw) > MAX_BYTES:
        raise ValueError("JSON grew beyond small-artifact limit")
    value = json.loads(raw)
    if not isinstance(value, dict):
        raise ValueError("expected a JSON object")
    return value


def verify_record(record: dict, repo: Path) -> dict:
    path = Path(record.get("path", ""))
    result = {"recorded": record, "status": "invalid_record", "verified": False}
    if (not record.get("path") or not isinstance(record.get("bytes"), int)
            or record["bytes"] < 0
            or not re.fullmatch(r"[0-9a-f]{64}", str(record.get("sha256", "")))):
        return result
    if not path.is_absolute():
        path = repo / path
    result["resolved_path"] = str(path)
    try:
        if record["bytes"] > MAX_BYTES or path.stat().st_size > MAX_BYTES:
            result["status"] = "over_read_limit"
        elif path.stat().st_size != record["bytes"]:
            result["status"] = "size_mismatch"
        else:
            actual = digest(path)
            result.update(actual_sha256=actual, verified=actual == record["sha256"])
            result["status"] = "verified" if result["verified"] else "hash_mismatch"
    except (OSError, ValueError) as error:
        result.update(status="unavailable", error=str(error))
    return result


def catalogue_records(value, role=""):
    """Extract archived metadata only; never stat/hash the large FITS files."""
    if isinstance(value, dict):
        if "path" in value and "columns" in value:
            yield {"role": role, **{key: value.get(key) for key in
                   ("path", "columns", "rows", "bytes", "sha256")},
                   "live_verified": False}
        else:
            for key, child in value.items():
                yield from catalogue_records(child, f"{role}/{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from catalogue_records(child, f"{role}/{index}")


def phase_ledger(p10: dict, e2e: dict) -> list:
    roles = p10["model_phase_contract"]
    phases = set(e2e["phase_roles"])
    for key in ("training", "validation_and_selection", "sealed_blind_test"):
        phases.update(roles[key])
    rows = []
    for phase in sorted(phases):
        p12_role = "not_in_original_p12_contract"
        if phase in roles["training"]:
            p12_role = "training"
        elif phase in roles["validation_and_selection"]:
            p12_role = "selection_repeatedly_inspected"
        elif phase in roles["sealed_blind_test"]:
            p12_role = "historical_blind_opened_once_20260905"
        e2e_role = e2e["phase_roles"].get(phase, "excluded_or_historical")
        rows.append({"phase": phase, "p12_role": p12_role, "e2e_role": e2e_role,
                     "reserved_confirmation": e2e_role == "confirmation",
                     "all_programme_access_history_verified": False,
                     "fresh_blind_eligible": False,
                     "payload_access_authorized_by_this_report": False,
                     "replication_candidate_pending_exposure_audit":
                         p12_role == "not_in_original_p12_contract" and e2e_role == "train"})
    return rows


def build_report(repo: Path) -> dict:
    inputs = {
        "candidate": repo / "docs/evidence/p12/P12A_PRODUCTION_CANDIDATE_FROZEN.json",
        "sources": repo / "configs/p10_response_sources_v1.json",
        "p10": repo / "configs/p10_phase_registry_v1.json",
        "e2e": repo / "configs/e2e_coupled_data_v1.json",
    }
    data = {name: read_json(path) for name, path in inputs.items()}
    candidate = data["candidate"]
    if candidate.get("schema_version") != "p12a-production-candidate-frozen-v1":
        raise ValueError("unsupported candidate schema")
    if candidate.get("pass") is not True:
        raise ValueError("candidate was not frozen successfully")
    if candidate.get("tempering") is not None or candidate.get("recalibration") is not None:
        raise ValueError("expected untempered, uncorrected P12-A")
    artifacts = candidate["artifacts"]
    missing = sorted(REQUIRED - artifacts.keys())
    checks = {name: verify_record(artifacts[name], repo)
              for name in sorted(REQUIRED & artifacts.keys())}
    blockers = [f"missing_artifact:{name}" for name in missing]
    blockers += [f"artifact:{name}:{value['status']}" for name, value in checks.items()
                 if not value["verified"]]
    if checks.get("completion", {}).get("verified"):
        completion = read_json(Path(checks["completion"]["resolved_path"]))
        if completion.get("conditioning_features") != FEATURES:
            blockers.append("conditioning_feature_order_mismatch")
        if completion.get("checkpoint_sha256") != artifacts["checkpoint"]["sha256"]:
            blockers.append("completion_checkpoint_binding_mismatch")
    if candidate.get("base_encoder", {}).get("checkpoint") != artifacts.get("base_checkpoint"):
        blockers.append("encoder_checkpoint_binding_mismatch")
    desi = data["sources"]["desi_candidate"]
    expected_release = "DA2/loa-v1/LSScats/v2.1/PIP"
    if desi.get("release") != expected_release or desi.get("deployment_family") != "loa-v1":
        blockers.append("unexpected_deployment_release")
    catalogues = list(catalogue_records(desi))
    allowed_root = "/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2.1/"
    if not catalogues or any(not str(row["path"]).startswith(allowed_root) for row in catalogues):
        blockers.append("catalogue_release_root_mismatch")
    # None of these scientific/deployment gates can be passed by metadata checks.
    outstanding = ["frozen_encoder_and_oof_transform_lineage",
                   "coordinate_and_native_host_label_closure",
                   "loa_live_source_refreeze_and_redshift_quality_policy",
                   "response_field_crosswalk", "final_checkpoint_context_parity",
                   "golden_mock_end_to_end_replay"]
    return {
        "schema_version": "p12a-vac-preflight-v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "bounded artifact verification and archived source metadata only",
        "inputs": {name: {"path": str(path), "sha256": digest(path)}
                   for name, path in inputs.items()},
        "artifact_checks": checks, "artifact_blockers": blockers,
        "artifact_inventory_pass": not blockers,
        "conditioning_features": FEATURES,
        "recorded_quality_bits": candidate["quality_bits"],
        "recorded_posterior_draws": candidate["posterior_draws"],
        "coordinate_audit_status": "unresolved",
        "catalogue_inventory": catalogues,
        "catalogue_hashes_are_archived_not_live_verified": True,
        "phase_exposure_ledger": phase_ledger(data["p10"], data["e2e"]),
        "outstanding_gates": outstanding,
        "ready_for_desi_canary": False, "science_release_ready": False,
        "truth_payloads_read": [], "catalogue_payloads_read": [],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--illustris-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = build_report(args.illustris_root.resolve())
    # Exclusive creation prevents overwriting an earlier audit. A failed write
    # leaves an invalid report, never a success marker or a launch license.
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({"output": str(args.output),
                      "artifact_inventory_pass": report["artifact_inventory_pass"],
                      "ready_for_desi_canary": False}))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
