"""Cross-check a compound table against CAS Common Chemistry, by CAS RN.

This is the primary validation step of the chemlib pipeline: CAS itself -
the body that assigns CAS numbers - looked up directly by CAS Registry
Number via:

    https://commonchemistry.cas.org/api/detail?cas_rn={cas}

This is a direct, unambiguous lookup by registry number, and CAS Common
Chemistry's response already includes synonyms - so it takes exactly one
HTTP request per distinct CAS number (pubchem_fallback needs two: one to
resolve the name, one for synonyms).

A row is CONFIRMED only if, for the substance CAS has on file for that RN:
  - name              (normalized) == CAS's name OR any of its synonyms
                        (normalized) - CAS's preferred name often differs
                        stylistically from the input's (same rationale as
                        the fallback step).
  - formula           (element counts, order-independent) == CAS's
                        molecularFormula.
  - molecular_weight  (within a small tolerance) == CAS's molecularMass.
The input's formula and weight must agree with what CAS reports for that
CAS number, every time - this is not optional/best-effort.

CAS Common Chemistry is a curated set of ~500k "common" substances (far
smaller than PubChem's ~100M+), so a meaningful not-found rate is expected,
especially for obscure compounds - that alone doesn't mean the input entry
is wrong, just unindexed here. Everything this step rejects (including
"not found") gets a second chance via pubchem_fallback, which reads this
step's output table and reprocesses the rows whose validation_status is
'rejected'.

Every input row survives into the output table: a rejected compound keeps
its row with `validation_status='rejected'` and a `reject_reason` list, so
"why did this compound disappear?" stays answerable without re-running the
lookups, and the rejection breakdown is a GROUP BY rather than a separate
file. Resumability comes from the on-disk JSONL cache (the HTTP lookups are
the expensive part), not from the output file growing row by row.

Requires a (free) API key from https://www.cas.org/services/commonchemistry-api,
read from the CAS_API_KEY environment variable (or passed to run()/--api-key).
Before the batch runs, one sanity-check lookup verifies the key works, so a
bad/missing key fails fast with a clear message instead of mislabeling every
row as "not found".
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import threading
import urllib.parse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from ardeco.chemlib import schema
from ardeco.chemlib._compound_check import (
    RateLimiter,
    formulas_match,
    get_json_with_retry,
    log_run,
    normalize_name,
    weights_match,
)
from ardeco.io import read_rows, write_table

API_BASE = "https://commonchemistry.cas.org/api"
USER_AGENT = "ardeco-chemlib/cas-validation (stdlib urllib)"
SANITY_CHECK_CAS = (
    "50-00-0"  # formaldehyde - used to verify the API key works before the batch runs
)

VALIDATION_SOURCE = "cas"


@dataclass
class CASResult:
    rn: str | None = None
    name: str | None = None
    formula: str | None = None
    molecular_weight: str | None = None
    smiles: str | None = None
    synonyms: list[str] | None = None
    error: str | None = None  # set when the lookup itself failed


def fetch_cas_detail(
    cas: str, api_key: str, timeout: float, max_retries: int, limiter: RateLimiter
) -> CASResult:
    url = f"{API_BASE}/detail?{urllib.parse.urlencode({'cas_rn': cas})}"
    headers = {"X-API-KEY": api_key, "User-Agent": USER_AGENT}
    data, error = get_json_with_retry(url, headers, timeout, max_retries, limiter)
    if error:
        return CASResult(
            error="cas_not_found_in_cas_common_chemistry"
            if error == "not_found"
            else error
        )

    d = data or {}
    return CASResult(
        rn=d.get("rn"),
        name=d.get("name"),
        formula=d.get("molecularFormula"),
        molecular_weight=d.get("molecularMass"),
        smiles=d.get("canonicalSmile"),
        synonyms=d.get("synonyms"),
    )


def verify_api_key(api_key: str, timeout: float) -> None:
    """Fail fast with a clear message if the key is missing/invalid, rather
    than silently mislabeling every row as "not found" for the whole run."""
    limiter = RateLimiter(1.0)
    result = fetch_cas_detail(
        SANITY_CHECK_CAS, api_key, timeout, max_retries=1, limiter=limiter
    )
    if result.error == "auth_error":
        raise SystemExit(
            "CAS Common Chemistry rejected the API key (401/403). Check CAS_API_KEY / --api-key."
        )
    if result.error and result.error != "cas_not_found_in_cas_common_chemistry":
        raise SystemExit(f"CAS Common Chemistry sanity check failed: {result.error}")


def load_cache(path: Path) -> dict[str, CASResult]:
    cache: dict[str, CASResult] = {}
    if not path.exists():
        return cache
    with path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            cas = rec.pop("cas")
            cache[cas] = CASResult(**rec)
    return cache


def name_matches(row_name: str, result: CASResult) -> bool:
    """True if `row_name` matches CAS's preferred name or any known synonym
    for this CAS RN."""
    target = normalize_name(row_name)
    candidates = [result.name, *(result.synonyms or [])]
    return any(c and normalize_name(c) == target for c in candidates)


def classify(row: dict[str, Any], result: CASResult) -> list[str]:
    """Return the list of reasons `row` is rejected; empty list = confirmed."""
    if result.error:
        return [result.error]
    reasons = []
    if not name_matches(row["name"] or "", result):
        reasons.append("name_mismatch")
    if not result.formula or not formulas_match(
        row["source_formula"] or "", result.formula
    ):
        reasons.append("formula_mismatch")
    if not result.molecular_weight or not weights_match(
        row["source_molecular_weight"], result.molecular_weight
    ):
        reasons.append("molecular_weight_mismatch")
    return reasons


def annotate(row: dict[str, Any], result: CASResult) -> dict[str, Any]:
    """Return `row` with this step's validation columns filled in.

    A confirmed row takes CAS's preferred name and canonical SMILES: the
    name has already been proven to match the input's (directly or through
    a synonym), so carrying both would just be a redundant second column.
    What the registry returned is kept in the `ref_*` columns either way, so
    a mismatch can be audited without re-querying.
    """
    reasons = classify(row, result)
    confirmed = not reasons
    return {
        **row,
        "name": (result.name or row["name"]) if confirmed else row["name"],
        "smiles": result.smiles if confirmed else None,
        "validation_source": VALIDATION_SOURCE if confirmed else None,
        "validation_status": "confirmed" if confirmed else "rejected",
        "reject_reason": schema.reasons(*reasons),
        "ref_id": result.rn,
        "ref_name": result.name,
        "ref_formula": result.formula,
        "ref_molecular_weight": result.molecular_weight,
        "ref_smiles": result.smiles,
    }


def process(
    rows: list[dict[str, Any]],
    cache_path: Path,
    cache: dict[str, CASResult],
    api_key: str,
    timeout: float,
    max_retries: int,
    rate: float,
    workers: int,
) -> list[dict[str, Any]]:
    """Resolve every distinct CAS number in `rows` against CAS Common
    Chemistry (or the on-disk cache) and return the rows annotated with this
    step's validation columns, in input order.

    Resolution is per distinct CAS number, annotation is per row: a registry
    number appearing on several rows costs one lookup. Progress through a
    long run is visible in two places - the counter on stderr, and the JSONL
    cache growing as lookups land, which is also what makes an interrupted
    run resume without re-fetching."""
    rows_by_cas: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        rows_by_cas.setdefault((row["cas_number"] or "").strip(), []).append(row)
    distinct_cas = list(rows_by_cas)
    n_to_fetch = sum(1 for cas in distinct_cas if cas not in cache)
    print(
        f"{len(distinct_cas)} distinct CAS numbers ({n_to_fetch} need a CAS Common Chemistry lookup, "
        f"{len(distinct_cas) - n_to_fetch} already cached) across {len(rows)} rows"
    )

    limiter = RateLimiter(rate)
    cache_lock = threading.Lock()

    def resolve(cas: str) -> CASResult:
        cached = cache.get(cas)
        if cached is not None:
            return cached
        result = fetch_cas_detail(cas, api_key, timeout, max_retries, limiter)
        with cache_lock:
            with cache_path.open("a", encoding="utf-8") as cache_file:
                cache_file.write(
                    json.dumps({"cas": cas, **asdict(result)}, ensure_ascii=False)
                    + "\n"
                )
        return result

    resolved: dict[str, CASResult] = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for i, (cas, result) in enumerate(
            zip(distinct_cas, pool.map(resolve, distinct_cas)), 1
        ):
            resolved[cas] = result
            if i % 100 == 0 or i == len(distinct_cas):
                print(f"[{i}/{len(distinct_cas)}] resolved", file=sys.stderr)

    annotated = [
        annotate(row, resolved[(row["cas_number"] or "").strip()]) for row in rows
    ]
    n_confirmed = sum(1 for r in annotated if r["validation_status"] == "confirmed")
    print(f"Done. confirmed={n_confirmed} rejected={len(annotated) - n_confirmed}")
    return annotated


def run(
    input_path: Path,
    output_path: Path,
    cache_path: Path,
    log_path: Path,
    api_key: str,
    *,
    rate: float = 4.0,
    workers: int = 4,
    timeout: float = 15.0,
    max_retries: int = 4,
    limit: int | None = None,
    argv: list[str] | None = None,
) -> tuple[int, int]:
    """Validate `input_path` against CAS Common Chemistry, writing the
    annotated table to `output_path`, and return (n_confirmed, n_rejected)."""
    verify_api_key(api_key, timeout)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    log_run(
        log_path,
        event="start",
        command=" ".join(argv if argv is not None else sys.argv),
        input=str(input_path),
        output=str(output_path),
        cache=str(cache_path),
        rate=rate,
        workers=workers,
        limit=limit,
    )

    rows = read_rows(input_path)
    if limit:
        rows = rows[:limit]

    cache = load_cache(cache_path)
    print(f"Loaded {len(cache)} cached CAS Common Chemistry lookups from {cache_path}")

    annotated = process(
        rows, cache_path, cache, api_key, timeout, max_retries, rate, workers
    )
    write_table(annotated, output_path, schema.VALIDATED, stage="02_cas")

    n_confirmed = sum(1 for r in annotated if r["validation_status"] == "confirmed")
    n_rejected = len(annotated) - n_confirmed
    log_run(
        log_path,
        event="done",
        command=" ".join(argv if argv is not None else sys.argv),
        output=str(output_path),
        n_confirmed=n_confirmed,
        n_rejected=n_rejected,
    )
    return n_confirmed, n_rejected


def add_subparser(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
) -> None:
    """Wire this module's CLI as the `cas_validate` subcommand of
    `python -m ardeco.chemlib` (see chemlib/__main__.py)."""
    parser = subparsers.add_parser(
        "cas_validate",
        help="validate a compound table against CAS Common Chemistry",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    _ = parser.add_argument("input", type=Path, help="raw compound Parquet table")
    _ = parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="output Parquet table (all rows, tagged)",
    )
    _ = parser.add_argument(
        "--cache",
        type=Path,
        required=True,
        help="JSONL cache of CAS lookups keyed by CAS number",
    )
    _ = parser.add_argument(
        "--log", type=Path, help="run_log.jsonl path (default: alongside --output)"
    )
    _ = parser.add_argument(
        "--api-key", help="CAS Common Chemistry API key (default: $CAS_API_KEY)"
    )
    _ = parser.add_argument(
        "--rate",
        type=float,
        default=4.0,
        help="max aggregate requests/second (default: 4.0)",
    )
    _ = parser.add_argument(
        "--workers", type=int, default=4, help="concurrent lookup threads (default: 4)"
    )
    _ = parser.add_argument(
        "--timeout", type=float, default=15.0, help="per-request timeout in seconds"
    )
    _ = parser.add_argument(
        "--max-retries", type=int, default=4, help="retries on 429/5xx/network errors"
    )
    _ = parser.add_argument(
        "--limit", type=int, help="only process the first N rows (for a dry run)"
    )
    parser.set_defaults(_run=_cli_run)


def _load_dotenv(path: Path = Path(".env")) -> None:
    """Load KEY=VALUE lines from `path` into the environment (existing
    env vars win). Minimal by design - CAS_API_KEY is the only secret this
    project reads today, so a dependency for a handful of lines isn't
    worth it. Expects to run from the repo root, same as Snakemake."""
    if not path.exists():
        return
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        os.environ.setdefault(key.strip(), value.strip().strip('"').strip("'"))


def _cli_run(args: argparse.Namespace) -> None:
    _load_dotenv()
    api_key = args.api_key or os.environ.get("CAS_API_KEY")
    if not api_key:
        raise SystemExit(
            "No API key. Set the CAS_API_KEY environment variable or pass --api-key "
            "(get one at https://www.cas.org/services/commonchemistry-api)."
        )
    log_path = args.log or args.output.parent / "run_log.jsonl"
    run(
        args.input,
        args.output,
        args.cache,
        log_path,
        api_key,
        rate=args.rate,
        workers=args.workers,
        timeout=args.timeout,
        max_retries=args.max_retries,
        limit=args.limit,
    )
    print(f"Output  -> {args.output}")
    print(f"Cache   -> {args.cache}")
    print(f"Run log -> {log_path}")
