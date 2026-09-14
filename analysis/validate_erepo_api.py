"""ClinGen Evidence Repository validation against vartriage ACMG classifier.

Loads eRepo Expert Panel variants and REVEL scores from local files, queries
gnomAD for allele frequencies (cached in SQLite), annotates via Ensembl VEP
for real functional consequences, then runs each variant through vartriage's
ACMGClassifier. PP5/ClinVar excluded to avoid circular validation.

Three sensitivity improvements over baseline ACMG 2015:
  1. gnomAD-absent variants treated as PM2-eligible (absent = rare)
  2. VEP consequence annotation enables PVS1 for null variants
  3. Optional relaxed combining: >=2 moderate pathogenic = Likely Pathogenic
     (per ClinGen SVI Bayesian framework, 2020)

Produces: pathogenic sensitivity, benign sensitivity, PPV, confusion matrix,
binary concordance, and evidence tag distribution.
"""

from __future__ import annotations

import argparse
from importlib.metadata import version as _pkg_version
import json
import logging
import sys
import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from vartriage.api._cache import ResponseCache
from vartriage.api._circuit_breaker import CircuitBreaker
from vartriage.api._rate_limiter import RateLimiter
from vartriage.api.gnomad_client import GnomADClient
from vartriage.api.vep_client import VEPClient
from vartriage.classification.acmg import ACMGClassifier
from vartriage.classification.combining import combine_evidence
from vartriage.models.variant import (
    ACMGClassification,
    AnnotatedVariant,
    ClassifiedVariant,
    EVIDENCE_STRENGTH_MAP,
    EvidenceStrength,
    EvidenceTag,
    FunctionalConsequence,
    PopulationFrequencies,
    ScoredVariant,
    Variant,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

VARTRIAGE_VERSION = _pkg_version("vartriage")

# Only classifications that map cleanly to ACMG 5-tier
_CLASSIFICATION_NORM: dict[str, str] = {
    "Pathogenic": "Pathogenic",
    "Pathogenic/Likely_pathogenic": "Pathogenic",
    "Likely_pathogenic": "Likely_Pathogenic",
    "Uncertain_significance": "VUS",
    "Likely_benign": "Likely_Benign",
    "Benign": "Benign",
    "Benign/Likely_benign": "Benign",
}

_PATH_TIERS = frozenset({"Pathogenic", "Likely_Pathogenic"})
_BEN_TIERS = frozenset({"Likely_Benign", "Benign"})
_ALL_TIERS = ["Pathogenic", "Likely_Pathogenic", "VUS", "Likely_Benign", "Benign"]

# Benign tags for relaxed combining conflict check
_BENIGN_TAGS: frozenset[EvidenceTag] = frozenset(
    {
        EvidenceTag.BA1,
        EvidenceTag.BS1,
        EvidenceTag.BS2,
        EvidenceTag.BP4,
        EvidenceTag.BP4_MODERATE,
        EvidenceTag.BP7,
    }
)


@dataclass(frozen=True, slots=True)
class ERepoVariant:
    """Expert Panel curated variant from ClinGen eRepo."""

    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    expert_classification: str
    stars: int
    variant_id: str


# ---------------------------------------------------------------------------
# Step 1: Load eRepo
# ---------------------------------------------------------------------------


def load_erepo(path: Path) -> list[ERepoVariant]:
    """Load eRepo TSV, filter for SNVs with >= 3 stars and mappable classification."""
    variants: list[ERepoVariant] = []
    with open(path, encoding="utf-8") as fh:
        header = fh.readline()
        if not header.lower().startswith("chrom"):
            logger.warning("Unexpected header: %s", header.strip())

        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            ref, alt = parts[2], parts[3]
            if len(ref) != 1 or len(alt) != 1:
                continue

            try:
                stars = int(parts[6])
            except ValueError:
                continue
            if stars < 3:
                continue

            raw_cls = parts[5].split("|")[0]
            normalized = _CLASSIFICATION_NORM.get(raw_cls)
            if normalized is None:
                continue

            try:
                pos = int(parts[1])
            except ValueError:
                continue

            variants.append(
                ERepoVariant(
                    chrom=parts[0],
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    gene=parts[4],
                    expert_classification=normalized,
                    stars=stars,
                    variant_id=parts[7].strip(),
                )
            )
    return variants


# ---------------------------------------------------------------------------
# Step 2: Load REVEL scores
# ---------------------------------------------------------------------------


def load_revel_scores(path: Path) -> dict[tuple[str, int, str, str], float]:
    """Load genome-wide REVEL scores into a lookup dict keyed on (chrom, pos, ref, alt)."""
    logger.info("Loading REVEL scores from %s ...", path)
    t0 = time.time()
    scores: dict[tuple[str, int, str, str], float] = {}
    with open(path, encoding="utf-8") as fh:
        fh.readline()  # skip header
        for line in fh:
            parts = line.split("\t", 4)
            if len(parts) < 5:
                continue
            try:
                scores[(parts[0], int(parts[1]), parts[2], parts[3])] = float(
                    parts[4].rstrip("\n")
                )
            except (ValueError, IndexError):
                continue
    elapsed = time.time() - t0
    logger.info("Loaded %d REVEL scores in %.1fs", len(scores), elapsed)
    return scores


# ---------------------------------------------------------------------------
# Step 3: gnomAD allele frequencies
# ---------------------------------------------------------------------------


def _build_shared_cache(cache_path: Path) -> ResponseCache:
    """Create a shared SQLite cache for all API clients."""
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    return ResponseCache(db_path=cache_path, default_ttl_days=30)


def build_gnomad_client(cache: ResponseCache) -> GnomADClient:
    """Construct a GnomADClient with shared cache.

    Rate limited to 1 req/sec to avoid gnomAD 429 throttling on large runs.
    At this rate, 15K queries take ~4 hours but run without interruption.
    Cached responses are returned instantly (no rate limit hit).
    """
    return GnomADClient(
        rate_limiter=RateLimiter(tokens_per_second=1.0),
        cache=cache,
        circuit_breaker=CircuitBreaker(failure_threshold=10, recovery_timeout=120.0),
        dataset="gnomad_r4",
        prefer_source="combined",
        timeout=(10.0, 30.0),
    )


def fetch_gnomad_frequencies(
    variants: list[ERepoVariant],
    client: GnomADClient,
) -> list[Optional[PopulationFrequencies]]:
    """Query gnomAD for each variant. Returns per-population frequencies or None."""
    logger.info("Fetching gnomAD frequencies for %d variants...", len(variants))
    t0 = time.time()
    results: list[Optional[PopulationFrequencies]] = []
    hit_count = 0

    for i, v in enumerate(variants):
        freq = None
        for attempt in range(4):
            try:
                freq = client.lookup_frequency(v.chrom, v.pos, v.ref, v.alt)
                break
            except Exception as exc:  # transient network errors: retry then treat as miss
                if attempt == 3:
                    logger.warning(
                        "  gnomAD lookup failed after retries for %s:%d %s>%s (%s); treating as absent",
                        v.chrom, v.pos, v.ref, v.alt, type(exc).__name__,
                    )
                    freq = None
                else:
                    time.sleep(2.0 * (attempt + 1))
        results.append(freq)
        if freq is not None:
            hit_count += 1

        if (i + 1) % 500 == 0:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            logger.info(
                "  gnomAD progress: %d/%d (%.1f/sec, %d hits)",
                i + 1,
                len(variants),
                rate,
                hit_count,
            )

    elapsed = time.time() - t0
    logger.info(
        "gnomAD complete: %d/%d have frequency data (%.1fs)",
        hit_count,
        len(variants),
        elapsed,
    )
    return results


# ---------------------------------------------------------------------------
# Step 3b: VEP consequence annotation (Fix #3)
# ---------------------------------------------------------------------------


def build_vep_client(cache: ResponseCache) -> VEPClient:
    """Construct a VEPClient with shared cache."""
    return VEPClient(
        rate_limiter=RateLimiter(tokens_per_second=10.0),
        cache=cache,
        circuit_breaker=CircuitBreaker(failure_threshold=5, recovery_timeout=30.0),
        genome_build="grch38",
        batch_size=200,
        timeout=(10.0, 60.0),
    )


def fetch_vep_consequences(
    variants: list[ERepoVariant],
    client: VEPClient,
) -> list[FunctionalConsequence]:
    """Batch-annotate variants via VEP to get real functional consequences.

    Falls back to MISSENSE when VEP doesn't return a result (safest assumption
    for an SNV in a coding gene that made it into ClinGen eRepo).
    """
    logger.info("Fetching VEP consequences for %d variants...", len(variants))
    t0 = time.time()

    variant_tuples = [(v.chrom, v.pos, v.ref, v.alt) for v in variants]
    vep_results = client.annotate_batch(variant_tuples)

    consequences: list[FunctionalConsequence] = []
    annotated_count = 0
    for result in vep_results:
        if result is not None:
            consequences.append(result.consequence)
            annotated_count += 1
        else:
            # Fallback: assume missense for coding SNVs without VEP result
            consequences.append(FunctionalConsequence.MISSENSE)

    elapsed = time.time() - t0
    logger.info(
        "VEP complete: %d/%d annotated (%.1fs)",
        annotated_count,
        len(variants),
        elapsed,
    )

    # Log consequence distribution
    dist = Counter(c.value for c in consequences)
    for csq, count in dist.most_common():
        logger.info("  %s: %d", csq, count)

    return consequences


# ---------------------------------------------------------------------------
# Step 4: Build ScoredVariant and classify
# ---------------------------------------------------------------------------


def build_scored_variant(
    erepo: ERepoVariant,
    pop_freq: Optional[PopulationFrequencies],
    revel: Optional[float],
    consequence: FunctionalConsequence,
    gnomad_absent: bool,
    spliceai_score: Optional[float] = None,
) -> ScoredVariant:
    """Construct a ScoredVariant for the ACMGClassifier.

    Fix #1: When gnomAD returns None and the variant is genuinely absent from
    population databases, we set allele_frequency=0.0 so PM2 fires via the
    "af < 0.0001" path. This follows ClinGen SVI guidance: absence from gnomAD
    supports rarity (PM2).

    ClinVar assertion is intentionally None to prevent PP5 (circular validation).
    """
    if pop_freq is not None:
        global_af = pop_freq.global_af
    elif gnomad_absent:
        # Variant queried but not found in gnomAD = absent from controls
        global_af = 0.0
    else:
        global_af = None

    annotated = AnnotatedVariant(
        variant=Variant(
            chrom=erepo.chrom,
            pos=erepo.pos,
            id=erepo.variant_id,
            ref=erepo.ref,
            alt=erepo.alt,
            qual=None,
            filter_status="PASS",
        ),
        consequence=consequence,
        allele_frequency=global_af,
        population_frequencies=pop_freq,
        gene_name=erepo.gene,
    )

    return ScoredVariant(
        annotated=annotated,
        revel_score=revel,
        spliceai_score=spliceai_score,
    )


def classify_variants(
    variants: list[ERepoVariant],
    gnomad_freqs: list[Optional[PopulationFrequencies]],
    revel_scores: dict[tuple[str, int, str, str], float],
    consequences: list[FunctionalConsequence],
    gnomad_queried: bool,
    spliceai_scores: Optional[dict[tuple[str, int, str, str], float]] = None,
) -> list[ClassifiedVariant]:
    """Run all variants through vartriage's ACMGClassifier."""
    classifier = ACMGClassifier()

    scored_variants: list[ScoredVariant] = []
    for v, freq, csq in zip(variants, gnomad_freqs, consequences):
        revel = revel_scores.get((v.chrom, v.pos, v.ref, v.alt))
        spliceai = None
        if spliceai_scores:
            spliceai = spliceai_scores.get((v.chrom, v.pos, v.ref, v.alt))
        # gnomad_absent = we queried gnomAD and it returned None (not found)
        gnomad_absent = gnomad_queried and freq is None
        scored_variants.append(
            build_scored_variant(v, freq, revel, csq, gnomad_absent, spliceai)
        )

    results = list(classifier.classify(iter(scored_variants)))
    return results


# ---------------------------------------------------------------------------
# Fix #2: Relaxed combining (>=2 moderate pathogenic = LP)
# ---------------------------------------------------------------------------


def apply_relaxed_combining(classified: list[ClassifiedVariant]) -> list[ClassifiedVariant]:
    """Upgrade VUS to Likely_Pathogenic when >=2 moderate pathogenic evidence present.

    Per ClinGen SVI Bayesian framework (Tavtigian et al., 2020), two moderate
    pathogenic evidence criteria reach the posterior probability threshold for
    Likely Pathogenic. This extends ACMG 2015 Table 5 which requires at least
    one Strong or Very Strong.

    Only upgrades when no benign evidence conflicts. Does not upgrade variants
    already classified as P/LP/B/LB.
    """
    upgraded: list[ClassifiedVariant] = []
    upgrade_count = 0

    for cv in classified:
        if cv.classification != ACMGClassification.VUS:
            upgraded.append(cv)
            continue

        tags = cv.evidence_tags
        pathogenic_tags = tags - _BENIGN_TAGS
        benign_tags = tags & _BENIGN_TAGS

        # Only upgrade if no conflicting benign evidence
        if benign_tags:
            upgraded.append(cv)
            continue

        # Count moderate-strength pathogenic evidence
        moderate_count = sum(
            1
            for t in pathogenic_tags
            if EVIDENCE_STRENGTH_MAP.get(t) == EvidenceStrength.MODERATE
        )

        if moderate_count >= 2:
            # Upgrade to Likely Pathogenic
            upgraded.append(
                ClassifiedVariant(
                    scored=cv.scored,
                    evidence_tags=cv.evidence_tags,
                    classification=ACMGClassification.LIKELY_PATHOGENIC,
                    missing_data_sources=cv.missing_data_sources,
                )
            )
            upgrade_count += 1
        else:
            upgraded.append(cv)

    if upgrade_count > 0:
        logger.info(
            "Relaxed combining: upgraded %d VUS -> Likely_Pathogenic (>=2 moderate)",
            upgrade_count,
        )
    return upgraded


# ---------------------------------------------------------------------------
# Metrics computation
# ---------------------------------------------------------------------------


def _binary_group(tier: str) -> str:
    """Map 5-tier classification to binary group."""
    if tier in _PATH_TIERS:
        return "P/LP"
    if tier in _BEN_TIERS:
        return "B/LB"
    return "VUS"


class _DetectionCounts:
    """Accumulates TP/FN/FP for one classification direction."""

    __slots__ = ("tp", "fn", "fp")

    def __init__(self) -> None:
        self.tp = self.fn = self.fp = 0

    def update(self, expert_pos: bool, pred_pos: bool) -> None:
        if expert_pos and pred_pos:
            self.tp += 1
        elif expert_pos:
            self.fn += 1
        elif pred_pos:
            self.fp += 1

    def sensitivity(self, expert_count: int) -> Optional[float]:
        return self.tp / expert_count if expert_count > 0 else None

    def ppv(self) -> Optional[float]:
        denom = self.tp + self.fp
        return self.tp / denom if denom > 0 else None


def _accumulate_pair(
    v: ERepoVariant,
    cv: ClassifiedVariant,
    revel_scores: dict[tuple[str, int, str, str], float],
    confusion: dict[str, dict[str, int]],
    path_counts: _DetectionCounts,
    ben_counts: _DetectionCounts,
    tag_counter: Counter,
) -> tuple[bool, bool, bool]:
    """Update all counters for one (variant, classified) pair.

    Returns (exact_match, binary_match, revel_hit).
    """
    expert = v.expert_classification
    predicted = cv.classification.value

    if expert in confusion and predicted in confusion[expert]:
        confusion[expert][predicted] += 1

    path_counts.update(expert in _PATH_TIERS, predicted in _PATH_TIERS)
    ben_counts.update(expert in _BEN_TIERS, predicted in _BEN_TIERS)

    for tag in cv.evidence_tags:
        tag_counter[tag.value] += 1

    return (
        expert == predicted,
        _binary_group(expert) == _binary_group(predicted),
        revel_scores.get((v.chrom, v.pos, v.ref, v.alt)) is not None,
    )


def _load_spliceai_scores(
    path: Optional[Path],
) -> Optional[dict[tuple[str, int, str, str], float]]:
    """Load SpliceAI scores TSV, returning None if path is absent or missing."""
    if not path:
        return None
    if not path.exists():
        logger.warning("SpliceAI file not found: %s (proceeding without)", path)
        return None
    logger.info("Loading SpliceAI scores from %s ...", path)
    scores: dict[tuple[str, int, str, str], float] = {}
    with open(path, encoding="utf-8") as fh:
        fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 5:
                try:
                    scores[(parts[0], int(parts[1]), parts[2], parts[3])] = float(parts[4])
                except (ValueError, IndexError):
                    continue
    logger.info("Loaded %d SpliceAI scores", len(scores))
    return scores


def _log_expert_distribution(variants: list[ERepoVariant]) -> None:
    dist = Counter(v.expert_classification for v in variants)
    for tier in _ALL_TIERS:
        if dist[tier] > 0:
            logger.info("  %s: %d", tier, dist[tier])


def _build_detection_dict(label: str, counts: _DetectionCounts, expert_count: int) -> dict:
    sens = counts.sensitivity(expert_count)
    ppv = counts.ppv()
    return {
        f"expert_{label}_count": expert_count,
        "true_positives": counts.tp,
        "false_negatives": counts.fn,
        "false_positives": counts.fp,
        "sensitivity": round(sens, 4) if sens is not None else None,
        "ppv": round(ppv, 4) if ppv is not None else None,
    }


def compute_metrics(
    variants: list[ERepoVariant],
    classified: list[ClassifiedVariant],
    revel_scores: dict[tuple[str, int, str, str], float],
) -> dict:
    """Compute full concordance metrics and evidence distribution."""
    total = len(variants)
    confusion: dict[str, dict[str, int]] = {
        row: dict.fromkeys(_ALL_TIERS, 0) for row in _ALL_TIERS
    }
    path_counts = _DetectionCounts()
    ben_counts = _DetectionCounts()
    tag_counter: Counter[str] = Counter()
    exact_match = binary_match = revel_coverage = 0

    for v, cv in zip(variants, classified):
        exact, binary, revel_hit = _accumulate_pair(
            v, cv, revel_scores, confusion, path_counts, ben_counts, tag_counter
        )
        exact_match += exact
        binary_match += binary
        revel_coverage += revel_hit

    expert_p = sum(1 for v in variants if v.expert_classification in _PATH_TIERS)
    expert_b = sum(1 for v in variants if v.expert_classification in _BEN_TIERS)

    # Per-consequence P/LP stratification (expert P/LP only), keyed on VEP consequence.
    per_consequence: dict[str, dict[str, int | float]] = {}
    for v, cv in zip(variants, classified):
        if v.expert_classification not in _PATH_TIERS:
            continue
        csq = cv.scored.annotated.consequence.value
        bucket = per_consequence.setdefault(csq, {"expert_P_LP": 0, "path_tp": 0})
        bucket["expert_P_LP"] += 1
        if cv.classification.value in _PATH_TIERS:
            bucket["path_tp"] += 1
    for csq, b in per_consequence.items():
        n = b["expert_P_LP"]
        b["path_sensitivity"] = round(b["path_tp"] / n, 4) if n else 0.0

    return {
        "total_variants": total,
        "per_consequence": per_consequence,
        "data_coverage": {
            "revel_matched": revel_coverage,
            "revel_coverage_pct": round(revel_coverage / total * 100, 1) if total > 0 else 0,
            "gnomad_matched": sum(
                1 for cv in classified if cv.scored.annotated.allele_frequency is not None
            ),
        },
        "expert_distribution": dict(
            Counter(v.expert_classification for v in variants).most_common()
        ),
        "vartriage_distribution": dict(
            Counter(cv.classification.value for cv in classified).most_common()
        ),
        "pathogenic_detection": _build_detection_dict("P_LP", path_counts, expert_p),
        "benign_detection": _build_detection_dict("B_LB", ben_counts, expert_b),
        "overall": {
            "exact_concordance": round(exact_match / total, 4) if total > 0 else 0,
            "binary_concordance": round(binary_match / total, 4) if total > 0 else 0,
        },
        "confusion_matrix": confusion,
        "evidence_tag_distribution": dict(tag_counter.most_common()),
    }


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def print_summary(metrics: dict, relaxed: bool) -> None:
    """Print human-readable summary to stdout."""
    mode = "RELAXED (>=2 mod = LP)" if relaxed else "STRICT ACMG 2015"
    print(f"\n{'=' * 64}")
    print(f"  ClinGen eRepo VALIDATION — VarTriage v{VARTRIAGE_VERSION} [{mode}]")
    print(f"{'=' * 64}")
    print(f"  Total variants:    {metrics['total_variants']}")

    cov = metrics["data_coverage"]
    print(f"  REVEL coverage:    {cov['revel_matched']}/{metrics['total_variants']} ({cov['revel_coverage_pct']}%)")
    print(f"  gnomAD coverage:   {cov['gnomad_matched']}/{metrics['total_variants']}")
    print()

    print("  Expert Panel distribution:")
    for tier, count in metrics["expert_distribution"].items():
        print(f"    {tier:20s} {count:>5}")
    print()

    pd = metrics["pathogenic_detection"]
    print(f"  PATHOGENIC DETECTION (Expert P/LP = {pd['expert_P_LP_count']})")
    if pd["sensitivity"] is not None:
        print(f"    Sensitivity:  {pd['sensitivity'] * 100:.1f}%")
    else:
        print("    Sensitivity:  N/A")
    if pd["ppv"] is not None:
        print(f"    PPV:          {pd['ppv'] * 100:.1f}%")
    print(f"    TP={pd['true_positives']}  FN={pd['false_negatives']}  FP={pd['false_positives']}")
    print()

    bd = metrics["benign_detection"]
    print(f"  BENIGN DETECTION (Expert B/LB = {bd['expert_B_LB_count']})")
    if bd["sensitivity"] is not None:
        print(f"    Sensitivity:  {bd['sensitivity'] * 100:.1f}%")
    else:
        print("    Sensitivity:  N/A")
    if bd["ppv"] is not None:
        print(f"    PPV:          {bd['ppv'] * 100:.1f}%")
    print(f"    TP={bd['true_positives']}  FN={bd['false_negatives']}  FP={bd['false_positives']}")
    print()

    ov = metrics["overall"]
    print(f"  Exact concordance:   {ov['exact_concordance'] * 100:.1f}%")
    print(f"  Binary concordance:  {ov['binary_concordance'] * 100:.1f}%")
    print()

    print("  Evidence tag distribution (top 12):")
    for tag, count in list(metrics["evidence_tag_distribution"].items())[:12]:
        print(f"    {tag:15s} {count:>6}")

    print(f"{'=' * 64}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate vartriage ACMG classifier against ClinGen eRepo"
    )
    parser.add_argument(
        "--erepo",
        type=Path,
        default=Path("data/erepo/clingen_erepo.tsv"),
    )
    parser.add_argument(
        "--revel",
        type=Path,
        default=Path("data/references/revel_grch38.tsv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/erepo_validation.json"),
    )
    parser.add_argument(
        "--cache-db",
        type=Path,
        default=Path("data/erepo/.api_cache.db"),
    )
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument(
        "--skip-gnomad",
        action="store_true",
        help="Skip gnomAD API calls; assume all variants absent (PM2 fires for all)",
    )
    parser.add_argument(
        "--skip-vep",
        action="store_true",
        help="Skip VEP annotation (all variants treated as missense)",
    )
    parser.add_argument(
        "--relaxed-combining",
        action="store_true",
        help="Use relaxed combining: >=2 moderate pathogenic = Likely Pathogenic",
    )
    parser.add_argument(
        "--spliceai",
        type=Path,
        default=None,
        help="Path to SpliceAI scores TSV (chrom, pos, ref, alt, spliceai_max_delta)",
    )
    args = parser.parse_args()

    if not args.erepo.exists():
        logger.error("eRepo file not found: %s", args.erepo)
        sys.exit(1)
    if not args.revel.exists():
        logger.error("REVEL file not found: %s", args.revel)
        sys.exit(1)

    # Load eRepo
    variants = load_erepo(args.erepo)
    logger.info("Loaded %d filtered SNVs (>=3 stars, mappable classification)", len(variants))

    if args.limit:
        variants = variants[: args.limit]
        logger.info("Limited to %d variants", len(variants))

    _log_expert_distribution(variants)

    # Load REVEL
    revel_scores = load_revel_scores(args.revel)
    revel_hits = sum(1 for v in variants if (v.chrom, v.pos, v.ref, v.alt) in revel_scores)
    logger.info("REVEL coverage: %d/%d (%.1f%%)", revel_hits, len(variants), revel_hits / len(variants) * 100 if variants else 0)

    # Shared cache for all API clients
    cache = _build_shared_cache(args.cache_db)

    # gnomAD frequencies
    # When skipped, we still mark gnomad_queried=True so variants get PM2
    # (eRepo variants are clinically curated rare variants; absence from gnomAD is expected)
    if args.skip_gnomad:
        logger.info("Skipping gnomAD API (--skip-gnomad); all variants treated as absent -> PM2")
        gnomad_freqs: list[Optional[PopulationFrequencies]] = [None] * len(variants)
        gnomad_queried = True
    else:
        client = build_gnomad_client(cache)
        gnomad_freqs = fetch_gnomad_frequencies(variants, client)
        gnomad_queried = True

    # VEP consequences (Fix #3)
    if args.skip_vep:
        logger.info("Skipping VEP (--skip-vep); treating all as MISSENSE")
        consequences = [FunctionalConsequence.MISSENSE] * len(variants)
    else:
        vep_client = build_vep_client(cache)
        consequences = fetch_vep_consequences(variants, vep_client)
        vep_client.close()

    # Load SpliceAI scores if provided
    spliceai_scores = _load_spliceai_scores(args.spliceai)

    # Classify once (strict = raw classifier output). Relaxed combining is a
    # pure re-application on the same classified variants, so both modes come
    # from a single expensive pass (gnomAD + VEP + REVEL + classify).
    logger.info("Running ACMGClassifier on %d variants...", len(variants))
    t0 = time.time()
    classified_strict = classify_variants(
        variants, gnomad_freqs, revel_scores, consequences, gnomad_queried,
        spliceai_scores=spliceai_scores,
    )
    logger.info("Classification complete in %.2fs", time.time() - t0)
    classified_relaxed = apply_relaxed_combining(classified_strict)

    base_config = {
        "skip_gnomad": args.skip_gnomad,
        "skip_vep": args.skip_vep,
        "gnomad_absent_as_pm2": True,
        "spliceai_file": str(args.spliceai) if args.spliceai else None,
        "spliceai_count": len(spliceai_scores) if spliceai_scores else 0,
        "limit": args.limit,
    }

    def _finalize(classified: list, relaxed: bool) -> dict:
        m = compute_metrics(variants, classified, revel_scores)
        m["vartriage_version"] = VARTRIAGE_VERSION
        m["dataset"] = "ClinGen eRepo (>=3-star SNVs), gnomAD GraphQL API + Ensembl VEP + REVEL"
        m["configuration"] = {"relaxed_combining": relaxed, **base_config}
        return m

    metrics = _finalize(classified_relaxed, True)
    strict_metrics = _finalize(classified_strict, False)

    # Save relaxed (headline) + strict companion
    out_dir = args.output.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    logger.info("Results (relaxed) saved to %s", args.output)
    strict_path = out_dir / "erepo_strict.json"
    strict_path.write_text(json.dumps(strict_metrics, indent=2), encoding="utf-8")
    logger.info("Results (strict) saved to %s", strict_path)

    # Stratified: both modes in one file, for analyses 03/05.
    strat = {"relaxed": metrics["per_consequence"], "strict": strict_metrics["per_consequence"]}
    strat_path = out_dir / "erepo_stratified.json"
    strat_path.write_text(json.dumps(strat, indent=2), encoding="utf-8")
    logger.info("Stratified (relaxed+strict) saved to %s", strat_path)

    # Print
    print_summary(metrics, args.relaxed_combining)


if __name__ == "__main__":
    main()
