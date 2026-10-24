"""mzIdentML numeric and modification-localization scores."""

from __future__ import annotations

import re

from qpx.core.cv_terms import CV_PEP, CV_PEP_GLOBAL, SITE_LOCALIZATION_ACCESSIONS, SKIP_SCORE_ACCESSIONS
from qpx.core.scores import is_higher_better, normalize_score_name

_SKIP_SCORE_ACCESSIONS = SKIP_SCORE_ACCESSIONS
_PEP_ACCESSIONS = frozenset({CV_PEP, CV_PEP_GLOBAL})


def _parse_position_scores(value_str: str) -> dict[int, float]:
    """Parse residue-position localization scores from one cvParam value."""
    pos_scores: dict[int, float] = {}
    for token in re.split(r"[,;]\s*", value_str):
        m = re.match(r"([A-Za-z])(\d+)\s*:\s*([\d.eE+-]+)", token.strip())
        if m:
            pos = int(m.group(2))
            try:
                pos_scores[pos] = float(m.group(3))
            except ValueError:
                pass

    return pos_scores


def _attach_position_scores(modifications, pos_scores, score_name) -> None:
    """Attach localization values to every matching modification position."""
    for mod in modifications:
        for pos_info in mod.get("positions", []):
            p = pos_info.get("position", 0)
            if p in pos_scores:
                scores = pos_info.get("scores") or []
                scores.append(
                    {
                        "score_name": score_name,
                        "score_value": pos_scores[p],
                        "higher_better": is_higher_better(score_name),
                    }
                )
                pos_info["scores"] = scores


def _attach_single_score(modifications, value_str, score_name) -> None:
    """Attach a scalar localization value to the first modification position."""
    try:
        val = float(value_str)
    except (ValueError, TypeError):
        return
    for mod in modifications:
        positions = mod.get("positions", [])
        if positions:
            scores = positions[0].get("scores") or []
            scores.append(
                {
                    "score_name": score_name,
                    "score_value": val,
                    "higher_better": is_higher_better(score_name),
                }
            )
            positions[0]["scores"] = scores
            break


def _attach_site_localization_scores(
    modifications: list[dict] | None,
    cv_params: list[dict],
) -> list[dict] | None:
    """Enrich modification positions with site localization scores from SII cvParams.

    In mzIdentML, site localization scores (phosphoRS, ptmRS, Ascore) appear as
    cvParams on the SpectrumIdentificationItem.  Their ``value`` field typically
    encodes per-position probabilities, e.g. ``"S3: 95.3, T5: 2.1, S7: 2.6"``.
    When a single float is present it is applied to the first modification position.
    """
    if not modifications:
        return modifications

    # Collect site-localization cvParams
    site_cvs = [cv for cv in cv_params if cv["accession"] in SITE_LOCALIZATION_ACCESSIONS]
    if not site_cvs:
        return modifications

    for cv in site_cvs:
        score_name = normalize_score_name(cv["name"])
        value_str = cv.get("value", "")
        if not value_str:
            continue

        # Try to parse position-keyed format: "S3: 95.3, T5: 2.1"
        pos_scores = _parse_position_scores(value_str)
        if pos_scores:
            # Attach scores to matching modification positions
            _attach_position_scores(modifications, pos_scores, score_name)
        else:
            # Single float value — apply to first modification position
            _attach_single_score(modifications, value_str, score_name)

    return modifications


def _extract_scores(cv_params: list[dict]) -> list[dict]:
    """Extract numeric scores from SII cvParams."""
    scores = []
    for cv in cv_params:
        acc = cv["accession"]
        if acc in _SKIP_SCORE_ACCESSIONS:
            continue
        value_str = cv.get("value", "")
        if not value_str:
            continue
        # Try to parse as float
        try:
            value = float(value_str)
        except (ValueError, TypeError):
            continue
        normalized = normalize_score_name(cv["name"])
        scores.append(
            {
                "score_name": normalized,
                "score_value": value,
                "higher_better": is_higher_better(normalized),
            }
        )
    return scores


def _extract_pep(cv_params: list[dict]) -> float | None:
    """Extract posterior error probability from SII cvParams.

    Looks for MS:1001493 (PEP) or MS:1002352 (global PEP).
    Returns float value or None.
    """
    for cv in cv_params:
        if cv["accession"] in _PEP_ACCESSIONS:
            try:
                return float(cv.get("value", ""))
            except (ValueError, TypeError):
                pass
    return None
