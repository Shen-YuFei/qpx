"""OpenMS modification and site-localisation fields."""

from __future__ import annotations

import ast


def _parse_site_scores(raw) -> dict[int, float]:
    """Parse a Luciphor ``Luciphor_site_scores`` UserParam into {position: score}.

    The value is a stringified dict (e.g. ``{1: -31.36, 9: 31.36}``) keyed by
    1-based amino-acid index. Positional ``None`` values are dropped so a site
    without a localisation score contributes no score entry.
    """
    if not raw:
        return {}
    try:
        parsed = ast.literal_eval(str(raw))
    except (ValueError, SyntaxError, TypeError):
        return {}  # malformed literal -> no scores
    if not isinstance(parsed, dict):
        return {}
    out = {}
    for k, v in parsed.items():
        try:
            out[int(k)] = float(v)
        except (TypeError, ValueError):
            continue  # unparseable position or a None/non-numeric score
    return out


def to_modifications(aa_sequence, site_scores: dict[int, list[dict]] | None = None) -> list[dict]:
    """Convert a pyopenms ``AASequence`` to the QPX ``modifications`` structure.

    Each modification is ``{name, accession, positions: [{position, amino_acid,
    scores}]}`` — the residue position (1-based) plus, when a site-localisation
    score dict is supplied, the per-site localisation score(s) (as QPX ``score``
    structs). N- and C-terminal modifications are reported at position 1 /
    ``len(sequence)`` with an ``amino_acid`` of ``None`` (no single residue
    carries them), mirroring ``to_proforma``'s ``[UNIMOD:N]-`` / ``-[UNIMOD:N]``
    rendering.
    """
    if not aa_sequence.isModified():
        return None
    mods: list[dict] = []
    seen: dict[tuple[str, str | None], dict] = {}

    def _add_mod(name: str, accession: str | None, position: int, amino_acid: str | None, scores: list[dict] | None) -> None:
        key = (name, accession)
        entry = seen.get(key)
        if entry is None:
            entry = {"name": name, "accession": accession, "positions": []}
            seen[key] = entry
            mods.append(entry)
        entry["positions"].append({"position": position, "amino_acid": amino_acid, "scores": scores})

    def _acc(accession) -> str | None:
        # Normalise pyopenms "UniMod:N" to QPX's canonical "UNIMOD:N".
        if not accession:
            return None
        return str(accession).replace("UniMod:", "UNIMOD:")

    for i in range(aa_sequence.size()):
        mod = aa_sequence[i].getModification()
        if not mod:
            continue
        _add_mod(
            mod.getId() or "",
            _acc(mod.getUniModAccession()),
            i + 1,
            aa_sequence[i].getOneLetterCode(),
            site_scores.get(i + 1) if site_scores else None,
        )

    nterm = aa_sequence.getNTerminalModification()
    if nterm:
        # N-terminal position is always 0 (amino_acid null); see the QPX
        # modifications spec.
        _add_mod(nterm.getId() or "", _acc(nterm.getUniModAccession()), 0, None, None)
    cterm = aa_sequence.getCTerminalModification()
    if cterm:
        # C-terminal position is always length + 1 (amino_acid null).
        _add_mod(cterm.getId() or "", _acc(cterm.getUniModAccession()), aa_sequence.size() + 1, None, None)
    return mods or None


# Site-localisation algorithms whose UserParams OpenMS may stamp onto peptide
# hits: ``<name>_pep_score`` (peptide-level score) and ``<name>_site_scores``
# (a ``{position: score}`` dict keyed by 1-based residue index). Multiple
# algorithms may be present in one file (e.g. a run localised by Ascore and
# another by PhosphoRS), so all matches are collected rather than a single one.
_LOCALIZATION_ALGORITHMS = ("Luciphor", "Ascore", "PhosphoRS", "AScore")


def localization_scores(hit) -> tuple[list[dict] | None, dict[int, list[dict]]]:
    """Return ``(additional_scores, site_scores)`` from a hit's site-localisation UserParams.

    For each known localisation algorithm present on the hit, ``<Algo>_pep_score``
    and ``<Algo>_global_flr`` (when present) become ``additional_scores`` entries
    (the QPX ``score`` struct); ``<Algo>_site_scores`` is parsed into
    ``{position: [score, ...]}`` so :func:`to_modifications` can attach the
    per-site localisation score(s) to the matching modification position.
    """
    additional: list[dict] = []
    site_scores: dict[int, list[dict]] = {}
    for algo in _LOCALIZATION_ALGORITHMS:
        pep_key = f"{algo}_pep_score"
        if hit.metaValueExists(pep_key):
            additional.append(
                {"score_name": pep_key.lower(), "score_value": float(hit.getMetaValue(pep_key)), "higher_better": None}
            )
        flr_key = f"{algo}_global_flr"
        if hit.metaValueExists(flr_key):
            additional.append(
                {"score_name": flr_key.lower(), "score_value": float(hit.getMetaValue(flr_key)), "higher_better": None}
            )
        site_key = f"{algo}_site_scores"
        if hit.metaValueExists(site_key):
            parsed = _parse_site_scores(hit.getMetaValue(site_key))
            for pos, score in parsed.items():
                site_scores.setdefault(pos, []).append(
                    {"score_name": f"{algo.lower()}_site_score", "score_value": score, "higher_better": None}
                )
    return (additional or None), site_scores
