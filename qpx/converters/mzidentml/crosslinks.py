"""mzIdentML crosslink groups and donor/acceptor positions."""

from __future__ import annotations

from typing import Optional

from qpx.core.cv_terms import CV_CROSSLINK_SII, CV_LOOPLINK_SII, CV_NONCOVALENT_SII


def _build_looplink(
    alpha_pep: dict,
    linker_info: dict,
) -> Optional[list[dict]]:
    """Build cross_links for a looplink (intra-peptide) PSM."""
    xl_mods = alpha_pep.get("xl_mods", {})
    if not xl_mods:
        return None

    links = []
    for xl_value, mods in xl_mods.items():
        donor = next((m for m in mods if m["role"] == "donor"), None)
        acceptor = next((m for m in mods if m["role"] == "acceptor"), None)
        if not donor or not acceptor:
            continue

        linker = linker_info.get(xl_value, linker_info.get("_default", {}))
        links.append(
            {
                "xl_type": "intra",
                "partner_sequence": None,
                "partner_peptidoform": None,
                "donor_position": donor["position"],
                "acceptor_position": acceptor["position"],
                "linker_name": linker.get("name", "unknown"),
                "linker_accession": linker.get("accession"),
                "linker_mass": donor["mass_delta"],
            }
        )

    return links or None


def _interlink_sites(a_mods, b_mods) -> tuple:
    """Resolve donor and acceptor in either peptide orientation."""
    donor = next((m for m in a_mods if m["role"] == "donor"), None)
    acceptor = next((m for m in b_mods if m["role"] == "acceptor"), None)

    # If donor not in alpha, check if donor is in beta (swapped orientation)
    if not donor:
        donor = next((m for m in b_mods if m["role"] == "donor"), None)
        acceptor = next((m for m in a_mods if m["role"] == "acceptor"), None)

    return donor, acceptor


def _build_interlink(
    alpha_pep: dict,
    beta_sii: Optional[dict],
    peptides: dict,
    linker_info: dict,
    xl_type: str,
) -> Optional[list[dict]]:
    """Build cross_links for an inter-peptide or noncovalent PSM."""
    if beta_sii is None:
        return None

    beta_pep = peptides.get(beta_sii["peptide_ref"], {})
    beta_sequence = beta_pep.get("sequence", "")

    # Find donor/acceptor positions from xl_mods
    alpha_xl = alpha_pep.get("xl_mods", {})
    beta_xl = beta_pep.get("xl_mods", {})

    # Match by shared xl_value — handle both orientations:
    # alpha=donor/beta=acceptor OR alpha=acceptor/beta=donor
    links = []
    all_xl_values = set(alpha_xl.keys()) | set(beta_xl.keys())
    for xl_value in all_xl_values:
        a_mods = alpha_xl.get(xl_value, [])
        b_mods = beta_xl.get(xl_value, [])

        donor, acceptor = _interlink_sites(a_mods, b_mods)

        if not donor:
            continue

        linker = linker_info.get(xl_value, linker_info.get("_default", {}))

        links.append(
            {
                "xl_type": xl_type,
                "partner_sequence": beta_sequence,
                "partner_peptidoform": beta_sequence,
                "donor_position": donor["position"],
                "acceptor_position": acceptor["position"] if acceptor else None,
                "linker_name": linker.get("name", "unknown"),
                "linker_accession": linker.get("accession"),
                "linker_mass": donor["mass_delta"],
            }
        )

    # If no xl_mods found (e.g., noncovalent without mods), create a minimal link
    if not links and xl_type == "noncovalent":
        links.append(
            {
                "xl_type": "noncovalent",
                "partner_sequence": beta_sequence,
                "partner_peptidoform": beta_sequence,
                "donor_position": 0,
                "acceptor_position": None,
                "linker_name": "noncovalent",
                "linker_accession": None,
                "linker_mass": 0.0,
            }
        )

    return links or None


def _sii_link_group(sii) -> tuple:
    """Read crosslink group and association flags from one SII."""
    xl_group = None
    is_looplink = False
    is_noncovalent = False

    for cv in sii["cv_params"]:
        if cv["accession"] == CV_CROSSLINK_SII:
            xl_group = cv["value"]
        elif cv["accession"] == CV_LOOPLINK_SII:
            is_looplink = True
        elif cv["accession"] == CV_NONCOVALENT_SII:
            is_noncovalent = True
            xl_group = cv["value"]

    return xl_group, is_looplink, is_noncovalent


def _paired_sii_groups(paired_groups, xl_type) -> list[dict]:
    """Pair the first two SIIs, retaining orphan items as regular PSMs."""
    groups = []
    for _group_val, paired in paired_groups.items():
        if len(paired) >= 2:
            # First SII is alpha (donor), second is beta (acceptor)
            groups.append({"xl_type": xl_type, "alpha": paired[0], "beta": paired[1]})
        else:
            # Orphan SII — treat as regular
            groups.append({"xl_type": "none", "alpha": paired[0]})
    return groups


def _group_siis(siis: list[dict]) -> list[dict]:
    """Group SII elements into classified cross-link groups.

    Returns list of dicts: ``{xl_type, alpha, beta?}``
    where ``xl_type`` is one of: ``"none"``, ``"inter"``, ``"intra"``, ``"noncovalent"``.
    """
    groups = []

    # Index SIIs by their crosslink group value
    xl_groups: dict[str, list[dict]] = {}  # group_value → [sii, ...]
    looplink_siis = []
    noncovalent_groups: dict[str, list[dict]] = {}
    regular_siis = []

    for sii in siis:
        xl_group, is_looplink, is_noncovalent = _sii_link_group(sii)

        if is_looplink:
            looplink_siis.append(sii)
        elif is_noncovalent and xl_group:
            noncovalent_groups.setdefault(xl_group, []).append(sii)
        elif xl_group:
            xl_groups.setdefault(xl_group, []).append(sii)
        else:
            regular_siis.append(sii)

    # Regular PSMs
    for sii in regular_siis:
        groups.append({"xl_type": "none", "alpha": sii})

    # Looplinks (intra-peptide)
    for sii in looplink_siis:
        groups.append({"xl_type": "intra", "alpha": sii})

    # Inter-peptide crosslinks (paired SIIs)
    groups.extend(_paired_sii_groups(xl_groups, "inter"))
    # Noncovalent associations
    groups.extend(_paired_sii_groups(noncovalent_groups, "noncovalent"))

    return groups
