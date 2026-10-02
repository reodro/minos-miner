"""Turn model scores into an answer set. Owner: answer_set."""


def assemble_records(keys, p_any, p_hom, thr, mode_giab, donor, allowed=None):
    """Return records (pos, ref, alt, zyg, synthetic, confidence).

    GIAB-person round: the donor's official alleles are the background, exactly as published
    (including two-allele 1/2 sites); model calls add injections at positions the donor does
    not already occupy. Synthetic round: everything comes from the model.
    Model calls keep at most two alleles per position (a person has two copies): the most
    confident ones. Two alleles at one position means genotype 1/2, so both are het.
    """
    recs = [(k[0], k[1], k[2], z, False, 2.0) for k, z in donor.items()]
    taken = {k[0] for k in donor}
    by_pos = {}
    for k, pa, ph in zip(keys, p_any, p_hom):
        if pa < thr or k in donor or k[0] in taken:
            continue
        if mode_giab and allowed is not None and not allowed(k):
            continue          # Minos never inserts outside the person's confident region
        by_pos.setdefault(k[0], []).append((k[0], k[1], k[2], "hom" if ph >= 0.5 else "het",
                                            bool(mode_giab), float(pa)))
    for alleles in by_pos.values():
        alleles.sort(key=lambda r: -r[5])
        top = alleles[:2]
        if len(top) == 2:
            top = [r[:3] + ("het",) + r[4:] for r in top]
        recs.extend(top)
    return recs
