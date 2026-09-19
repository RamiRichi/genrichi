"""
GenRichi Phase 5.2/5.3 — inspection of reference / annotation resources.

Pure standard-library helpers that read *observed facts* from the resource
files themselves (FASTA index, sequence dictionary, BWA index, VCF headers,
VEP cache info.txt, index magic bytes) rather than trusting configured paths
or filenames. Used by preflight.py (fail-fast validation) and by other
pipeline helpers that need the same observations.

Nothing here changes, thresholds or reinterprets any scientific value.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
from datetime import datetime, timezone

# Files larger than this are not hashed on every run (dbSNP, gnomAD and the
# genome FASTA are multi-GB); their size and mtime are recorded instead.
SHA256_MAX_BYTES = 512 * 1024 * 1024

BWA_INDEX_SUFFIXES = (".amb", ".ann", ".bwt", ".pac", ".sa")

# An index older than its data by less than this is normal (files are published
# / copied minutes apart); only a larger gap is reported as possibly stale.
STALE_TOLERANCE_SECONDS = 24 * 3600

_BUILD_ALIASES = {"grch38": "GRCh38", "hg38": "GRCh38", "grch37": "GRCh37", "hg19": "GRCh37"}

# Keys of the VEP cache info.txt that describe what the cache contains.
_CACHE_INFO_KEYS = (
    "species", "assembly", "source_assembly", "source_gencode", "source_genebuild",
    "source_ClinVar", "source_COSMIC", "source_dbSNP", "source_gnomADe", "source_gnomADg",
    "source_1000genomes", "source_HGMD-PUBLIC", "source_sift", "source_polyphen",
    "source_regbuild", "sift", "polyphen", "regulatory",
)


class ResourceError(ValueError):
    """A resource file is present but malformed / unreadable."""


# ── Generic file facts ────────────────────────────────────────────────────
def _utc_iso(ts: float) -> str:
    return datetime.fromtimestamp(ts, tz=timezone.utc).isoformat(timespec="seconds")


def sha256_file(path, chunk_size: int = 1 << 20):
    """SHA-256 of a file, or None if missing/unreadable."""
    if not path:
        return None
    try:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(chunk_size), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError:
        return None


def file_facts(path, sha256_max_bytes: int = SHA256_MAX_BYTES) -> dict:
    """Observed size / mtime / SHA-256 of a file (SHA-256 skipped, with a reason, for very large files)."""
    facts = {"path": str(path) if path else None, "exists": False, "size_bytes": None,
             "mtime_utc": None, "sha256": None, "sha256_skipped_reason": None}
    if not path:
        return facts
    try:
        st = os.stat(path)
    except OSError:
        return facts
    if not os.path.isfile(path):
        return facts
    facts.update(exists=True, size_bytes=st.st_size, mtime_utc=_utc_iso(st.st_mtime))
    if st.st_size > sha256_max_bytes:
        facts["sha256_skipped_reason"] = f"size > {sha256_max_bytes // (1024 * 1024)} MiB"
    else:
        facts["sha256"] = sha256_file(path)
    return facts


def vep_cache_version_from_extra(vep_extra):
    """`--cache_version N` from a VEP extra-args string, or None if not pinned."""
    if not vep_extra:
        return None
    m = re.search(r"--cache_version[= ](\S+)", str(vep_extra))
    return m.group(1) if m else None


def normalize_build(build) -> str | None:
    if build is None:
        return None
    return _BUILD_ALIASES.get(str(build).strip().lower())


def build_family_from_text(text) -> str | None:
    """'GRCh38' / 'GRCh37' if `text` unambiguously names one build, else None."""
    t = str(text or "").lower()
    has38 = "grch38" in t or "hg38" in t
    has37 = "grch37" in t or "hg19" in t or "b37" in t
    if has38 and not has37:
        return "GRCh38"
    if has37 and not has38:
        return "GRCh37"
    return None


# ── FASTA / .fai / .dict ─────────────────────────────────────────────────
def read_fai(path) -> list:
    """[(name, length, offset, linebases, linewidth), ...]; raises ResourceError if malformed."""
    rows = []
    try:
        with open(path, encoding="utf-8") as fh:
            for lineno, line in enumerate(fh, 1):
                line = line.rstrip("\n\r")
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 5:
                    raise ResourceError(f"{path}: line {lineno} has {len(parts)} column(s), expected >= 5")
                try:
                    rows.append((parts[0], int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])))
                except ValueError:
                    raise ResourceError(f"{path}: line {lineno} has non-numeric length/offset fields")
    except OSError as exc:
        raise ResourceError(f"{path}: cannot read ({exc})")
    if not rows:
        raise ResourceError(f"{path}: contains no sequences")
    return rows


def read_dict(path) -> list:
    """[(name, length, m5_or_None), ...] from the @SQ lines; raises ResourceError if malformed/empty."""
    rows = []
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                if not line.startswith("@SQ"):
                    continue
                fields = dict(f.split(":", 1) for f in line.rstrip("\n\r").split("\t")[1:] if ":" in f)
                if "SN" not in fields or "LN" not in fields:
                    raise ResourceError(f"{path}: @SQ line without SN/LN: {line.strip()[:80]}")
                try:
                    rows.append((fields["SN"], int(fields["LN"]), fields.get("M5")))
                except ValueError:
                    raise ResourceError(f"{path}: @SQ line with non-numeric LN: {line.strip()[:80]}")
    except OSError as exc:
        raise ResourceError(f"{path}: cannot read ({exc})")
    if not rows:
        raise ResourceError(f"{path}: contains no @SQ sequence records")
    return rows


def dict_path_for(genome_path) -> str:
    return os.path.splitext(str(genome_path))[0] + ".dict"


def check_fasta_indexes(genome_path):
    """
    Verify FASTA, .fai and .dict are mutually consistent and that the .fai
    describes this FASTA (extents inside the file, no unindexed trailing data).

    Returns (errors, warnings, facts) — lists of message strings + observed facts.
    Missing files are NOT reported here (validate_resources already does).
    """
    errors: list = []
    warnings: list = []
    facts: dict = {}
    genome_path = str(genome_path)
    fai_path = genome_path + ".fai"
    dpath = dict_path_for(genome_path)
    if not (os.path.isfile(genome_path) and os.path.isfile(fai_path)):
        return errors, warnings, facts

    try:
        fai = read_fai(fai_path)
    except ResourceError as exc:
        return [f"Unusable FASTA index: {exc}"], warnings, facts

    facts.update(
        contig_count=len(fai),
        total_bases=sum(r[1] for r in fai),
        first_contig=fai[0][0],
        fasta_size_bytes=os.path.getsize(genome_path),
        fai_sha256=sha256_file(fai_path),
    )
    names = [r[0] for r in fai]
    if len(set(names)) != len(names):
        errors.append(f"{fai_path}: duplicate contig names in the FASTA index")

    # FASTA <-> fai
    try:
        with open(genome_path, "rb") as fh:
            head = fh.read(4096)
    except OSError as exc:
        errors.append(f"Cannot read FASTA {genome_path}: {exc}")
        head = b""
    if head:
        if not head.startswith(b">"):
            errors.append(f"{genome_path}: does not look like a FASTA file (first byte is not '>')")
        else:
            tokens = head.split(b"\n", 1)[0][1:].split()
            first_name = tokens[0].decode("utf-8", "replace") if tokens else ""
            if first_name != fai[0][0]:
                errors.append(
                    f"FASTA first contig '{first_name}' does not match .fai first contig '{fai[0][0]}' "
                    "(index belongs to a different FASTA?)"
                )
    max_end = 0
    for name, length, offset, linebases, linewidth in fai:
        if linebases <= 0:
            errors.append(f"{fai_path}: contig '{name}' has invalid line length {linebases}")
            break
        end = offset + length + ((length + linebases - 1) // linebases - 1) * (linewidth - linebases)
        max_end = max(max_end, end)
    size = facts["fasta_size_bytes"]
    if max_end > size:
        errors.append(
            f"The .fai describes sequence beyond the end of the FASTA (index end byte {max_end} > file size {size}): "
            "FASTA truncated or .fai belongs to a different file"
        )
    elif size - max_end > 64:
        errors.append(
            f"FASTA has {size - max_end} bytes after the last contig indexed by the .fai: "
            "the .fai is stale (FASTA contains contigs the index does not list)"
        )

    # dict <-> fai
    facts["dict_path"] = dpath
    if os.path.isfile(dpath):
        facts["dict_sha256"] = sha256_file(dpath)
        try:
            dic = read_dict(dpath)
        except ResourceError as exc:
            errors.append(f"Unusable sequence dictionary: {exc}")
            dic = None
        if dic is not None:
            facts["dict_contig_count"] = len(dic)
            fai_pairs = [(r[0], r[1]) for r in fai]
            dic_pairs = [(r[0], r[1]) for r in dic]
            if fai_pairs != dic_pairs:
                if len(fai_pairs) != len(dic_pairs):
                    errors.append(
                        f"Sequence dictionary and .fai disagree: {len(dic_pairs)} dictionary contigs vs "
                        f"{len(fai_pairs)} .fai contigs"
                    )
                else:
                    idx = next(i for i, (a, b) in enumerate(zip(fai_pairs, dic_pairs)) if a != b)
                    errors.append(
                        f"Sequence dictionary and .fai disagree at contig #{idx + 1}: "
                        f".fai has {fai_pairs[idx]}, .dict has {dic_pairs[idx]}"
                    )
        if os.path.getmtime(dpath) < os.path.getmtime(genome_path) - STALE_TOLERANCE_SECONDS:
            warnings.append(f"Sequence dictionary is older than the FASTA (possibly stale): {dpath}")
    if os.path.getmtime(fai_path) < os.path.getmtime(genome_path) - STALE_TOLERANCE_SECONDS:
        warnings.append(f"FASTA index is older than the FASTA (possibly stale): {fai_path}")
    return errors, warnings, facts


def check_bwa_index(genome_path, expected_total_bases=None, expected_contigs=None):
    """
    Verify the classic BWA index (.amb/.ann/.bwt/.pac/.sa — required by `bwa mem`)
    exists, is non-empty, and was built from a FASTA of the same size:
    .ann first line is 'l_pac n_seqs seed'; .pac has ceil(l_pac/4)+1 bytes.
    Returns (errors, warnings, facts).
    """
    errors: list = []
    warnings: list = []
    facts: dict = {}
    genome_path = str(genome_path)
    missing = [s for s in BWA_INDEX_SUFFIXES if not os.path.isfile(genome_path + s)]
    if missing:
        errors.append(
            f"BWA index incomplete for {genome_path}: missing {', '.join(missing)} "
            "(the pipeline aligns with `bwa mem`, which needs all of .amb/.ann/.bwt/.pac/.sa)"
        )
        return errors, warnings, facts
    empty = [s for s in BWA_INDEX_SUFFIXES if os.path.getsize(genome_path + s) == 0]
    if empty:
        errors.append(f"BWA index file(s) are empty for {genome_path}: {', '.join(empty)}")
        return errors, warnings, facts

    try:
        with open(genome_path + ".ann", encoding="utf-8") as fh:
            l_pac, n_seqs = (int(x) for x in fh.readline().split()[:2])
    except (OSError, ValueError):
        errors.append(f"{genome_path}.ann: cannot parse BWA index header")
        return errors, warnings, facts
    facts.update(bwa_index_l_pac=l_pac, bwa_index_n_seqs=n_seqs)
    if expected_total_bases is not None and l_pac != expected_total_bases:
        errors.append(
            f"BWA index was built from a different FASTA: index covers {l_pac} bases, "
            f"the .fai lists {expected_total_bases}"
        )
    if expected_contigs is not None and n_seqs != expected_contigs:
        errors.append(
            f"BWA index was built from a different FASTA: index has {n_seqs} contigs, "
            f"the .fai lists {expected_contigs}"
        )
    pac_expected = (l_pac + 3) // 4 + 1
    pac_size = os.path.getsize(genome_path + ".pac")
    if pac_size != pac_expected:
        errors.append(
            f"{genome_path}.pac has {pac_size} bytes, expected {pac_expected} for a {l_pac}-base index "
            "(index incomplete or corrupted)"
        )
    fasta_mtime = os.path.getmtime(genome_path)
    stale = [s for s in BWA_INDEX_SUFFIXES if os.path.getmtime(genome_path + s) < fasta_mtime - STALE_TOLERANCE_SECONDS]
    if stale:
        warnings.append(f"BWA index file(s) older than the FASTA (possibly stale): {', '.join(stale)}")
    return errors, warnings, facts


# ── Reference manifest ("expected resource set") ─────────────────────────
_MANIFEST_COMPARE_KEYS = ("contig_count", "total_bases", "fai_sha256", "dict_sha256")


def build_reference_manifest(genome_path, genome_build, description=None) -> dict:
    """Fingerprint of the current reference (contig set, .fai and .dict hashes) for later comparison."""
    errors, _warnings, facts = check_fasta_indexes(genome_path)
    if errors or not facts:
        raise ResourceError("Cannot build a manifest from an inconsistent reference: " + "; ".join(errors or ["missing index"]))
    return {
        "schema_version": "1.0",
        "description": description,
        "genome_build": normalize_build(genome_build) or genome_build,
        "reference": {
            "fasta_basename": os.path.basename(str(genome_path)),
            "contig_count": facts["contig_count"],
            "total_bases": facts["total_bases"],
            "first_contig": facts["first_contig"],
            "fai_sha256": facts["fai_sha256"],
            "dict_sha256": facts.get("dict_sha256"),
        },
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }


def compare_to_manifest(manifest: dict, facts: dict, genome_build) -> list:
    """Error strings for every way the observed reference differs from the manifest."""
    errors = []
    ref = (manifest or {}).get("reference") or {}
    exp_build = (manifest or {}).get("genome_build")
    if exp_build and normalize_build(genome_build) and normalize_build(exp_build) != normalize_build(genome_build):
        errors.append(f"manifest expects genome build {exp_build}, config declares {genome_build}")
    for key in _MANIFEST_COMPARE_KEYS:
        expected = ref.get(key)
        if expected is not None and facts.get(key) != expected:
            errors.append(f"{key}: expected {expected}, observed {facts.get(key)}")
    return errors


# ── VCF resources ────────────────────────────────────────────────────────
def is_bgzf(path) -> bool:
    try:
        with open(path, "rb") as fh:
            head = fh.read(18)
    except OSError:
        return False
    return len(head) >= 14 and head[:4] == b"\x1f\x8b\x08\x04" and head[12:14] == b"BC"


def _index_candidates(vcf_path):
    return [(vcf_path + ".tbi", b"TBI\x01"), (vcf_path + ".csi", b"CSI\x01")]


def check_vcf_index(vcf_path):
    """
    A bgzip-compressed VCF with a *usable* index: bgzf framing, index has the
    correct magic bytes, index not older than the data.
    Returns (errors, warnings). Presence alone is reported elsewhere.
    """
    errors: list = []
    warnings: list = []
    vcf_path = str(vcf_path)
    if not os.path.isfile(vcf_path):
        return errors, warnings
    if vcf_path.endswith(".gz") and not is_bgzf(vcf_path):
        errors.append(f"{vcf_path}: not BGZF-compressed (plain gzip cannot be tabix-indexed / randomly accessed)")
    found = [(p, m) for p, m in _index_candidates(vcf_path) if os.path.isfile(p)]
    if not found:
        return errors, warnings
    idx_path, magic = found[0]
    try:
        with gzip.open(idx_path, "rb") as fh:
            got = fh.read(4)
    except (OSError, EOFError):
        got = b""
    if got != magic:
        errors.append(f"{idx_path}: not a valid {magic[:3].decode()} index (bad or unreadable header)")
    elif os.path.getmtime(idx_path) < os.path.getmtime(vcf_path) - STALE_TOLERANCE_SECONDS:
        warnings.append(f"{idx_path} is older than {os.path.basename(vcf_path)} (index possibly stale)")
    return errors, warnings


def read_vcf_header(path, max_lines: int = 20000) -> dict:
    """
    Observed VCF meta-information: fileformat, fileDate, source, reference,
    dbSNP_BUILD_ID, contig header count. Empty dict + 'error' on failure.
    """
    out: dict = {"fileformat": None, "file_date": None, "source": None, "reference": None,
                 "dbsnp_build_id": None, "contig_header_count": 0}
    keymap = {"fileformat": "fileformat", "fileDate": "file_date", "source": "source",
              "reference": "reference", "dbSNP_BUILD_ID": "dbsnp_build_id"}
    try:
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            for n, line in enumerate(fh):
                if not line.startswith("##"):
                    break
                if n >= max_lines:
                    break
                if line.startswith("##contig="):
                    out["contig_header_count"] += 1
                    continue
                m = re.match(r"##([A-Za-z_]+)=(.*)$", line.rstrip("\n\r"))
                if m and m.group(1) in keymap and out[keymap[m.group(1)]] is None:
                    out[keymap[m.group(1)]] = m.group(2)
    except (OSError, EOFError) as exc:
        out["error"] = f"cannot read VCF header: {exc}"
    return out


# ── VEP cache ────────────────────────────────────────────────────────────
def read_cache_info_txt(path) -> dict:
    info = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n\r").split("\t", 1)
            if len(parts) == 2 and parts[0] in _CACHE_INFO_KEYS:
                info[parts[0]] = parts[1]
    return info


def inspect_vep_cache(cache_dir, cache_version, genome_build, species: str = "homo_sapiens"):
    """
    Observe the VEP cache itself: which `<version>_<assembly>` directory
    exists and what its info.txt says about assembly and sources.

    Returns (errors, warnings, facts). `cache_version` may be None (unpinned).
    """
    errors: list = []
    warnings: list = []
    facts: dict = {"cache_dir": str(cache_dir) if cache_dir else None}
    expected_asm = normalize_build(genome_build)
    if not cache_dir or not os.path.isdir(cache_dir):
        return errors, warnings, facts  # reported by validate_resources
    if expected_asm is None:
        return errors, warnings, facts  # unknown build is reported by validate_config_ranges

    species_dir = os.path.join(cache_dir, species)
    candidates = []
    if os.path.isdir(species_dir):
        candidates = sorted(d for d in os.listdir(species_dir)
                            if re.fullmatch(rf"\d+_{re.escape(expected_asm)}", d)
                            and os.path.isdir(os.path.join(species_dir, d)))
    if cache_version:
        version_dir = f"{cache_version}_{expected_asm}"
        if version_dir not in candidates:
            errors.append(
                f"VEP cache directory not found: {os.path.join(species_dir, version_dir)} "
                f"(available for {expected_asm}: {', '.join(candidates) or 'none'}) — "
                f"--cache_version {cache_version} / genome_build {genome_build} do not match the installed cache"
            )
            return errors, warnings, facts
    elif len(candidates) == 1:
        version_dir = candidates[0]
        warnings.append(
            "VEP --cache_version is not pinned in annotation.vep.extra; the only installed "
            f"{expected_asm} cache is {version_dir}; pin it with --cache_version to make the "
            "executable/cache pairing explicit"
        )
    elif not candidates:
        errors.append(f"No VEP cache for {expected_asm} found under {species_dir}")
        return errors, warnings, facts
    else:
        errors.append(
            f"VEP --cache_version is not pinned and several {expected_asm} caches are installed "
            f"({', '.join(candidates)}); pin one with --cache_version in annotation.vep.extra"
        )
        return errors, warnings, facts

    cache_path = os.path.join(species_dir, version_dir)
    facts.update(cache_path=cache_path, version_dir=version_dir, cache_version=version_dir.split("_", 1)[0])
    info_path = os.path.join(cache_path, "info.txt")
    if not os.path.isfile(info_path):
        errors.append(f"VEP cache {cache_path} has no info.txt — cache is incomplete, its assembly cannot be verified")
        return errors, warnings, facts
    try:
        info = read_cache_info_txt(info_path)
    except OSError as exc:
        errors.append(f"Cannot read {info_path}: {exc}")
        return errors, warnings, facts
    facts["info_txt_sha256"] = sha256_file(info_path)
    facts["info"] = info
    observed_asm = normalize_build(info.get("assembly"))
    facts["assembly"] = info.get("assembly")
    if observed_asm is None:
        errors.append(f"{info_path}: missing or unrecognised 'assembly' entry ({info.get('assembly')!r})")
    elif observed_asm != expected_asm:
        errors.append(
            f"VEP cache assembly mismatch: info.txt says {info.get('assembly')}, "
            f"config declares {genome_build} ({cache_path})"
        )
    return errors, warnings, facts
