"""
Shared builders for tiny but structurally real resource fixtures used by the
tests: BGZF-compressed VCFs with tabix-magic indexes, a consistent
FASTA + .fai + .dict + classic BWA index, and a VEP cache directory with an
info.txt.
"""

import gzip
import struct
import zlib
from pathlib import Path

BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def bgzf_bytes(text: str) -> bytes:
    """A valid BGZF file (one data block + the standard EOF block) holding `text`."""
    data = text.encode("utf-8")
    comp = zlib.compressobj(6, zlib.DEFLATED, -15)
    payload = comp.compress(data) + comp.flush()
    header = (b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
              + struct.pack("<H", 18 + len(payload) + 8 - 1))
    trailer = struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data))
    return header + payload + trailer + BGZF_EOF


def write_vcf_with_index(path: Path, header_lines):
    """A BGZF-compressed VCF (header only) plus a tabix-magic .tbi."""
    body = "\n".join(header_lines) + "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    path.write_bytes(bgzf_bytes(body))
    Path(str(path) + ".tbi").write_bytes(gzip.compress(b"TBI\x01" + b"\x00" * 32))


def write_reference(tmp: Path, name="hg38.fa"):
    """Tiny but fully consistent reference: FASTA + .fai + .dict + classic BWA index."""
    genome = tmp / name
    genome.write_bytes(b">chr1\nACGT\n")
    (tmp / (name + ".fai")).write_bytes(b"chr1\t4\t6\t4\t5\n")
    (tmp / (Path(name).stem + ".dict")).write_bytes(b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:4\n")
    (tmp / (name + ".amb")).write_bytes(b"4 1 0\n")
    (tmp / (name + ".ann")).write_bytes(b"4 1 11\n0 chr1 (null)\n0 4 0\n")
    (tmp / (name + ".bwt")).write_bytes(b"bwt-bytes")
    (tmp / (name + ".pac")).write_bytes(b"\x00\x00")  # ceil(4/4)+1 bytes
    (tmp / (name + ".sa")).write_bytes(b"sa-bytes")
    return genome


def write_vep_cache(tmp: Path, version="113", assembly="GRCh38", info_assembly=None):
    d = tmp / "vep_cache" / "homo_sapiens" / f"{version}_{assembly}"
    d.mkdir(parents=True)
    (d / "info.txt").write_text(
        f"species\thomo_sapiens\nassembly\t{info_assembly or assembly}\nsource_gencode\tGENCODE 47\n"
        "source_ClinVar\t202404\nsource_COSMIC\t99\nsource_dbSNP\t156\n", encoding="utf-8")
    return tmp / "vep_cache"
