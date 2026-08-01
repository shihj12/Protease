"""Stream and reduce the Nature study's deposited MaxQuant results.

The source is PXD024364 / MSV000086944. The 19.9 GB ZIP is never downloaded:
only the compressed bytes of the two required members are read over FTP using
ZIP offsets. Compact, derived caches are retained locally for later runs.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import os
import struct
import zlib
from dataclasses import dataclass
from ftplib import FTP_TLS
from pathlib import Path
from statistics import median
from typing import Callable, Iterator

from coverage_analysis.protease_coverage import NATURE_ENZYMES, interval_mask


FTP_HOST = "massive-ftp.ucsd.edu"
FTP_ROOT = "/v03/MSV000086944"
ZIP_PATH = f"{FTP_ROOT}/search/txt.zip"
FASTA_NAMES = (
    "UP000005640_9606.fasta",
    "UP000005640_9606_additional.fasta",
    "Homo_sapiens.GRCh38.pep.all.fa",
    "Homo_sapiens.GRCh38.pep.abinitio.fa",
)
CACHE_VERSION = 1


@dataclass(frozen=True)
class ZipEntry:
    name: str
    compressed_size: int
    uncompressed_size: int
    compression: int
    local_header_offset: int


def _connect() -> FTP_TLS:
    ftp = FTP_TLS(FTP_HOST, timeout=120)
    ftp.login("anonymous", "codex-analysis@openai.com")
    ftp.prot_p()
    return ftp


def _remote_size(path: str) -> int:
    ftp = _connect()
    try:
        size = ftp.size(path)
        if size is None:
            raise OSError(f"FTP server did not report a size for {path}")
        return size
    finally:
        ftp.close()


def _fetch_range(path: str, offset: int, length: int) -> bytes:
    ftp = _connect()
    socket = None
    try:
        socket = ftp.transfercmd(f"RETR {path}", rest=offset)
        output = bytearray()
        while len(output) < length:
            chunk = socket.recv(min(1024 * 1024, length - len(output)))
            if not chunk:
                break
            output.extend(chunk)
        return bytes(output)
    finally:
        if socket is not None:
            socket.close()
        ftp.close()


def _zip_entries(path: str) -> dict[str, ZipEntry]:
    file_size = _remote_size(path)
    tail_length = min(file_size, 131_072)
    tail = _fetch_range(path, file_size - tail_length, tail_length)
    eocd_position = tail.rfind(b"PK\x05\x06")
    if eocd_position < 0:
        raise ValueError("ZIP end-of-central-directory record not found")
    eocd = struct.unpack("<4s4H2LH", tail[eocd_position : eocd_position + 22])
    _, _, _, _, entry_count, central_size, central_offset, _ = eocd
    if entry_count == 0xFFFF or central_size == 0xFFFFFFFF or central_offset == 0xFFFFFFFF:
        locator_position = tail.rfind(b"PK\x06\x07", 0, eocd_position)
        if locator_position < 0:
            raise ValueError("ZIP64 locator not found")
        _, _, zip64_offset, _ = struct.unpack(
            "<4sLQL", tail[locator_position : locator_position + 20]
        )
        zip64 = _fetch_range(path, zip64_offset, 56)
        (
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            entry_count,
            central_size,
            central_offset,
        ) = struct.unpack("<4sQ2H2L4Q", zip64)

    central = _fetch_range(path, central_offset, central_size)
    entries: dict[str, ZipEntry] = {}
    position = 0
    for _ in range(entry_count):
        if central[position : position + 4] != b"PK\x01\x02":
            raise ValueError("Invalid central directory entry")
        values = struct.unpack(
            "<4s6H3L5H2L", central[position : position + 46]
        )
        (
            _,
            _,
            _,
            _,
            compression,
            _,
            _,
            _,
            compressed_size,
            uncompressed_size,
            name_length,
            extra_length,
            comment_length,
            disk_number,
            _,
            _,
            local_offset,
        ) = values
        name_start = position + 46
        name = central[name_start : name_start + name_length].decode("utf-8")
        extra = central[
            name_start + name_length : name_start + name_length + extra_length
        ]
        if (
            uncompressed_size == 0xFFFFFFFF
            or compressed_size == 0xFFFFFFFF
            or local_offset == 0xFFFFFFFF
            or disk_number == 0xFFFF
        ):
            extra_position = 0
            while extra_position + 4 <= len(extra):
                extra_id, extra_size = struct.unpack(
                    "<HH", extra[extra_position : extra_position + 4]
                )
                data = extra[
                    extra_position + 4 : extra_position + 4 + extra_size
                ]
                extra_position += 4 + extra_size
                if extra_id != 1:
                    continue
                value_position = 0
                if uncompressed_size == 0xFFFFFFFF:
                    uncompressed_size = struct.unpack_from("<Q", data, value_position)[0]
                    value_position += 8
                if compressed_size == 0xFFFFFFFF:
                    compressed_size = struct.unpack_from("<Q", data, value_position)[0]
                    value_position += 8
                if local_offset == 0xFFFFFFFF:
                    local_offset = struct.unpack_from("<Q", data, value_position)[0]
                break
        entries[name] = ZipEntry(
            name,
            int(compressed_size),
            int(uncompressed_size),
            int(compression),
            int(local_offset),
        )
        position += 46 + name_length + extra_length + comment_length
    return entries


def _remote_zip_lines(path: str, entry: ZipEntry) -> Iterator[str]:
    local_header = _fetch_range(path, entry.local_header_offset, 30)
    values = struct.unpack("<4s5H3L2H", local_header)
    if values[0] != b"PK\x03\x04":
        raise ValueError(f"Invalid local ZIP header for {entry.name}")
    name_length, extra_length = values[-2:]
    data_offset = entry.local_header_offset + 30 + name_length + extra_length

    ftp = _connect()
    socket = ftp.transfercmd(f"RETR {path}", rest=data_offset)
    decompressor = zlib.decompressobj(-15) if entry.compression == 8 else None
    if entry.compression not in (0, 8):
        socket.close()
        ftp.close()
        raise ValueError(f"Unsupported ZIP compression method {entry.compression}")
    remaining = entry.compressed_size
    buffer = b""
    try:
        while remaining:
            compressed = socket.recv(min(1024 * 1024, remaining))
            if not compressed:
                raise OSError(f"Unexpected end of FTP stream for {entry.name}")
            remaining -= len(compressed)
            data = decompressor.decompress(compressed) if decompressor else compressed
            buffer += data
            while b"\n" in buffer:
                line, buffer = buffer.split(b"\n", 1)
                yield line.rstrip(b"\r").decode("utf-8", "replace")
        if decompressor:
            buffer += decompressor.flush()
        if buffer:
            yield buffer.rstrip(b"\r").decode("utf-8", "replace")
    finally:
        socket.close()
        ftp.close()


def _download_ftp_file(remote_path: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    ftp = _connect()
    try:
        with partial.open("wb") as handle:
            ftp.retrbinary(f"RETR {remote_path}", handle.write, blocksize=1024 * 1024)
        os.replace(partial, destination)
    finally:
        ftp.close()


def prepare_fastas(cache_dir: Path, progress: Callable[[str], None]) -> list[Path]:
    fasta_dir = cache_dir / "nature_fastas"
    paths: list[Path] = []
    for name in FASTA_NAMES:
        destination = fasta_dir / name
        if not destination.exists():
            progress(f"Downloading Nature reference FASTA: {name}")
            _download_ftp_file(f"{FTP_ROOT}/sequence/{name}", destination)
        paths.append(destination)
    return paths


def _enzyme_from_experiment(column: str) -> str | None:
    for enzyme in NATURE_ENZYMES:
        token = enzyme.replace("-", "")
        if f"_{token}_" in column:
            return enzyme
    return None


def prepare_peptide_cache(
    cache_dir: Path,
    entries: dict[str, ZipEntry],
    progress: Callable[[str], None],
) -> Path:
    destination = cache_dir / f"nature_peptides_v{CACHE_VERSION}.tsv.gz"
    if destination.exists():
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    lines = _remote_zip_lines(ZIP_PATH, entries["txt/peptides.txt"])
    header = next(lines).split("\t")
    index = {name: i for i, name in enumerate(header)}
    experiment_columns: list[tuple[int, int]] = []
    enzyme_bit = {enzyme: 1 << i for i, enzyme in enumerate(NATURE_ENZYMES)}
    for column_index, column in enumerate(header):
        if column.startswith("Experiment "):
            enzyme = _enzyme_from_experiment(column)
            if enzyme:
                experiment_columns.append((column_index, enzyme_bit[enzyme]))

    kept = 0
    with gzip.open(partial, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("id", "sequence", "enzyme_bits"))
        for row_number, line in enumerate(lines, start=1):
            fields = line.split("\t")
            if fields[index["Reverse"]] or fields[index["Potential contaminant"]]:
                continue
            bits = 0
            for column_index, bit in experiment_columns:
                if fields[column_index] not in ("", "0"):
                    bits |= bit
            if not bits:
                continue
            writer.writerow((fields[index["id"]], fields[index["Sequence"]], bits))
            kept += 1
            if row_number % 100_000 == 0:
                progress(f"Reduced {row_number:,} peptide rows ({kept:,} retained)")
    os.replace(partial, destination)
    progress(f"Retained {kept:,} target peptides from the Nature MaxQuant output")
    return destination


def _fasta_identifiers(header: str) -> list[str]:
    first = header.split()[0]
    identifiers = [first]
    parts = first.split("|")
    if len(parts) >= 3:
        identifiers.extend((parts[1], parts[2]))
    return identifiers


def load_reference_sequences(paths: list[Path]) -> dict[str, str]:
    sequences: dict[str, str] = {}
    for path in paths:
        header: str | None = None
        parts: list[str] = []

        def flush() -> None:
            if header is None:
                return
            sequence = "".join(parts).upper()
            for identifier in _fasta_identifiers(header):
                sequences.setdefault(identifier, sequence)

        with path.open("r", encoding="utf-8") as handle:
            for raw in handle:
                line = raw.strip()
                if line.startswith(">"):
                    flush()
                    header = line[1:]
                    parts = []
                elif line:
                    parts.append(line)
        flush()
    return sequences


def load_peptides(path: Path) -> tuple[list[str | None], bytearray]:
    sequences: list[str | None] = []
    bits = bytearray()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            peptide_id = int(row["id"])
            if peptide_id >= len(sequences):
                extension = peptide_id + 1 - len(sequences)
                sequences.extend([None] * extension)
                bits.extend(b"\x00" * extension)
            sequences[peptide_id] = row["sequence"]
            bits[peptide_id] = int(row["enzyme_bits"])
    return sequences, bits


def _slug(name: str) -> str:
    return (
        name.lower()
        .replace("trypsin/lys-c", "trypsin_lysc")
        .replace(" + ", "__")
        .replace("-", "")
        .replace(" ", "_")
    )


def prepare_per_protein_coverage(
    cache_dir: Path,
    entries: dict[str, ZipEntry],
    fasta_paths: list[Path],
    peptide_cache: Path,
    designs: list[tuple[str, tuple[str, ...]]],
    progress: Callable[[str], None],
) -> Path:
    destination = cache_dir / f"nature_per_protein_v{CACHE_VERSION}.tsv.gz"
    if destination.exists():
        return destination

    progress("Loading the Nature reference sequences and reduced peptide table")
    reference = load_reference_sequences(fasta_paths)
    peptide_sequences, peptide_bits = load_peptides(peptide_cache)
    enzyme_bit = {enzyme: 1 << i for i, enzyme in enumerate(NATURE_ENZYMES)}
    design_columns = [(_slug(name), name, enzymes) for name, enzymes in designs]

    lines = _remote_zip_lines(ZIP_PATH, entries["txt/proteinGroups.txt"])
    header = next(lines).split("\t")
    index = {name: i for i, name in enumerate(header)}
    partial = destination.with_suffix(destination.suffix + ".part")
    retained = missing_sequence = missing_peptide_match = 0
    with gzip.open(partial, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        columns = ["group_id", "major_protein", "length"]
        for slug, _, _ in design_columns:
            columns.extend((f"{slug}_covered", f"{slug}_pct"))
        writer.writerow(columns)

        for line in lines:
            fields = line.split("\t")
            if (
                fields[index["Only identified by site"]]
                or fields[index["Reverse"]]
                or fields[index["Potential contaminant"]]
            ):
                continue
            majority = fields[index["Majority protein IDs"]].split(";")
            major = majority[0]
            sequence = reference.get(major)
            if sequence is None:
                for candidate in majority[1:]:
                    sequence = reference.get(candidate)
                    if sequence is not None:
                        major = candidate
                        break
            if not sequence:
                missing_sequence += 1
                continue

            masks = {enzyme: 0 for enzyme in NATURE_ENZYMES}
            peptide_ids = fields[index["Peptide IDs"]]
            for token in peptide_ids.split(";") if peptide_ids else ():
                peptide_id = int(token)
                if peptide_id >= len(peptide_sequences):
                    continue
                peptide = peptide_sequences[peptide_id]
                bits = peptide_bits[peptide_id]
                if not peptide or not bits:
                    continue
                start = sequence.find(peptide)
                if start < 0:
                    missing_peptide_match += 1
                    continue
                peptide_mask = 0
                while start >= 0:
                    peptide_mask |= interval_mask(start, start + len(peptide))
                    start = sequence.find(peptide, start + 1)
                for enzyme, bit in enzyme_bit.items():
                    if bits & bit:
                        masks[enzyme] |= peptide_mask

            output: list[object] = [fields[index["id"]], major, len(sequence)]
            for _, _, enzymes in design_columns:
                combined = 0
                for enzyme in enzymes:
                    combined |= masks[enzyme]
                covered = combined.bit_count()
                output.extend((covered, 100.0 * covered / len(sequence)))
            writer.writerow(output)
            retained += 1
            if retained % 2_000 == 0:
                progress(f"Calculated observed coverage for {retained:,} protein groups")
    os.replace(partial, destination)
    progress(
        f"Nature coverage cache: {retained:,} protein groups; "
        f"{missing_sequence:,} groups lacked a source sequence; "
        f"{missing_peptide_match:,} peptide/group mappings did not localize"
    )
    return destination


def aggregate_per_protein(
    path: Path, designs: list[tuple[str, tuple[str, ...]]]
) -> list[dict[str, object]]:
    accumulators = {
        name: {"covered": 0, "total": 0, "percentages": []}
        for name, _ in designs
    }
    n_proteins = 0
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            n_proteins += 1
            length = int(row["length"])
            for name, _ in designs:
                slug = _slug(name)
                accumulator = accumulators[name]
                accumulator["covered"] += int(row[f"{slug}_covered"])
                accumulator["total"] += length
                accumulator["percentages"].append(float(row[f"{slug}_pct"]))

    results: list[dict[str, object]] = []
    for name, enzymes in designs:
        accumulator = accumulators[name]
        covered = int(accumulator["covered"])
        total = int(accumulator["total"])
        percentages = accumulator["percentages"]
        results.append(
            {
                "combination": name,
                "enzymes": ";".join(enzymes),
                "n_proteins": n_proteins,
                "covered_residues": covered,
                "total_residues": total,
                "residue_weighted_pct": 100.0 * covered / total,
                "median_pct": median(percentages),
            }
        )
    return results


def nature_coverage(
    cache_dir: Path,
    designs: list[tuple[str, tuple[str, ...]]],
    progress: Callable[[str], None] = print,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    fasta_paths = prepare_fastas(cache_dir, progress)
    per_protein_path = cache_dir / f"nature_per_protein_v{CACHE_VERSION}.tsv.gz"
    peptide_path = cache_dir / f"nature_peptides_v{CACHE_VERSION}.tsv.gz"
    if not per_protein_path.exists():
        progress("Reading ZIP directory for PXD024364/MSV000086944")
        entries = _zip_entries(ZIP_PATH)
        peptide_path = prepare_peptide_cache(cache_dir, entries, progress)
        per_protein_path = prepare_per_protein_coverage(
            cache_dir, entries, fasta_paths, peptide_path, designs, progress
        )
    rows = aggregate_per_protein(per_protein_path, designs)
    provenance = {
        "dataset": "PXD024364 / MSV000086944",
        "source_zip": f"ftp://{FTP_HOST}{ZIP_PATH}",
        "maxquant_members": ["txt/peptides.txt", "txt/proteinGroups.txt"],
        "fasta_files": [path.name for path in fasta_paths],
        "fasta_sha256": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in fasta_paths
        },
        "peptide_cache": str(peptide_path),
        "per_protein_cache": str(per_protein_path),
    }
    return rows, provenance


def _label_cache_path(
    cache_dir: Path,
    designs: list[tuple[str, tuple[str, ...]]],
    labels: dict[str, frozenset[str]],
) -> Path:
    identity = repr(
        (
            CACHE_VERSION,
            designs,
            [(name, "".join(sorted(residues))) for name, residues in labels.items()],
        )
    )
    digest = hashlib.sha256(identity.encode("utf-8")).hexdigest()[:12]
    return cache_dir / f"nature_label_coverage_v{CACHE_VERSION}_{digest}.tsv.gz"


def _label_slug(label: str) -> str:
    return label.lower().replace("+", "_")


def prepare_label_coverage_cache(
    cache_dir: Path,
    entries: dict[str, ZipEntry],
    fasta_paths: list[Path],
    peptide_cache: Path,
    designs: list[tuple[str, tuple[str, ...]]],
    labels: dict[str, frozenset[str]],
    progress: Callable[[str], None],
) -> Path:
    """Map observed label-bearing peptides to each Nature protein group."""

    destination = _label_cache_path(cache_dir, designs, labels)
    if destination.exists():
        return destination

    progress("Loading Nature sequences and observed peptides for SILAC masks")
    reference = load_reference_sequences(fasta_paths)
    peptide_sequences, peptide_bits = load_peptides(peptide_cache)
    peptide_labels: list[tuple[str, ...] | None] = [None] * len(peptide_sequences)
    for peptide_id, peptide in enumerate(peptide_sequences):
        if peptide:
            residues = frozenset(peptide)
            peptide_labels[peptide_id] = tuple(
                label
                for label, label_residues in labels.items()
                if not residues.isdisjoint(label_residues)
            )

    enzyme_bit = {enzyme: 1 << i for i, enzyme in enumerate(NATURE_ENZYMES)}
    design_columns = [(_slug(name), name, enzymes) for name, enzymes in designs]
    lines = _remote_zip_lines(ZIP_PATH, entries["txt/proteinGroups.txt"])
    header = next(lines).split("\t")
    index = {name: i for i, name in enumerate(header)}
    partial = destination.with_suffix(destination.suffix + ".part")
    retained = missing_sequence = missing_peptide_match = 0

    with gzip.open(partial, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        columns = ["group_id", "major_protein", "length"]
        for slug, _, _ in design_columns:
            for label in labels:
                token = f"{slug}_{_label_slug(label)}"
                columns.extend((f"{token}_covered", f"{token}_pct"))
        writer.writerow(columns)

        for line in lines:
            fields = line.split("\t")
            if (
                fields[index["Only identified by site"]]
                or fields[index["Reverse"]]
                or fields[index["Potential contaminant"]]
            ):
                continue
            majority = fields[index["Majority protein IDs"]].split(";")
            major = majority[0]
            sequence = reference.get(major)
            if sequence is None:
                for candidate in majority[1:]:
                    sequence = reference.get(candidate)
                    if sequence is not None:
                        major = candidate
                        break
            if not sequence:
                missing_sequence += 1
                continue

            masks = {
                enzyme: {label: 0 for label in labels}
                for enzyme in NATURE_ENZYMES
            }
            peptide_ids = fields[index["Peptide IDs"]]
            for token in peptide_ids.split(";") if peptide_ids else ():
                peptide_id = int(token)
                if peptide_id >= len(peptide_sequences):
                    continue
                peptide = peptide_sequences[peptide_id]
                bits = peptide_bits[peptide_id]
                carried_labels = peptide_labels[peptide_id]
                if not peptide or not bits or not carried_labels:
                    continue
                start = sequence.find(peptide)
                if start < 0:
                    missing_peptide_match += 1
                    continue
                peptide_mask = 0
                while start >= 0:
                    peptide_mask |= interval_mask(start, start + len(peptide))
                    start = sequence.find(peptide, start + 1)
                for enzyme, bit in enzyme_bit.items():
                    if bits & bit:
                        for label in carried_labels:
                            masks[enzyme][label] |= peptide_mask

            output: list[object] = [fields[index["id"]], major, len(sequence)]
            for _, _, enzymes in design_columns:
                for label in labels:
                    combined = 0
                    for enzyme in enzymes:
                        combined |= masks[enzyme][label]
                    covered = combined.bit_count()
                    output.extend((covered, 100.0 * covered / len(sequence)))
            writer.writerow(output)
            retained += 1
            if retained % 2_000 == 0:
                progress(
                    f"Calculated observed SILAC coverage for {retained:,} protein groups"
                )

    os.replace(partial, destination)
    progress(
        f"Nature SILAC cache: {retained:,} protein groups; "
        f"{missing_sequence:,} lacked a source sequence; "
        f"{missing_peptide_match:,} peptide/group mappings did not localize"
    )
    return destination


def aggregate_label_coverage(
    path: Path,
    designs: list[tuple[str, tuple[str, ...]]],
    labels: dict[str, frozenset[str]],
) -> list[dict[str, object]]:
    accumulators = {
        (name, label): {"covered": 0, "total": 0, "percentages": []}
        for name, _ in designs
        for label in labels
    }
    n_proteins = 0
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            n_proteins += 1
            length = int(row["length"])
            for name, _ in designs:
                slug = _slug(name)
                for label in labels:
                    token = f"{slug}_{_label_slug(label)}"
                    accumulator = accumulators[(name, label)]
                    accumulator["covered"] += int(row[f"{token}_covered"])
                    accumulator["total"] += length
                    accumulator["percentages"].append(float(row[f"{token}_pct"]))

    rows: list[dict[str, object]] = []
    for name, enzymes in designs:
        for label in labels:
            accumulator = accumulators[(name, label)]
            covered = int(accumulator["covered"])
            total = int(accumulator["total"])
            rows.append(
                {
                    "combination": name,
                    "enzymes": ";".join(enzymes),
                    "label": label,
                    "n_proteins": n_proteins,
                    "covered_residues": covered,
                    "total_residues": total,
                    "residue_weighted_pct": 100.0 * covered / total,
                    "median_pct": median(accumulator["percentages"]),
                }
            )
    return rows


def nature_label_coverage(
    cache_dir: Path,
    designs: list[tuple[str, tuple[str, ...]]],
    labels: dict[str, frozenset[str]],
    progress: Callable[[str], None] = print,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Observed residue coverage from peptides carrying each SILAC label."""

    cache_dir.mkdir(parents=True, exist_ok=True)
    fasta_paths = prepare_fastas(cache_dir, progress)
    peptide_path = cache_dir / f"nature_peptides_v{CACHE_VERSION}.tsv.gz"
    destination = _label_cache_path(cache_dir, designs, labels)
    if not peptide_path.exists() or not destination.exists():
        progress("Reading ZIP directory for Nature SILAC reference")
        entries = _zip_entries(ZIP_PATH)
        peptide_path = prepare_peptide_cache(cache_dir, entries, progress)
        destination = prepare_label_coverage_cache(
            cache_dir,
            entries,
            fasta_paths,
            peptide_path,
            designs,
            labels,
            progress,
        )
    rows = aggregate_label_coverage(destination, designs, labels)
    return rows, {
        "dataset": "PXD024364 / MSV000086944",
        "label_cache": str(destination),
        "peptide_cache": str(peptide_path),
        "search_specificity": "specific cleavage",
    }
