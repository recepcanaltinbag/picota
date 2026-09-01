"""
Deterministic synthetic assembly graphs for PICOTA regression tests.

The graphs built here reproduce the topology that composite transposons (CTs)
create in a short-read assembly graph: a repeated IS element collapses into a
single node, and every distinct cargo it flanks becomes a separate branch that
returns to that node. A genome carrying N CTs that share one IS therefore yields
N bubbles through a single shared node -- the scenario PICOTA must resolve.

All sequences are generated from a fixed seed so tests are reproducible.
Segments carry LN/KC/dp tags so the same fixtures can be reused once depth
parsing lands (roadmap phase 1).
"""

import random

DNA = "ACGT"


def make_dna(length, seed):
    """Generate a deterministic pseudo-random DNA sequence."""
    rng = random.Random(seed)
    return "".join(rng.choice(DNA) for _ in range(length))


def mutate(seq, divergence, seed):
    """Return a copy of `seq` with each base substituted at probability `divergence`."""
    rng = random.Random(seed)
    out = []
    for base in seq:
        if rng.random() < divergence:
            out.append(rng.choice([b for b in DNA if b != base]))
        else:
            out.append(base)
    return "".join(out)


def write_gfa(path, segments, links, depths=None, kmer_size=95):
    """
    Write a GFA1 file.

    segments : dict of {segment_name: sequence}
    links    : list of (from_name, from_orient, to_name, to_orient) tuples
    depths   : optional dict of {segment_name: float coverage}
    """
    depths = depths or {}
    with open(path, "w") as handle:
        for name, seq in segments.items():
            length = len(seq)
            fields = [f"S\t{name}\t{seq}", f"LN:i:{length}"]
            if name in depths:
                depth = depths[name]
                # KC is the k-mer count; assemblers report it instead of depth.
                kmer_count = int(round(depth * max(length - kmer_size + 1, 1)))
                fields.append(f"KC:i:{kmer_count}")
                fields.append(f"dp:f:{depth}")
            handle.write("\t".join(fields) + "\n")
        for src, src_orient, dst, dst_orient in links:
            handle.write(f"L\t{src}\t{src_orient}\t{dst}\t{dst_orient}\t0M\n")
    return path


def shared_repeat_two_cargos(path, repeat_len=2500, cargo_a_len=900,
                             cargo_b_len=1000, repeat_depth=60.0,
                             cargo_depth=30.0, seed=1):
    """
    Build a graph with one repeated element flanked by two distinct cargos.

    Topology (both cargos return to the same repeat node)::

            .--> cargoA --.
        repeat             repeat
            `--> cargoB --'

    Biologically this is one IS present in several genomic copies, two of which
    carry different cargo. Ground truth: two distinct composite transposons.
    `repeat_depth` is set above `cargo_depth` because a multi-copy IS collapses
    into a single node whose read depth is the sum of its copies.
    """
    segments = {
        "repeat": make_dna(repeat_len, seed),
        "cargoA": make_dna(cargo_a_len, seed + 1),
        "cargoB": make_dna(cargo_b_len, seed + 2),
    }
    links = [
        ("repeat", "+", "cargoA", "+"),
        ("cargoA", "+", "repeat", "+"),
        ("repeat", "+", "cargoB", "+"),
        ("cargoB", "+", "repeat", "+"),
    ]
    depths = {
        "repeat": repeat_depth,
        "cargoA": cargo_depth,
        "cargoB": cargo_depth,
    }
    return write_gfa(path, segments, links, depths)


def shared_repeat_n_cargos(path, n_cargos, repeat_len=2500, cargo_len=900,
                           cargo_len_step=0, repeat_depth=None,
                           cargo_depth=30.0, seed=1):
    """
    Build a graph where `n_cargos` distinct cargos all return to one repeat node.

    This is the reviewer scenario in miniature: one IS present in many genomic
    copies, several of which flank different cargo. Ground truth is `n_cargos`
    distinct composite transposons.

    `cargo_len_step` grows each successive cargo, which is what makes the
    repeat/cargo length ratio -- and therefore the deduplication outcome -- vary
    across candidates. Leave it at 0 for the hardest (equal-length) case.

    The repeat node's depth defaults to n_cargos x cargo_depth, mirroring an
    assembler collapsing every copy of the IS into one node.
    """
    if repeat_depth is None:
        repeat_depth = n_cargos * cargo_depth

    segments = {"repeat": make_dna(repeat_len, seed)}
    depths = {"repeat": repeat_depth}
    links = []
    for i in range(n_cargos):
        name = f"cargo{i}"
        segments[name] = make_dna(cargo_len + i * cargo_len_step, seed + 100 + i)
        depths[name] = cargo_depth
        links.append(("repeat", "+", name, "+"))
        links.append((name, "+", "repeat", "+"))
    return write_gfa(path, segments, links, depths)
