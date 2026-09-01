"""
Deduplication of candidate cycles (roadmap phase 2).

The legacy deduplication deletes real composite transposons. Two stages are
responsible, and both encode the same wrong assumption -- that sharing a
repeated node makes two candidates duplicates -- when a shared repeat is exactly
what defines a composite transposon:

  cycle_match_based_on_contig_id()  strips strand for its exact-duplicate test
      but keeps strand for its similarity test, and calls >70% shared node
      length a duplicate. When N distinct CTs share one IS node the reported
      output saturates at two, whatever N is (defect D2).

  filter_cycles_with_kmer()  compares k-mer *sets*, so an IS-cargo-IS composite
      transposon looks identical to the IS-cargo fragment inside it (D3), and
      divides by the new candidate's k-mer count, an asymmetric containment
      measure that makes the outcome depend on DFS traversal order (D4).

This module replaces both with measures that cannot delete a distinct CT:

  * Paths are duplicates only when they are the *same cycle* -- the same nodes
    in the same cyclic order, up to rotation and reverse-complement traversal.
    Two bubbles through a shared repeat are different cycles, so both survive.

  * Sequences are compared by multiset Jaccard over canonical, circular k-mers.
    Counting k-mers instead of setting them keeps copy number visible; the
    symmetric denominator removes the order dependence.

See docs/ROADMAP.md and tests/test_known_defects.py.
"""

import math
from collections import Counter, defaultdict

from src.cycle_kmer_hash import print_progress_bar

_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G",
               "a": "t", "t": "a", "g": "c", "c": "g",
               "N": "N", "n": "n"}


def reverse_complement(seq):
    return "".join(_COMPLEMENT.get(base, "N") for base in reversed(seq))


# ─── Path identity ───────────────────────────────────────────────────────────

def flip_orientation(node):
    """'14349+' -> '14349-'. Nodes without a sign are returned unchanged."""
    if node.endswith("+"):
        return node[:-1] + "-"
    if node.endswith("-"):
        return node[:-1] + "+"
    return node


def reverse_traversal(path):
    """
    The same cycle walked the other way round.

    A cycle a+ -> b+ -> c+ traversed in reverse is c- -> b- -> a-: the order
    reverses and every node is read on the opposite strand.
    """
    return [flip_orientation(node) for node in reversed(path)]


def canonical_path_key(path):
    """
    A key that is identical for every representation of one cycle.

    Invariant under rotation (a cycle has no start) and under reverse-complement
    traversal (the same physical cycle read the other way). Two paths share a key
    if and only if they are the same cycle -- which is the only condition under
    which discarding one of them is safe.
    """
    if not path:
        return ()

    candidates = []
    for variant in (list(path), reverse_traversal(path)):
        for offset in range(len(variant)):
            candidates.append(tuple(variant[offset:] + variant[:offset]))
    return min(candidates)


def paths_are_duplicates(path_a, path_b):
    """True only when both paths describe the same cycle."""
    return canonical_path_key(path_a) == canonical_path_key(path_b)


def dedup_paths(paths, progress=False):
    """
    Keep one representative per distinct cycle, preserving input order.

    Unlike the legacy path filter this never drops a path for merely sharing
    nodes with an accepted one, so N composite transposons built on the same IS
    yield N candidates.
    """
    seen = set()
    kept = []
    total = len(paths)
    for index, path in enumerate(paths, start=1):
        key = canonical_path_key(path)
        if key not in seen:
            seen.add(key)
            kept.append(path)
        if progress and total:
            print_progress_bar(index, total, prefix="Processing:", suffix="Complete")
    return kept


# ─── Sequence similarity ─────────────────────────────────────────────────────

def canonical_kmers(seq, k, circular=True):
    """
    Count the k-mers of a sequence, each folded onto its reverse complement.

    Folding makes the count strand-agnostic without building two separate sets.
    `circular` wraps the end onto the start, because a cycle has no origin: two
    rotations of one cycle must produce the same counts, and without wrapping
    they differ by the k-1 k-mers that straddle the cut.
    """
    if k <= 0 or len(seq) < k:
        return Counter()

    window = seq + seq[:k - 1] if circular and len(seq) >= k else seq
    counts = Counter()
    for i in range(len(window) - k + 1):
        kmer = window[i:i + k]
        rc = reverse_complement(kmer)
        counts[kmer if kmer <= rc else rc] += 1
    return counts


def multiset_jaccard(counts_a, counts_b):
    """
    Multiset Jaccard: sum of per-k-mer minima over sum of per-k-mer maxima.

    Symmetric, so the outcome cannot depend on which candidate was seen first,
    and multiplicity-aware, so a sequence containing a repeat twice does not
    look identical to one containing it once.
    """
    if not counts_a or not counts_b:
        return 0.0

    intersection = 0
    union = 0
    for kmer in counts_a.keys() | counts_b.keys():
        a = counts_a.get(kmer, 0)
        b = counts_b.get(kmer, 0)
        intersection += min(a, b)
        union += max(a, b)
    return intersection / union if union else 0.0


def estimated_ani(jaccard, k):
    """
    Convert a k-mer Jaccard index into an estimated nucleotide identity.

    The Mash transform. Raw k-mer overlap is a cliff whose position depends
    entirely on k -- two sequences at 99.5% identity share 81% of their 21-mers
    but only 48% of their 80-mers -- so a threshold on it is an artefact of the
    parameter rather than a statement about the sequences. The transform undoes
    that: the same pair estimates 99.49% and 99.47% identity at k=21 and k=80,
    which is a quantity a biologist can actually choose a threshold for.

    Returns a fraction in [0, 1].
    """
    if jaccard <= 0:
        return 0.0
    if jaccard >= 1:
        return 1.0
    return max(0.0, 1.0 + math.log(2 * jaccard / (1 + jaccard)) / k)


def length_ratio(length_a, length_b):
    """Shorter length over longer, or 0 when either is empty."""
    longest = max(length_a, length_b)
    return min(length_a, length_b) / longest if longest else 0.0


def filter_cycles_multiset(cycle_info_list, k_mer_sim, threshold_sim,
                           name_prefix_cycle, min_ani=None,
                           min_length_ratio=0.95):
    """
    Drop cycles that are near-identical to one already accepted.

    Drop-in replacement for filter_cycles_with_kmer with the same signature and
    naming behaviour; `threshold_sim` is still a percentage. What changes is the
    measure: multiset Jaccard over canonical circular k-mers rather than
    directional containment over k-mer sets.

    When `min_ani` is given (as a fraction, e.g. 0.99) the decision switches
    from raw Jaccard to two criteria that must both hold: the estimated
    nucleotide identity is at least `min_ani`, AND the two lengths are within
    `min_length_ratio` of each other. Both are needed. Identity alone would
    merge a composite transposon with the fragment inside it -- IS-cargo-IS
    against IS-cargo estimates 99.65% identity at k=80, comfortably above any
    sensible threshold -- because a copy-number difference is not divergence,
    and the transform cannot tell them apart. Length separates them: those two
    differ by 38% in size, while a genuine duplicate does not differ at all.

    As in the legacy filter, the inverted index only surfaces candidates that
    share at least one k-mer, so cycles with nothing in common are never
    compared no matter how low the threshold is set.
    """
    if not cycle_info_list:
        print('Print empty Cycle List')
        return cycle_info_list

    accepted_counts = []
    kmer_index = defaultdict(list)
    kept = []
    name_it = 0
    threshold = threshold_sim / 100.0

    for index, cycle_el in enumerate(cycle_info_list, start=1):
        reverse_ori = 'reverseoriented_' if cycle_el.reverseOriented else ''
        counts = canonical_kmers(cycle_el.sequence, k_mer_sim)

        is_duplicate = False
        if counts:
            # Only cycles sharing at least one k-mer can clear the threshold.
            candidate_indices = set()
            for kmer in counts:
                candidate_indices.update(kmer_index[kmer])
            for candidate in candidate_indices:
                other_counts, other_length = accepted_counts[candidate]
                jaccard = multiset_jaccard(counts, other_counts)
                if min_ani is None:
                    if jaccard >= threshold:
                        is_duplicate = True
                        break
                elif (estimated_ani(jaccard, k_mer_sim) >= min_ani
                      and length_ratio(cycle_el.length, other_length) >= min_length_ratio):
                    is_duplicate = True
                    break

        if not is_duplicate:
            name_it += 1
            cycle_el.name = f"{name_prefix_cycle}_{reverse_ori}{name_it}"
            kept.append(cycle_el)

            position = len(accepted_counts)
            accepted_counts.append((counts, cycle_el.length))
            for kmer in counts:
                kmer_index[kmer].append(position)

        print_progress_bar(index, len(cycle_info_list),
                           prefix='Processing:', suffix='Complete')

    print('\n')
    return kept
