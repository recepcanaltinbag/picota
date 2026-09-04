"""
Deciding whether a previous run's output can be reused.

The pipeline skips work whose output file already exists. That is what makes a
68-sample run resumable, and it is also how a result outlives the thing that
produced it: `run_blast` trusted any file at its output path, so a search
interrupted halfway left a truncated table that every later run accepted, and
re-running after a database was rebuilt or the scoring code changed quietly
reported the old answer. Nothing failed, which is why it survived.

Existence is the wrong question. The right one is whether the inputs that
produced the file are still the inputs today, so each output carries a sidecar
naming them:

    picota_final_tab
    picota_final_tab.picota-sig.json

The sidecar is written only after the work succeeds, so an interrupted run
leaves an output with no sidecar and the next run redoes it. Contents are
hashed rather than stat'ed: a checkout, a copy or an rsync moves mtime without
changing a byte, and recomputing a whole phase because git touched a file is
its own kind of wrong answer.
"""

import hashlib
import json
import os

SIGNATURE_SUFFIX = '.picota-sig.json'

# Above this, hashing costs more than the mistake it prevents, and a file this
# large is a database rather than something edited by hand: size and mtime
# carry the change. The bundled references are all far below it.
MAX_HASH_BYTES = 256 * 1024 * 1024


def file_signature(path):
    """
    What identifies this file's content, or None when it is not there.

    A missing input is itself part of the signature -- a search run when a
    database was absent must not be reused once the database appears.
    """
    try:
        stat = os.stat(path)
    except OSError:
        return None

    if stat.st_size > MAX_HASH_BYTES:
        return {'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns}

    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return {'size': stat.st_size, 'sha256': digest.hexdigest()}


def signature(inputs=(), extra=None):
    """
    The signature for one unit of work.

    `inputs` are the files it reads -- queries, databases, and the source of the
    code that interprets them. `extra` carries the scalars that change a result
    without changing a file: thresholds, the program invoked, the output format.
    """
    return {
        'inputs': {path: file_signature(path) for path in inputs},
        'extra': dict(extra or {}),
    }


def signature_path(output_path):
    return output_path + SIGNATURE_SUFFIX


def is_current(output_path, sig):
    """
    True when `output_path` exists and was produced from exactly `sig`.

    An output with no sidecar is never current. That deliberately re-does the
    first run after this check is introduced: those files were produced by code
    that recorded nothing, so there is no evidence they match anything.
    """
    if not os.path.exists(output_path):
        return False
    try:
        with open(signature_path(output_path)) as handle:
            return json.load(handle) == sig
    except (OSError, ValueError):
        return False


def record(output_path, sig):
    """Write the sidecar. Call only once the work has actually succeeded."""
    with open(signature_path(output_path), 'w') as handle:
        json.dump(sig, handle, indent=2, sort_keys=True)


def discard(output_path):
    """Drop an output and its sidecar, so the next run cannot half-trust it."""
    for path in (output_path, signature_path(output_path)):
        try:
            os.remove(path)
        except OSError:
            pass
