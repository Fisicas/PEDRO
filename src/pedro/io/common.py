from __future__ import annotations

import hashlib


class ReadError(Exception):
    """The input file cannot be used; the message says why and what to do instead."""


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()
