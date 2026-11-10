"""Preserve generated products before replacement, including failed attempts.

These are independent copies: hard links would allow a tool truncating a current
file to change its history. Failure to preserve a prior product aborts publication.
"""
from contextlib import contextmanager
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import uuid


def new_attempt(root, stage):
    if not stage or any(c not in "abcdefghijklmnopqrstuvwxyz0123456789-" for c in stage):
        raise ValueError("Invalid retained stage name")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    path = Path(root).resolve() / "retained" / stage / (stamp + "-" + uuid.uuid4().hex)
    path.mkdir(parents=True, exist_ok=False)
    return path


def write_receipt(path, value):
    with Path(path).open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


@contextmanager
def retained_attempt(root, stage):
    path = new_attempt(root, stage)
    write_receipt(path / "started.json", {"stage": stage, "started_utc": datetime.now(timezone.utc).isoformat()})
    try:
        yield path
    except BaseException as error:
        write_receipt(path / "failed.json", {"status": "failed", "error_type": type(error).__name__, "error": str(error)})
        raise
    else:
        write_receipt(path / "completed.json", {"status": "completed"})


def _copy(source, target):
    if source.is_symlink():
        target.symlink_to(os.readlink(source), target_is_directory=source.is_dir())
    elif source.is_dir():
        shutil.copytree(source, target, symlinks=True)
    else:
        shutil.copy2(source, target)


def publish(source, destination):
    """Publish a durable attempt's product, preserving any current destination."""
    source, destination = Path(source), Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with retained_attempt(destination.parent, "publication") as attempt:
        incoming = attempt / "incoming"
        _copy(source, incoming)
        if destination.exists() or destination.is_symlink():
            # Rename rather than delete directories; no historical bytes are lost.
            destination.rename(attempt / "previous")
        try:
            incoming.rename(destination)
        except BaseException:
            if (attempt / "previous").exists() or (attempt / "previous").is_symlink():
                (attempt / "previous").rename(destination)
            raise
        write_receipt(attempt / "publication.json", {"source": str(source.resolve()), "destination": str(destination.resolve())})
