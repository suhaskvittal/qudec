"""Stable sparse placements and compact upper-triangle features for affinity inputs.

Matrix indices are anonymous feature slots, not detector IDs. Physical detector
slots occupy ``0 .. size-3``; the two boundary anchors are always the final two
slots. The default 48-slot layout has eleven four-slot buckets and two shared
stash slots, giving every detector ten placement candidates.
"""

import hashlib
import struct

import numpy as np


class LayoutOverflow(RuntimeError):
    """Raised when active detectors cannot be assigned to non-boundary slots."""


def _layout_parameters(size):
    size = int(size)
    if size < 12 or (size - 4) % 4:
        raise ValueError("size must reserve two boundaries, two stash slots, and four-slot buckets")
    bucket_count = (size - 4) // 4
    return size, size - 2, bucket_count


def build_candidates(coords, relevant, size=48, hashseed=17):
    """Build an ``N x 10`` candidate-slot table from ``(x,y,t,relevant)`` keys.

    Each coordinate key hashes to two distinct four-slot buckets. The last two
    physical slots are shared stash candidates. Only the first ``size-2`` slots
    can hold detectors; the final two matrix indices are boundary anchors.
    """
    size, physical_count, bucket_count = _layout_parameters(size)
    coords = np.asarray(coords, dtype=np.float64)
    relevant = np.asarray(relevant, dtype=bool)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError("coords must have shape (N, 3) with x, y, t columns")
    if relevant.shape != (len(coords),):
        raise ValueError("relevant must have shape (N,)")
    if not np.isfinite(coords).all():
        raise ValueError("all detector coordinates must be finite")
    candidates = np.empty((len(coords), 10), dtype=np.int32)
    seed_bytes = str(int(hashseed)).encode("ascii")
    for detector, (coordinate, is_relevant) in enumerate(zip(coords, relevant)):
        # Explicit little-endian float packing makes placement independent of
        # host byte order and avoids Python's randomized hash implementation.
        key = (seed_bytes + b"\0" + struct.pack("<dddB", *coordinate, int(is_relevant)))
        buckets = []
        for counter in (0, 1):
            digest = hashlib.blake2b(key + bytes([counter]), digest_size=8,
                                     person=b"affmat-v1").digest()
            bucket = int.from_bytes(digest, "little") % bucket_count
            if buckets and bucket == buckets[0]:
                bucket = (bucket + 1) % bucket_count
            buckets.append(bucket)
        first, second = buckets
        candidates[detector, :4] = np.arange(first * 4, first * 4 + 4)
        candidates[detector, 4:8] = np.arange(second * 4, second * 4 + 4)
        candidates[detector, 8:] = (physical_count - 2, physical_count - 1)
    return candidates


def deterministic_placement(ids, candidates, size=48):
    """Assign active IDs to physical slots with deterministic augmenting paths.

    IDs are processed in ascending order. Existing occupants can move only to
    another candidate slot of their own detector. The returned slot array aligns
    with the caller's original ``ids`` order. No detector is overwritten or
    dropped; impossible placements raise :class:`LayoutOverflow`.
    """
    size, physical_count, _ = _layout_parameters(size)
    ids = np.asarray(ids, dtype=np.int64)
    candidates = np.asarray(candidates)
    if ids.ndim != 1 or len(np.unique(ids)) != len(ids):
        raise ValueError("ids must be a one-dimensional list of unique detector IDs")
    if candidates.ndim != 2 or candidates.shape[1] != 10:
        raise ValueError("candidates must be an N x 10 array from build_candidates")
    if len(ids) and (ids.min() < 0 or ids.max() >= len(candidates)):
        raise ValueError("detector ID is outside the candidate table")
    if len(ids) > physical_count:
        raise LayoutOverflow(f"{len(ids)} active detectors exceed {physical_count} physical slots")
    rows = candidates[ids]
    if np.any(rows < 0) or np.any(rows >= physical_count):
        raise ValueError("candidate slots must be physical slots, before the boundary anchors")

    owner = np.full(physical_count, -1, dtype=np.int64)
    placement = np.full(len(candidates), -1, dtype=np.int64)

    def augment(detector, visited_slots):
        for slot in candidates[detector]:
            slot = int(slot)
            if slot in visited_slots:
                continue
            visited_slots.add(slot)
            previous = int(owner[slot])
            if previous == -1 or augment(previous, visited_slots):
                owner[slot] = detector
                placement[detector] = slot
                return True
        return False

    for detector in np.sort(ids):
        if not augment(int(detector), set()):
            raise LayoutOverflow(f"no candidate-slot assignment exists for detector {detector}")
    return placement[ids].astype(np.int32, copy=False)


def packed_placement(ids, size=48):
    """Return sequential physical slots aligned with IDs; reserve two final anchors."""
    size, physical_count, _ = _layout_parameters(size)
    ids = np.asarray(ids)
    if ids.ndim != 1:
        raise ValueError("ids must be one-dimensional")
    if len(ids) > physical_count:
        raise LayoutOverflow(f"{len(ids)} active detectors exceed {physical_count} physical slots")
    return np.arange(len(ids), dtype=np.int32)


def triangle_features(ids, affinity, boundary, scale, size=48, slots=None):
    """Return path-strength then connected-mask upper triangles as float32.

    ``affinity`` is the fixed symmetric-ish raw affinity table with a shared
    boundary in its final row/column. Interior values are symmetrized by the
    arithmetic mean before conversion to the existing clipped log-affinity
    strength. ``boundary`` must contain two rows of detector-to-physical-boundary
    natural-log distances, shape ``(2, N)``. Distinct boundary anchors are placed
    at ``size-2`` and ``size-1`` and are present even for the empty syndrome.

    The result length is ``2 * size * (size + 1) // 2``. The first half is the
    strength upper triangle including the diagonal; the second half is its
    connected/present mask in the same order.
    """
    size, physical_count, _ = _layout_parameters(size)
    ids = np.asarray(ids, dtype=np.int64)
    affinity = np.asarray(affinity)
    boundary = np.asarray(boundary, dtype=np.float64)
    if ids.ndim != 1 or len(np.unique(ids)) != len(ids):
        raise ValueError("ids must be a one-dimensional list of unique detector IDs")
    if affinity.ndim != 2 or affinity.shape[0] != affinity.shape[1] or affinity.shape[0] < 2:
        raise ValueError("affinity must be a square table with a shared-boundary row/column")
    detector_count = affinity.shape[0] - 1
    if len(ids) and (ids.min() < 0 or ids.max() >= detector_count):
        raise ValueError("detector ID is outside the affinity table")
    if boundary.shape != (2, detector_count):
        raise ValueError(f"boundary must have shape (2, {detector_count})")
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("scale must be finite and positive")
    if slots is None:
        slots = packed_placement(ids, size)
    else:
        slots = np.asarray(slots, dtype=np.int64)
    if slots.shape != ids.shape or (len(slots) and
            (slots.min() < 0 or slots.max() >= physical_count or len(np.unique(slots)) != len(slots))):
        raise ValueError("slots must be unique physical slots aligned with ids")

    strength = np.zeros((size, size), dtype=np.float64)
    mask = np.zeros((size, size), dtype=np.float64)
    # Boundary anchors and active detector nodes remain explicit for empty and
    # singleton syndromes, independent of any numerical detector ID.
    b0, b1 = size - 2, size - 1
    strength[b0, b0] = strength[b1, b1] = 1.
    mask[b0, b0] = mask[b1, b1] = 1.
    if len(ids):
        strength[slots, slots] = 1.
        mask[slots, slots] = 1.

        raw = np.asarray(affinity[np.ix_(ids, ids)], dtype=np.float64)
        raw = (raw + raw.T) * .5
        connected = raw > 0
        pair_strength = np.maximum(0., 1. + np.log(np.maximum(raw, 1e-300)) / scale)
        pair_strength = np.minimum(pair_strength, 1.)
        for i, slot_i in enumerate(slots):
            for j, slot_j in enumerate(slots):
                strength[slot_i, slot_j] = pair_strength[i, j] if connected[i, j] else 0.
                mask[slot_i, slot_j] = float(connected[i, j])
        # Raw affinity tables use a zero diagonal; active-node presence is an
        # independent mask signal, as in the original dense feature layout.
        strength[slots, slots] = 1.
        mask[slots, slots] = 1.

        distances = boundary[:, ids]
        boundary_mask = np.isfinite(distances)
        boundary_strength = np.zeros_like(distances)
        boundary_strength[boundary_mask] = np.clip(
            1. - distances[boundary_mask] / scale, 0., 1.)
        for i, slot in enumerate(slots):
            for b, boundary_slot in enumerate((b0, b1)):
                value = boundary_strength[b, i]
                present = float(boundary_mask[b, i])
                strength[slot, boundary_slot] = strength[boundary_slot, slot] = value
                mask[slot, boundary_slot] = mask[boundary_slot, slot] = present

    upper = np.triu_indices(size)
    return np.concatenate((strength[upper], mask[upper])).astype(np.float32, copy=False)
