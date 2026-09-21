"""Preliminary glycerol-stock storage suggestions."""

import re

from ..models import Box, GlycerolStock


BOX_NAME_PATTERN = re.compile(r"^L(?P<level>\d+)\s+B(?P<number>\d+)$", re.IGNORECASE)
BOX_ROWS = tuple("ABCDEFGHI")
BOX_COLUMNS = tuple(range(1, 10))


def storage_box_identity(box):
    """Return the assembly level and box number encoded in a box name."""
    match = BOX_NAME_PATTERN.match(box.name.strip())
    if not match:
        return None
    return int(match.group("level")), int(match.group("number"))


def suggest_storage_positions(plasmids, boxes, occupied_positions=()):
    """Suggest free positions for plasmids ordered by ID within each level.

    Box names follow the current convention, e.g. ``L2 B1``. Suggestions are
    deliberately not reservations: the caller must still save each GS and the
    existing position validation remains authoritative.
    """
    boxes_by_level = {}
    for box in boxes:
        identity = storage_box_identity(box)
        if identity is None:
            continue
        level, box_number = identity
        boxes_by_level.setdefault(level, []).append((box_number, box))

    for level_boxes in boxes_by_level.values():
        level_boxes.sort(key=lambda item: item[0])

    occupied = {
        (str(item[0]), str(item[1]), int(item[2]))
        for item in occupied_positions
    }
    sorted_plasmids = sorted(plasmids, key=lambda plasmid: (plasmid.idx is None, plasmid.idx or 0))
    suggestions = []
    positions_by_level = {}
    for plasmid in sorted_plasmids:
        level = plasmid.level
        level_positions = positions_by_level.setdefault(level, [])
        if not level_positions:
            for _, box in boxes_by_level.get(level, []):
                for row in BOX_ROWS:
                    for column in BOX_COLUMNS:
                        position = (str(box.id), row, column)
                        if position not in occupied:
                            level_positions.append((box, row, column, position))

        if level_positions:
            box, row, column, position = level_positions.pop(0)
            occupied.add(position)
            suggestions.append({
                "plasmid": plasmid,
                "level": level,
                "box": box,
                "box_row": row,
                "box_column": column,
                "position": f"{row}{column}",
                "available": True,
                "message": "Suggested position",
            })
        else:
            suggestions.append({
                "plasmid": plasmid,
                "level": level,
                "box": None,
                "box_row": None,
                "box_column": None,
                "position": None,
                "available": False,
                "message": (
                    f"No available box for L{level}" if level is not None
                    else "Assembly level is not defined"
                ),
            })
    return suggestions


def suggest_storage_positions_for_plasmids(plasmids):
    """Build suggestions using all existing glycerol-stock occupancy."""
    boxes = list(Box.objects.all())
    box_ids = [box.id for box in boxes]
    occupied_positions = GlycerolStock.objects.filter(box_id__in=box_ids).values_list(
        "box_id", "box_row", "box_column"
    )
    return suggest_storage_positions(plasmids, boxes, occupied_positions)
