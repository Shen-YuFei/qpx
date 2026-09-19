"""Encode numeric native spectrum identifiers without discarding components."""

import re

_NATIVE_COMPONENT = re.compile(r"([A-Za-z][A-Za-z0-9_]*)=(\d+)(?=\s|$)")


def scan_from_native_id(native_id: str) -> list[int]:
    """Keep numeric components in native order, including repeated values.

    Default Thermo controller components are omitted according to the QPX scan
    convention. An identifier without numeric components returns an empty list.
    """
    components = [(key.lower(), int(value)) for key, value in _NATIVE_COMPONENT.findall(native_id)]
    fields = dict(components)
    if len(components) == 3 and fields.get("controllertype") == 0 and fields.get("controllernumber") == 1 and "scan" in fields:
        return [fields["scan"]]
    return [value for _, value in components]
