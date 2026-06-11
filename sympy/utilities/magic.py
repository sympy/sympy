"""Functions that involve magic. """
from __future__ import annotations

def pollute(names: list[str], objects: list[object]) -> None:
    """Pollute the global namespace with symbols -> objects mapping. """
    from inspect import currentframe
    frame = currentframe()

    # Go back two frames to find the caller
    frame = frame.f_back if frame is not None else None
    frame = frame.f_back if frame is not None else None

    if frame is None:
        raise RuntimeError("Unable to get stack frame.")

    try:
        for name, obj in zip(names, objects):
            frame.f_globals[name] = obj
    finally:
        del frame  # break cyclic dependencies as stated in inspect docs
