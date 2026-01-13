"""Functions that involve magic. """

def pollute(names: list[str], objects: list[object]) -> None:
    """Pollute the global namespace with symbols -> objects mapping. """
    from inspect import currentframe
    frame = currentframe().f_back.f_back

    try:
        for name, obj in zip(names, objects):
            frame.f_globals[name] = obj
    finally:
        del frame  # break cyclic dependencies as stated in inspect docs
      