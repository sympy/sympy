import cProfile

def profiling(fname):
    """ Profile the functions.
    Display its statistics.

    Parameters
    ==========

    fname : string
            It is the name of the function to be profiled.
            Pass the parameter similar to how the function 
            is invoked.
            Example: str(<class_instance>.<function_name(<parameters>))

    Returns
    =======

    void function : none type

    """

    cProfile.run(fname)