"""
General SymPy exceptions and warnings.
"""

class SymPyDeprecationWarning(DeprecationWarning):
    r"""A warning for deprecated features of SymPy.

    This class is expected to be used with the warnings.warn function
    (note that one has to explicitly turn on deprecation warnings):

    >>> import warnings
    >>> from sympy.utilities.exceptions import SymPyDeprecationWarning
    >>> warnings.simplefilter("always", SymPyDeprecationWarning)
    >>> warnings.warn("Don't do this, it's deprecated", SymPyDeprecationWarning) #doctest:+SKIP
    __main__:1: SymPyDeprecationWarning: "Don't do this, it's deprecated"

    To provide additional information, do

    >>> warnings.warn(SymPyDeprecationWarning(feature="such and such",
    ...     last_supported_version="1.2.3",
    ...     useinstead="this other feature")) #doctest:+SKIP
    __main__:3: SymPyDeprecationWarning: The feature such and such is
    deprecated.  It will be last supported in SymPy version 1.2.3.
    Use this other feature instead.

    Either (or both) of the arguments last_supported_version and
    useinstead can be omitted.  In this case the corresponding
    sentence will not be shown:

    >>> warnings.warn(SymPyDeprecationWarning(feature="such and such",
    ...     useinstead="this other feature")) #doctest:+SKIP
    __main__:2: SymPyDeprecationWarning: The feature such and such is
    deprecated.  Use this other feature instead.

    You can still provide the argument value.  If it is a string, it
    will be appended to the end of the message:

    >>> warnings.warn(SymPyDeprecationWarning(feature="such and such",
    ...     useinstead="this other feature",
    ...     value="Contact the developers for further information.")) #doctest:+SKIP
    __main__:3: SymPyDeprecationWarning: The feature such and such is
    deprecated.  Use this other feature instead.  Contact the developers
    for further information.

    If, however, the argument value does not hold a string, a string
    representation of the object will be appended to the message:

    >>> warnings.warn(SymPyDeprecationWarning(feature="such and such",
    ...     useinstead="this other feature",
    ...     value=[1,2,3])) #doctest:+SKIP
    __main__:3: SymPyDeprecationWarning: The feature such and such is
    deprecated.  Use this other feature instead.  ([1, 2, 3])

    To mark a function as deprecated, you can use the decorator
    @deprecated.

    See Also
    ========
    sympy.core.decorators.deprecated
    """

    def __init__(self, value=None, feature=None, last_supported_version=None,
                 useinstead=None):
        self.fullMessage=""

        if feature:
            self.fullMessage = "The feature %s is deprecated." % feature

            if last_supported_version:
                self.fullMessage += "  It will be last supported in SymPy version %s." % last_supported_version
            if useinstead:
                self.fullMessage += "  Use %s instead." % useinstead

        if self.fullMessage:
            # We should also handle a non-string "value".
            if isinstance(value, str):
                if value:
                    self.fullMessage += "  " + value
            elif value:
                self.fullMessage += "  (%s)" % repr(value)
        else:
            # No extended arguments; replicate the original behaviour.
            self.fullMessage = repr(value)

    def __str__(self):
        return self.fullMessage
