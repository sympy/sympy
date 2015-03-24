"""
General SymPy exceptions and warnings.
"""

from __future__ import print_function, division

import warnings

from sympy.utilities.misc import filldedent


class SymPyDeprecationWarning(DeprecationWarning):
    r"""A warning for deprecated features of SymPy.

    This class is expected to be used with the warnings.warn function (note
    that one has to explicitly turn on deprecation warnings):

    >>> import warnings
    >>> from sympy.utilities.exceptions import SymPyDeprecationWarning
    >>> warnings.simplefilter(
    ...     "always", SymPyDeprecationWarning)
    >>> warnings.warn(
    ...     "Don't do this, it's deprecated",
    ...     SymPyDeprecationWarning) #doctest:+SKIP
    __main__:1: SymPyDeprecationWarning: "Don't do this, it's deprecated"

    The recommended way to use this class is, however, is by calling
    the warn method after constructing the message:

        >>> SymPyDeprecationWarning("Don't do this, it's deprecated.").warn() #doctest:+SKIP
        __main__:1: SymPyDeprecationWarning:

        Don't do this, it's deprecated.

          warning (see_above, SymPyDeprecationWarning)

    To provide additional information, create an instance of this
    class in this way:

    >>> SymPyDeprecationWarning(
    ...     feature="Such and such",
    ...     last_supported_version="1.2.3",
    ...     useinstead="this other feature")
    Such and such has been deprecated. It will be last supported in SymPy
    version 1.2.3. Use this other feature instead.

    Note that the text in ``feature`` begins a sentence, so if it begins with
    a plain English word, the first letter of that word should be capitalized.

    Either (or both) of the arguments ``last_supported_version`` and
    ``useinstead`` can be omitted. In this case the corresponding sentence
    will not be shown:

    >>> SymPyDeprecationWarning(feature="Such and such",
    ...     useinstead="this other feature")
    Such and such has been deprecated. Use this other feature instead.

    You can still provide the argument value.  If it is a string, it
    will be appended to the end of the message:

    >>> SymPyDeprecationWarning(
    ...     feature="Such and such",
    ...     useinstead="this other feature",
    ...     value="Contact the developers for further information.")
    Such and such has been deprecated. Use this other feature instead.
    Contact the developers for further information.

    If, however, the argument value does not hold a string, a string
    representation of the object will be appended to the message:

    >>> SymPyDeprecationWarning(
    ...     feature="Such and such",
    ...     useinstead="this other feature",
    ...     value=[1,2,3])
    Such and such has been deprecated. Use this other feature instead.
    ([1, 2, 3])

    To associate an issue with a deprecation, use the ``issue`` flag.

    >>> SymPyDeprecationWarning(
    ...    feature="Old feature",
    ...    useinstead="new feature",
    ...    issue=5241)
    Old feature has been deprecated. Use new feature instead. See
    https://github.com/sympy/sympy/issues/5241 for more info.

    Every formal deprecation should have an associated issue in the GitHub
    issue tracker.  All such issues should have the DeprecationRemoval
    tag.

    Additionally, each formal deprecation should mark the first release for
    which it was deprecated.  Use the ``deprecated_since_version`` flag for
    this.

    >>> SymPyDeprecationWarning(
    ...    feature="Old feature",
    ...    useinstead="new feature",
    ...    deprecated_since_version="0.7.2")
    Old feature has been deprecated since SymPy 0.7.2. Use new feature
    instead.

    Note that it may be necessary to go back through all the deprecations
    before a release to make sure that the version number is correct.  So just
    use what you believe will be the next release number (this usually means
    bumping the minor number by one).

    To mark a function as deprecated, you can use the decorator
    @deprecated.

    See Also
    ========
    sympy.core.decorators.deprecated

    """

    def __init__(self, value=None, feature=None, last_supported_version=None,
                 useinstead=None, issue=None, deprecated_since_version=None):
        self.fullMessage = ""

        if feature:
            if deprecated_since_version:
                self.fullMessage = "%s has been deprecated since SymPy %s. " % \
                                   (feature, deprecated_since_version)
            else:
                self.fullMessage = "%s has been deprecated. " % feature

        if last_supported_version:
            self.fullMessage += ("It will be last supported in SymPy "
                "version %s. ") % last_supported_version
        if useinstead:
            self.fullMessage += "Use %s instead. " % useinstead
        if issue:
            self.fullMessage += ("See "
                "https://github.com/sympy/sympy/issues/%d for more "
                "info. ") % issue

        if value:
            if not isinstance(value, str):
                value = "(%s)" % repr(value)
            value = " " + value
        else:
            value = ""

        self.fullMessage += value

    def __str__(self):
        return '\n%s\n' % filldedent(self.fullMessage)

    def warn(self, stacklevel=2):
        see_above = self.fullMessage
        # the next line is what the user would see after the error is printed
        # if stacklevel was set to 1. If you are writting a wrapper around this,
        # increase the stacklevel accordingly.
        warnings.warn(see_above, SymPyDeprecationWarning, stacklevel=stacklevel)

# Python by default hides DeprecationWarnings, which we do not want.
warnings.simplefilter("once", SymPyDeprecationWarning)
