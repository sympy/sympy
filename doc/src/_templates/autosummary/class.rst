{{ fullname | escape | underline}}

{# Ideally we would be able to use the autosummary for class methods as well,
but unfortunately, autosummary does not allow just documenting :members:
using the same methodology as autodoc, so if we were to use autosummary, we
would get documentation for a bunch of methods that we don't want to include.
#}

.. autoclass:: {{ fullname }}
   :members:
