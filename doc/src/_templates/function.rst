{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autofunction:: {{ objname }}

Examples
--------

.. plot::
   :context: close-figs

   >>> from sympy import *
   >>> from sympy.abc import x, y, z
   >>> init_printing(use_unicode=True)

.. include:: /{{ fullname }}_examples.rst

See Also
--------

.. seealso::

   {% for obj in see_also %}
   - :func:`~{{ obj }}`
   {% endfor %}

References
----------

.. rubric:: References

{% if references %}
{% for ref in references %}
.. [{{ ref.key }}] {{ ref.text }}
{% endfor %}
{% endif %} 