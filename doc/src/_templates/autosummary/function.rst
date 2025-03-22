{{ fullname | escape | underline }}

.. currentmodule:: {{ module }}

.. autofunction:: {{ objname }}

{% block examples %}
{% endblock %}

.. minigallery:: {{ fullname }}
    :add-heading: Examples using ``{{ objname }}``

{% block seealso %}
{% endblock %}

{%- if sourcename %}
.. raw:: html

   <div class="clearer"></div>
   <div class="viewcode-block" id="{{ fullname }}"><a class="viewcode-back" href="{{ backlink }}">{{ _('Return to documentation') }}</a></div>
{%- endif %}
