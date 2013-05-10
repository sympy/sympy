How to Build Documentation
==========================

To make the html documentation, install the prerequisites, e.g. on
Debian/Ubuntu (similarly for other distributions)::

    apt-get install python-sphinx texlive-latex-recommended dvipng

and do::

    make html

and to view it, do::

    epiphany _build/html/index.html


About the Translated Tutorial
=============================

Translated versions of the tutorial are generated when building the html
documentation using ``make html``. If you quickly want to check only the
translations and not the whole documentation, just run::

    make htmli18n

This will create ``tutorial.cs.html``, ``tutorial.ru.html`` (and so on for all
languages) in the ``_build/html/`` directory. The input is the English tutorial
``tutorial.en.rst`` and the .po files ``tutorial.cs.po``, ``tutorial.ru.po``, etc.
You can freely change the English tutorial -- sentences that are not translated
will remain in English in the translated versions.


How to Update Translations
==========================

In order to update translations, you first need to make sure that the
``tutorial.pot`` template is up-to-date by running::

    make gettext

If you are creating a translation for a new language, copy the generated
``tutorial.pot`` to a new file ``tutorial.??.po`` where ``??`` is the
two-character language code for your language. Also add the language
code to the LANGUAGES macro in the Makefile. When the translation work
for a new language has reached 90% or more, a link to the new translation
should be added at the bottom of tutorial.en.rst.

If you are just updating a translation, for example the
``tutorial.cs.po``, just do::

    make update-po

This will create a new ``tutorial.cs.po`` by using the template
``tutorial.pot`` and reusing the old translations from old ``tutorial.cs.po``
(if they still work) and leaving the rest untranslated.

Update your ``tutorial.??.po`` file with your translations, then just build it
using ``make html`` (see the previous section). When you are done, use
``git add`` to add your changes to the repository and submit a pull request.
