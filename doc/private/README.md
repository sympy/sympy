# The Internal API reference

This is a subproject that should be build independently. It contains ``autosummary`` directives to 
generate documentation for all modules and classes.

## Genarating the Documentation

Generate the documents for each module and class in your local copy of SymPy by:

    $ make html

## Status

The subproject is still work-in-progress. When the definition of what constitutes the ``private`` and ``public`` API is settled, it will be modified accordingly.

## Future Work

### Linking

- The subproject needs a way to link to the main docs project.
- It is better to use ``intersphinx`` though using an index page is easier and fine too. The later may require moving the subproject to ``/src/`` instead from experience.
