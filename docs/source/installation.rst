.. include:: ./substitutions.rst

Installation
============

In a terminal, type:

.. code-block:: R

   devtools::install_github("mannkendall/R")

And that will take care of things. |name| uses `semantic versioning <https://semver.org/>`_.
The latest version is |version|.

The most recent release of |name| is also available for download from its
`Github repository <https://github.com/mannkendall/R/releases/latest/>`_.

Requirements
------------
|name| is compatible with the following R versions:

.. literalinclude:: ../../DESCRIPTION
    :language: R
    :lines: 7

Furthermore, |name| relies on a the following external modules, which will be automatically
installed:

.. literalinclude:: ../../DESCRIPTION
    :language: R
    :lines: 9

Testing the installation
------------------------

The most basic check to see if the installation was successful consists in checking the version of
the package:

.. code-block:: R

    library(mannkendall)
    packageVersion("mannkendall")
