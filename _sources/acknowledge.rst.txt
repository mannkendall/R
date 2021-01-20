
.. include:: ./substitutions.rst

Acknowledging |name|
====================

1. Only use lower case letters when mentioning |name|. Always include the language and version
   number. Ideally, you should also include the Digital Object Identifier (DOI) associated to the
   specific release you have been using:

   |name| |github| |DOI|

2. If |name| was useful for your research, please cite the dedicated article:

   `Collaud Coen et al., 2020 <https://amt.copernicus.org/articles/13/6945/2020/>`__

3. |name| relies on external R libraries that require & deserve to be acknowledged in their own
   right. The following LaTeX blurb is one way to do so:

   .. code-block:: latex

        This research has made use of \textit{mannkendall v1.0.0} \citep[DOI:10.5194/amt-13-6945-2020][]{CollaudCoen2020}
        R code. \textit{mannkendall} relies on the following R packages: \textit{signal}, \textit{magrittr}
