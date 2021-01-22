.. include:: ./substitutions.rst

Running |name|
==============

The main function to compute the Mann-Kendall test with the desired prewhitening method, temporal
segmentation and confidence limit is:

.. code-block:: R

   MK.tempAggr <- function(data, PW.method = "3PW", resolution, alpha.mk = 95, alpha.cl = 90, alpha.xhomo = 90, alpha.ak = 95, seasonal = FALSE, seasons = NULL)

For example:

.. code-block:: R

    result <- MK.tempAggr(data = input, PW.method = '3PW', resolution = 0.01)
    result <- MK.tempAggr(data = input, PW.method = '3PW', resolution = 0.01, alpha.mk = 90, alpha.ak = 90, alpha.cl = 95, alpha.xhomo = 95)


Function description:
---------------------

The MK test and the Sen slope are applied on the given time granularity, temporal aggregation and
prewhitening (PW) method. Five prewhitening methods can be chosen, 3PW being the default option:

  - ``3PW`` (Collaud Coen et al., 2020): 3 prewhitening methods are applied (PW and TFPW_Y to determine the statistic significance (ss) of the MK test and the VCTFPW method to compute the Sen's slope
  - ``PW`` (prewhitened, Kulkarni and von Storch, 1995)
  - ``TFPW.Y`` (trend free PW,Yue et al., 2001)
  - ``TFPW.WS`` (trend free PW, Wang and Swail, 2001)
  - ``VCTFPW`` (variance corrected trend free PW, Wang et al., 2015)

For the PW, only statistically significant (ss) autocorrelation are taken into account. The default ss for the MK test is taken
at 95% confidence limit. The default ss for upper and lower confidence limits is 90% of the all
intervals differences distribution. The default ss for the autocorrelation coefficient is 95%. The
default ss for the homogeneity test between temporal aggregation of the MK test is 90%.
If seasonal Mann-Kendall is applied, the yearly trend is assigned only if the results of the seasonal test are homogeneous. The default ss for the homogeneity test between temporal aggregation of the seasonal MK test is 90%.


INPUT:
******
  - ``data`` (data.frame)  a data.frame with the first column being a **POSIXct** object with tz = "UTC" and the second column the variable to be analysed. In case a seasonal MK test has to be applied, then the ``seasonal`` flag must be set to ``TRUE`` and an extra vector for the split into seasons should be submitted (see ``seasons``) 
  - ``PW_method`` (character) =  the PW method to be used among 3PW, PW, TFPW.Y, TFPW.WS, VCTFPW. The default is ``3PW``
  - ``resolution`` (numeric) = interval to determine the number of ties. It should be similar to the resolution of the instrument. There is no default
  - ``alpha.mk`` (numeric) = confidence limit for Mk test in %. The default value is 95
  - ``alpha.cl`` (numeric) = confidence limit for the confidence limits of the Sen's slope in %. The default value is 90
  - ``alpha.xhomo`` (numeric) = confidence limit for the homogeneity between seasons in %. The default value is 90
  - ``alpha.ak`` (numeric) = confidence limit for the first lag autocorrelation in %. The default value is 95
  - ``seasonal`` (boolean) = ``TRUE`` if the analysis needs to be performed over user-defined seasons (default is ``FALSE``)
  - ``seasons`` (vector) = a vector of numerics, or characters, or factors needed to ``split`` the ``data`` into seasons. It is used only if ``seasonal = TRUE``


OUTPUT:
*******

``result`` (data.frame) comprises the following 5 variables:

  - ``slope`` (numeric) = Sen's slope in units/y
  - ``UCL`` (numeric) = upper confidence level in units/y
  - ``LCL`` (numeric) = lower confidence level in units/y
  - ``P`` (numeric) =  probability for the statistical significance. If 3PW is applied, P = max(P_PW, P_TFPW_Y)
  - ``ss`` (numeric) = statistical significance
        * ``alpha.mk`` if the test is ss at the alpha confidence level. Default = 95%
        * ``0`` if the test is not ss at the ``alpha.mk`` confidence level
        * ``-1`` if the test is a ``TFPW.Y`` fals epositive at ``alpha.mk`` confidence level
        * ``-2`` if the test is a PW false positive at ``alpha.mk`` confidence level

**Sources:**

  - Collaud Coen et al., Effects of the prewhitening method, the time granularity and the time segmentation on the Mann-Kendall trend detection and the associated Sen's slope,  Atmos. Meas. Tech., https://doi.org/10.5194/amt-13-6945-2020, 2020.

  - Sirois, A.: A brief and biased overview of time-series analysis of how to find that evasive trend, WMO/EMEP Workshop on Advanced Statistical Methods and Their Application to Air Quality Data Sets, Annex E., Global Atmosphere Watch No. 133, TD- No. 956, World Meteorological Organization, Geneva, Switzerland, 1998. annexe E, p. 26

  - Gilbert, R.: Statistical Methods for Environmental Pollution Monitoring, Van Nostrand Reinhold Company, New York, 1987.
    and the explanations about MULTMK/PARTMK de C. Libiseller
