BayesCTDesign v0.6.0 (Release date: 2019-06-17)
==============

Changes:

* Added print() and plot() methods which have the same functionality and interface as print_table() and plot_table().  Methods print_table() and plot_table() will be removed prior to v1.0.0.
* Corrected the calculation of variance.  Except for Gaussian outcomes, variance is currently calculated using delta method and the variance of log effect.  In the previous submission, the reported variance was exactly exp(log effect variance).
