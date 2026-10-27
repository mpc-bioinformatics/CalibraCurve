# CalibraCurve 


# 0.99.0
Initial submission to BioConductor.


# 1.1.1
Bugfixes

# 1.1.2
- Allow calculation of calibration curves with only one replicate per concentration level.
- Fix bug when input data is not ordered by concentration column
- Allow legend and add title to response factor plots
- Add new Logo to readme
- Fix bug when unequal numbers of replicates per concentration level are present

# 1.1.3
- small fix in documentation

# 1.1.4
- readDataTable() can now import a file containing data for multiple substances

# 1.1.5
- add warning if minReplicates is set to 1
- allows turning off log10 transformation of x-axis for plots
- allow usage of column names (as alternative to column numbers) when reading data
- option "multiplot" and "all_in_one" now also available for response factor plot
- add legend for response factor plot
- fix bug when multiple input files are used and at least one result cannot be computed
