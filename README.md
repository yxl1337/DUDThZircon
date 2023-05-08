# DUDThZircon

This repository was created to support the GCA article submission entitled "Temperature dependence of U and Th in igneous zircons" by Yuanyuan Liang et al.

---
### MATLAB code description

The data file is DvsTSiO2_source.xlsx.  This file contains our experimental data and a compilation of natural zircon DU/DTh estimates.

The main script is `linearRegressionCodeNM.m` which will call these other scripts:
- `importDataFromSpreadsheet.m` imports data from the spreadsheet
- `lnDvsTxal.m` calculates zircon saturation temperatures
- `InputMultipleRegression.m` adds estimated uncertainties for linear regression

The script `linearRegressionCodeNM.m` also calls a suite of MATLAB codes that perform linear regression through data with correlated uncertainties.  This list inclues
- `McLeanLinearRegression.m`
- `covmats4xl_v1.m`
- `myaa.m`
- `plotlinregNDxl_v3.m`
- `plotlinregxl_EnvelopeUnctMinMax.m`
- `plotlinregxl_XYLim.m`
- `plotlinregxl_tLim.m`
