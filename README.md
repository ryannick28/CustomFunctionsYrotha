# CustomFunctionsYrotha

This is a collection of my self-made R-functions. Mostly for creating nice univariate, pairs and 3D-plots.

Install with:
> devtools::install_github(repo = "ryannick28/CustomFunctionsYrotha", ref = 'main')

Because the used 'rgl' package can sometimes be difficult to install, I moved it in the Description from 'Imports' to 'Suggests'. It is only needed for creating 3D-plots. There is also a newer 3d-plot function which is based on the 'plotly' package and, therefore, does not rely on the 'rgl' package.
