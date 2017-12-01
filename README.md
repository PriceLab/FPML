# Installation Instructions

The easiest way to install the package and all of its dependencies is through `devtools::install_github`. Basically, do this:

1. Install the `devtools` package

```
library(BiocInstaller)
biocLite("devtools")
```

2. Load the `devtools` package

```
library(devtools)
```

3. Run the `install_github` command using this repository

```
install_github("PriceLab/FPML")
```

That should install the package and its dependencies