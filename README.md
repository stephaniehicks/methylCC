# methylCC

### Citation 

Hicks SC, Irizarry RA. (2019). _Genome Biology_ (accepted and in press)

### Why use methylCC? 

This is a package to estimate the cell composition 
    of whole blood in DNA methylation measured on any 
    platform technology (e.g. Illumina 450K microarray, 
    whole genome bisulfite sequencing (WGBS) and 
    reduced representation bisulfite sequencing (RRBS)). 

For help with the **methylCC** R-package, there is a vignette available in the /vignettes folder.
  
### Installing methylCC

The R-package **methylCC** can be installed from Github using the R 
package [devtools](https://github.com/hadley/devtools): 


Use  to install the latest version of **methylCC** from Github:
```s
library(devtools)
install_github("stephaniehicks/methylCC")
```

It can also be installed using Bioconductor: 

```s
# install BiocManager from CRAN (if not already installed)
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

# install methylCC package
BiocManager::install("methylCC")
```

After installation, the package can be loaded into R.
```s
library(methylCC)
```

# Bug reports
Report bugs as issues on the [GitHub repository](https://github.com/stephaniehicks/methylCC)


# Authors

* [Stephanie C. Hicks](https://github.com/stephaniehicks)
* [Rafael A. Irizarry](https://github.com/ririzarr)
