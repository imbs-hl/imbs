
##' Migrates installed packages from old version
##' 
##' @title Migrate installed packages
##' @description Uses Bioconductor to install R packages, skips already installed packages
##' 
##' @param old_dir path to old package dir
##' @return Vector of packages which are not on CRAN or Bioconductor
##' @export
##' 
migrate_installed_packages <- function(old_dir){
  
  ## Get list of installed packages -----------------------------------
  to_migrate_pkgs <- dir(old_dir)
  
  ## Get info for CRAN and Bioc packages ------------------------------
  biocLite <- NULL
  biocVersion <- NULL
  message(date(), ' Querying CRAN')
  cran_pkgs <-   available.packages(
    contriburl = contrib.url('http://cran.r-project.org'))[,'Package']
      
  message(date(), ' Querying Bioconductor')  
  source("http://bioconductor.org/biocLite.R")
  biocLite()
  bioc_pkgs_all <-   available.packages(
    contriburl = contrib.url(paste0('http://bioconductor.org/packages/',
                                    biocVersion(), '/bioc')))[,'Package']

  bioc_pkgs <- setdiff(bioc_pkgs_all, 'BiocInstaller')
  
  ## Build list of packages to install from Bioconductor or CRAN
  to_install_pkgs <- intersect(to_migrate_pkgs, union(cran_pkgs, bioc_pkgs))
  
  ## Install packages
  lapply(to_install_pkgs, function(x){
    inst <- rownames(installed.packages())
    if (length(setdiff(x, inst) == 1)){
          biocLite(x, type = 'source')
    }
  })
  
  ## Not installed
  message('Not installed:')
  setdiff(to_migrate_pkgs,
          rownames(installed.packages())
          )
}
  
## old_pkg_dir <- '/usr/local/Cellar/r/3.1.3/R.framework/Resources/library'
## migrateInstalledPackages(old_pkg_dir) 
