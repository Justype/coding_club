# install.packages("pak")
# See this page for more details:
#   https://r-in-production.org/packages.html

options(repos = c(
  P3M  = "https://p3m.dev/cran/__linux__/focal/latest", # cat /etc/os-release VERSION_CODENAME
  CRAN = "https://cloud.r-project.org"
))

pak::repo_add(CRAN = "PPM@latest") # override CRAN
pak::repo_add(P3M = "PPM@latest") # supplement P3M (prebuilt)

pak::pak("languageserver")
pak::pak("tidyverse")
