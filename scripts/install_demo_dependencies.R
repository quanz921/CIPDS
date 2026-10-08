pkgs<-c("data.table","pROC","MASS","survey","survival","jsonlite")
missing<-pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))install.packages(missing,repos="https://cloud.r-project.org")
if(!requireNamespace("catboost",quietly=TRUE))stop("Install CatBoost R 1.2.5 using https://catboost.ai/en/docs/installation/r-installation-binary-installation ; then rerun.")
cat("Demo dependencies present.\n")

