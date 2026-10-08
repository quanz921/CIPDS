suppressPackageStartupMessages(library(jsonlite))
files<-list.files(".",pattern="\\.[Rr]$",recursive=TRUE,full.names=TRUE)
files<-files[!grepl("demo_work",files)]
errors<-list()
for(f in files)tryCatch(parse(f,encoding="UTF-8"),error=function(e){errors[[f]]<<-conditionMessage(e)})
report<-list(scope="R syntax only; no research pipelines executed",files=length(files),errors=errors,passed=length(errors)==0)
write_json(report,"provenance/r_syntax_validation.json",pretty=TRUE,auto_unbox=TRUE)
cat("R syntax:",length(files),"files;",length(errors),"errors\n")
if(length(errors))quit(status=1)
