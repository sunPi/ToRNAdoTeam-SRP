# using package_info for the extra info delivered
### -> you might need to change the lib.loc according to system/version
allpkgs <- (devtools:::package_info(rownames(installed.packages(lib.loc = "/Library/Frameworks/R.framework/Versions/3.5/Resources/library")),include_base = T))
allpkgs[grepl("Bioconductor",allpkgs$source),]

# in a similar fashion, for cran/CRAN packages, and leading edge packages you might have got from Github repos
allpkgs[grepl("cran",allpkgs$source,ignore.case = TRUE),]
allpkgs[grepl("Github",allpkgs$source),]
