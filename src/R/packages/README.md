# Building Drop-seq R libraries
For each library in src/R/packages, you can install using `devtools::install_github`, e.g.

```
devtools::install_github("broadinstitute/Drop-seq", subdir="/src/R/packages/DropSeq.utilities")  
```

By default this will use the latest code in the master branch.  You may prefer to specify a tagged version with `ref=<tag>`. 