# Drop-seq
Java and R tools for analyzing [Drop-seq](http://mccarrolllab.com/dropseq/) data

Drop-seq questions may be directed to [dropseq@gmail.com](mailto:dropseq@gmail.com).  
You may also use this address to be added to the Drop-seq Google group.

See [Releases](https://github.com/broadinstitute/Drop-seq/releases) to download binaries.

See [Drop-seq alignment cookbook](doc/Drop-seq_Alignment_Cookbook.pdf) for detailed usage of these tools.

See [Census-seq computational protocols](doc/Census-seq_Computational_Protcools.pdf) for detailed usage of Census-Seq tools.

# Building from source and installing

Download source:
```
git clone https://github.com/broadinstitute/Drop-seq.git 
cd Drop-seq
```

There are two options for building and installing Java code.  

- Build the executable jarfile, a single wrapper script from which all command-line programs can be invoked, installed in the cloned git sandbox.
- Build a zipfile containing the executable jarfile and wrapper scripts for all command-line programs.  The zipfile can be installed in an arbitrary location.

## Building the exectable jarfile and single wrapper script
`./gradlew installDist`

The wrapper script will be `dropseq/build/install/dropseq/bin/dropseq`.  `dropseq.bat` for Windows.
## Build zipfile containing executable jarfile and wrapper scripts
`./gradlew distZip`

### Installing
```
unzip -d <install-location> ./dropseq/build/distributions/dropseq-<version>.zip
```
Note that the name of the zipfile will be based on the tagged version and state of your git sandbox.

The files will be deployed into a subdirectory of `<install-location>` above, with name based on the name of the zipfile. 

## Building Drop-seq R libraries

You are encouraged to use the pre-built binaries in [Releases](https://github.com/broadinstitute/Drop-seq/releases), but 
if you want to build from (possibly unstable) sources, you can build from github source as follows:

For each library in src/R/packages, you can install using `devtools::install_github`, e.g.

```
devtools::install_github("broadinstitute/Drop-seq", subdir="/src/R/packages/DropSeq.utilities")  
```

Alternately, you can download sources and use `R CMD BUILD` and `R CMD INSTALL` as usual.

