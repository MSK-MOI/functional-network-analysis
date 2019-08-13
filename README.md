<http://github.com/jimmymathews/grape-plots>

#Gene Reduction/Annotation Plots with Enrichment

This repository contains the source code for an R package that generates weighted network models giving hierarchical reductions of a given network of features of a sample data set. For gene networks the included [Gephi](https://gephi.org) plugin can be used to explore the models and related Gene Ontology annotations.

  1. **[Build workflow](#BuildWorkflow)** (skip if you want to install from archive files)
  2. **[Installation](#Installation)**
  3. **[Dependency notes](#DependencyNotes)**

Build workflow <a name="BuildWorkflow"></a>
--------------

Open up one *R session* and one *terminal session*.

In R,
```
install.packages(c("devtools", "roxygen2")) # If necessary
library(devtools)
library(roxygen2)
document('grapeplots/') # Repeat this before each rebuild
```

You might already have all the dependencies installed, but if not go to [Dependency notes](#DependencyNotes) first. In the terminal session:
```
R CMD INSTALL grapeplots
R CMD build grapeplots
R CMD check grapeplots_1.0.0.tar.gz
```

or just
```
./build.sh
```

(The version number might vary.)

For the Gephi plugin, clone and set up the plugin-development repository from their [website](https://gephi.org). Choose plugin type 'Filter' when prompted. Under `src/main/java/` make directory `grapeplots/` and copy in all the `.java` files from `gephi_plugin_code`. Then run

```
mvn clean package
```

For the above you need Maven to be installed.

The `.nbm` file created somewhere in a `target` folder can be loaded directly into Gephi as a new plugin. You can select 'GRAPE Plot' from the 'Topology' category of filters.

Installation <a name="Installation"></a>
------------
Assuming you built/checked the R package:

```
install.packages("grapeplots_1.0.0.tar.gz", repos=NULL)
library(grapeplots)

help(generate_reduction) # To see usage of the main function
```

Dependency notes <a name="DependencyNotes"></a>
----------------

```
install.packages(c("igraph","emdist","mclust","pbmcapply"))
```

If you want support for the GCT file format for gene expression data:

```
install.packages('BiocManager')
BiocManager::install("CePa")
```

For the GO annotation functionality:

```
BiocManager::install("rols")
```


