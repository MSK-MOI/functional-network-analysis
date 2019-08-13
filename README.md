<http://github.com/jimmymathews/grape-plots>

Gene Reduction/Annotation Plots with Enrichment
===============================================

This repository contains the source code for an R package that generates weighted network models giving hierarchical reductions of a given network of features of a sample data set. For gene networks the included [Gephi](https://gephi.org) plugin can be used to explore the models and related Gene Ontology annotations.

  1. **[Build workflow](#BuildWorkflow)** (skip if you want to install from archive files)
  2. **[Installation](#Installation)**
  3. **[Dependency notes](#DependencyNotes)**
  4. **[Example usage](#Example)**

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

For the GO annotation functionality you will need

```
BiocManager::install("rols")
```

and you will also need the `goa_human.gaf` annotation file available from the [EBI](https://www.ebi.ac.uk/GOA/downloads). For the most updated information, the script uses web API calls for the term definitions rather than a term definition file. However a partial cache file system is used to speed up repeated lookups.

Example Usage <a name="Example"></a>
-------------

Open up an R session or RStudio, then:

```
source('lung_gtex_run.R')
```

This example uses lung tissue RNA expression data from GTEx. The file `example_data/lung_tissue_expression_gtex_abridges.csv` is abridged to the 2000 genes with the most variance in the original dataset available from the [GTEx portal](https://gtexportal.org/) to save time. The main computations take about 7 minutes on an 8-core workstation, and the annotation lookup API calls take about 2 minutes. You should see output like the following:

```
Calculating weighted network reduction based on

node_data_file:             example_data/lung_tissue_expression_gtex_abridged.csv
topology_file:              NA
correlation_transformation: none
normalization:              none
method:                     gmt

[1/6] Loading data and loading/building network
      Inferring feature network from data using Pearson correlation. (No topology_file supplied).
      Data set 427x2000 all numeric.
      Calculated correlations (2000x2000)
      Cutoff value for correlation: 0.7

      2000 nodes processed of 2000 (Elapsed 0.04 mins of expected total 0.04 mins)
      Inferred 11361 edges, out of possible 1999000, with connectivity 0.00568334167083542
      Support of inferred graph contains 1475 nodes.
      Average degree 15.4047457627119.
      Diameter 16.

[2/6] Integrating sample set and network data
      Number of nodes in network not in the data set: 0
      Number of nodes in the data set not in the network: 525
      A2M ABTB1 ACADVL ACE ACVRL1 ...
      Number of nodes in common: 1475
      11361 edges after integrating dataset and feature network.

[3/6] Fitting Gaussian mixture models
      Number of cores according to parallel::detectCores(): 8
      Trying 7.
  |=====================================================================| 100%, Elapsed 02:02
      11361 2D models with 3 populations.

[4/6] Comparing all adjacent models
  |=====================================================================| 100%, Elapsed 01:31
      400710 adjacent edge pairs considered.

[5/6] Averaging over intermediating triangles to get virtual edges with weights
      400710 of 400710 edge pairs collated.
      53346 total virtual edges.

DONE
6.685659 minutes
```

For GO annotation:

```
gaf_annotate_graphml("example_data/lung_gtex_hierarchy.graphml")
gaf_annotate_graphml("example_data/lung_gtex_weighted.graphml")
```

Gephi hints:

  1. Install the plugin NBM file in Gephi with Tools > Plugins > Downloaded > Add Plugin. Find it under Filters > Topology > GRAPE Plot.
  2. Open the GraphML file `lung_gtex_weighted_annotated.graphml` or `lung_gtex_hierarchy_annotated.graphml`.
  3. Copy the `average_distance` edge attribute data to the `Weight` column, so that the weights calculated above, and not just the topology, influence the layout.
  4. Use the Force Directed Layout 2.
  5. Zoom out with the scroll wheel.
  6. For the hierarchy graph, set the node size to be a function of the `absorption_time` attribute.
  7. Turn on labels coming from the `name` attribute.

![alt text](example_data/lunggtex.png)











