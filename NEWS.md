# dasper 1.3.8

NEW FEATURES

* Fix bug related to `pkgdown::deploy_to_branch()`. 

# dasper 1.3.6

NEW FEATURES

* Add functionality to annotate junctions with `gene_name`/`symbol` using `EnsDb` inputted into `junction_annot()`. 

# dasper 1.3.2

NEW FEATURES

* Fix bugs within `plot_sashimi()` and enable the visualization of raw junction counts. 

# dasper 1.1.1

NEW FEATURES

* Use of testthat edition 3 and parrallel running of tests. 

# dasper 0.99.2

NEW FEATURES

* Merge documentation into one man page for junction, coverage and outlier processing functions to reduce runtime of roxygen examples. 

# dasper 0.99.1

NEW FEATURES

* Change `outlier_detect()` to using `basilisk` for interfacing into python replacing `reticulate`.

# dasper 0.99.0

NEW FEATURES

* Converted `dasper` into a Bioconductor-friendly format using [biocthis](https://lcolladotor.github.io/biocthis/).
* Added `junction_load()`, which loads raw junction data from RNA-sequencing into an `RangedSummarizedExperiment` object. Includes an option to allow download of user-specified control junctions.
* Added `junction_annot()`, which uses information from reference annotation and the strand of a junction to classify junctions as "annotated", "novel_acceptor", "novel_donor", "novel_exon_skip", "novel_combo", "ambig_gene" and "unannotated".
* Added `junction_filter()`, which filters junctions by their count, width, annotation or if they overlap a set of user-defined regions.
* Added `junction_norm()`, which normalises raw junction counts (into a proportion-spliced-in) by dividing the counts of each junction by the total number of counts in it's associated cluster.
* Added `junction_process()`, a wrapper function for all "junction_" prefixed functions except `junction_load()`. 
* Added `junction_score()`, which scores patient junctions based on the extent their counts deviate from a control count distribution.
* Added `coverage_norm()`, which will load and normalise coverage for exonic/intronic regions corresponding to each junction.
* Added `coverage_score()`, which scores coverage associated with each junction based on it's deviation from control coverage distributions. 
* Added `coverage_process()`, a wrapper function for all "coverage_" prefixed functions. 
* Added `outlier_detect()`, which uses the junction scores and coverage scores as input into an unsupervised outlier detection algorithm to find the most outlier-looking junctions in each sample. 
* Added `outlier_aggregate()`, which aggregates the junction-level outlier data to a cluster-level. 
* Added `outlier_process()`, a wrapper function for all "outlier_" prefixed functions. 
* Added `plot_sashimi()`, which enables the visualisation of junction data across genes/transcripts or regions of interest.

