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
