# dasper 0.99.0

NEW FEATURES

* Converted `dasper` into a Bioconductor-friendly format using [biocthis](https://lcolladotor.github.io/biocthis/)
* Added `junction_load()`, which loads raw junction data outputted from RNA-sequencing alignment into an `RangedSummarizedExperiment` object. Includes an option to allow download of user-specified control junctions. 
* Added `junction_annot()`, which uses information from reference annotation and the strand of a junction to classify junctions as "annotated", "novel_acceptor", "novel_donor", "novel_exon_skip", "novel_combo", "ambig_gene" and "unannotated". 

