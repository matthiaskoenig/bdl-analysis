# Experimental data
The original data sets provided by the experimental partners are available in the folders
* antibodies (Antibodies for 3 samples per timepoint)
* bile_infarcts (Bile infarcts)
* fluidigm (Gene expression data based on Fluidigm platform)
* histology (Histological measurements)

The data was cleaned and combined in `raw_data.xlsx`. From this file the `.csv` are generated which are used for the `R` analysis.

The Fluidigm probes were mapped manually based on the gene identifiers to UniProt ids `probe_mapping.xlsx` and `probe_mapping.csv`.

Available `*.Rdata` files are data dumps of certain parts of the analysis.
