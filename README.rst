.. |date| date::

*************************
OpenPedCan-exprs-shinyapp
*************************

:authors: Komal S Rathi
:contact: rathik@email.chop.edu
:organization: D3B, CHOP
:status: This is "work in progress"
:date: |date|

Introduction
============

Mock app built using R shiny to create Tumor vs Normal and Tumors-only plots to visualize TPM expression.

Structure
=========

.. code-block:: bash

    .
	├── R
	│   ├── pubTheme.R
	│   ├── tumor_plot.R
	│   ├── tumor_vs_normal_plot.R
	│   └── viewDataTable.R
	├── README.rst
	├── code
	│   └── subset_data.R 
	├── data
	│   ├── gene-expression-rsem-tpm-collapsed-subset.rds 
	│   └── histologies.tsv
	├── server.R
	├── ui.R
	└── www


Functions:
==========

1. **code/subset_data.R**: This script is to subset the original expression matrix to a smaller subset for mock application.
2. **R/viewDataTable.R**: `DT::datatable` with various add-on, buttons etc
3. **R/pubTheme.R**: ggplot2 publication quality themes
4. **R/tumor_plot.R**: function for tumors-only plots
5. **R/tumor_vs_normal_plot.R**: function for tumor vs normal plots

Input:
======

Input files used are:

1. Expression matrix: `input/gene-expression-rsem-tpm-collapsed-subset.rds` which is a subset of gene-expression-rsem-tpm-collapsed.rds created using `code/subset_data.R`
2. Histologies file: `input/histologies.tsv`

Output:
=======

Two types of output created using this app: one is dynamic (which can be seen on the app as `plotly` plots and `DT` tables) and the corresponding static version is written to the `www/` folder as `.pdf` and `.tsv` files.

Note: data/ (input directory) and www/ (output directory) have been added to .gitignore.
