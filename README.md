# monocle3_shiny
Shiny app to perform trajectory analysis on scRNAseq data

Monocle3 documentation can be found [here](https://cole-trapnell-lab.github.io/monocle3/)

Requires Rdata or Rds file containing Seurat object

Available at: https://saezlab.shinyapps.io/monocle3_shiny/

## Instructions

Simple to use shiny application for trajectory analysis from Monocle3.

1. Upload your Rdata or Rds file containing Seurat object by clicking browse or drag and drop.
2. Once file upload bar says "Upload Complete", hit "Upload" to push this data into the app.
3. Select an assay, identity and progenitor gene to work with, then hit the "Run Monocle3" button.
4. If you want to download any plots, select them from the dropdown menu and hit "Download Plot".
5. Download your results in csv form by hitting the "Download Results" button.
