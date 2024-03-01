# RichStudio
R Shiny app for integrative enrichment analysis and visualization of multiple gene datasets

How to use the app:
1. Visit hosted version on http://hurlab.med.und.edu:3838/RichStudio
2. Clone the repository and running the following in R console:

```R
renv::restore()
devtools::load_all()
RichStudio::launch_RichStudio()
```

In the future, we hope to publish the app as a package.
   
## App Overview:
### Upload files
- Upload either DEG sets or enrichment results to the app
  - Supports file or text uploads
- Keep track of all uploaded/created files
  - Rename/remove files
  - View as searchable table
  - Export as .txt, .csv, or .tsv
 
### Enrichment
- Perform gene enrichment from DEG sets or gene lists
  - Supports GO, KEGG, and Reactome annotation databases
  - Supports BP, MF, CC ontologies
  - 20 different supported species *(human, mouse, pig, etc.)*

### Cluster
- Uses kappa clustering to group enriched terms by gene overlap
 
## Visualizations
### Enrichment results
- Bar plot
- Dot plot
- Network plot
- Heatmap

### Cluster results
- Comprehensive cluster heatmap
- Individual term heatmap
