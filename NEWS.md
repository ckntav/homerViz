# homerViz 0.99.0

* Initial Bioconductor submission.
* Parse HOMER `findMotifsGenome.pl` output directories into tidy data frames
  (`read_homer_output()`).
* Dot plots and scatter plots for known and de novo motif enrichment
  (`plot_known_dotplot()`, `plot_denovo_dotplot()`,
   `plot_known_scatter()`, `plot_denovo_scatter()`).
* Motif logo panels with enrichment statistics
  (`plot_known_logos()`, `plot_denovo_logos()`).
* Interactive HTML reports with embedded logos and sortable tables
  (`save_known_html()`, `save_denovo_html()`, `save_homer_html()`).
* Save plots to PDF, PNG or SVG with `saveViz()`, re-exported from gVenn.
