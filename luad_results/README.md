I was able to reproduce the results highlighting the PSG genes in a module in the context of lung cancer.

The steps were:

1. Using cBioPortal I downloaded the large zip file of all data associated with the TCGA LUAD (lung adenocarcinoma) PanCancer Atlas.
2. The file in that zip archive called `data_RNA_Seq_v2_mRNA_median_Zscores.txt` has expression data for about 500 samples and 20000 genes. I converted it to a `csv` file in the format of the other files in the `example_data` folder from this repository. That is, I transposed it, converted the tabs to commas, removed the extra of 2 columns describing the feature/gene set, removed NA columns, so that the first part looks like this:
```
,UBE2Q2P2,LOC155060,SSX9,SDR16C6P
TCGA.05.4244.01,-0.428,0.9965,-0.1387,-0.3334
TCGA.05.4249.01,-0.3333,-0.3213,-0.1387,-0.3334
TCGA.05.4250.01,-0.8573,-0.9371,-0.1387,-0.3334
```
3. In R, I calculated the Pearson correlation matrix, and ranked the genes by average absolute value of Pearson correlation with other genes.
4. Feature selection is done by taking top 3000 with respect to the above average.
5. The feature matrix with 3000 features is input into the GMT pipeline by `luad_run.R`, using a Pearson correlation cutoff of 0.55 .
6. The verbose output describing the GMT process being carried out is saved to `run_stdout.txt`.
7. I opened the resulting graph, `luad_weighted_3000.graphml`, in Gephi. The GMT edge weights contribute to the layout I used, as well as playing around with the "Feature Network Reduction" Gephi plugin (the source code provided as part of this repository). The result is saved to `tcga_luad_gmt.gephi`.

Note: To decide if the GMT has really contributed to this result, try to reproduce the result by using Gephi or some other analysis method on the topology of `luad_weighted_3000.graphml` alone (i.e. ignoring the GMT edge weights). It's possible that the correlation-based feature selection and network inference has already done much of the work here.

-Jimmy

November 22 2020
