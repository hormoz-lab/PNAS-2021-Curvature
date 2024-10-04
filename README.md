# Manifold Curvature Paper (PNAS 2021)

This repository provides instructions and code to reproduce all results, numerics and figures from the [Manifold Curvature paper](https://doi.org/10.1073/pnas.2100473118):

> D. Sritharan*, S. Wang*, S. Hormoz. "Computing the Riemannian curvature of image patch and single-cell RNA sequencing data manifolds using extrinsic differential geometry." PNAS (2021), https://doi.org/10.1073/pnas.2100473118.

The instructions will walk you through how to (i) download relevant code and data, (ii) run simulations, and (iii) regenerate figures and reproduce numerics from the paper.

For convenience, there are four main paths referred to below that appear in multiple steps:

	CODE_PATH	- directory where the ManifoldCurvature code is cloned
	PAPER_PATH	- directory where this repository is cloned	
	DATA_PATH	- directory where image and SC data is downloaded
	RESULTS_PATH	- output from simulations and custom scripts for manuscript preparation	

## Download Software

1. If you haven't done so, first [download and install the ManifoldCurvature code](https://github.com/hormoz-lab/ManifoldCurvature) according to the instructions, making sure to update the MATLAB path as instructed.

2. Download the code from this repository into the directory given by PAPER_PATH:

	```bash
	$ git clone https://github.com/hormoz-lab/PNAS-2021-Curvature.git ${PAPER_PATH}
	```

3. [Download Seurat V3.1.2](https://satijalab.org/seurat/install.html) for working with single-cell datasets.

4. [Download dynamo V0.96.0](https://dynamo-release.readthedocs.io/en/latest/ten_minutes_to_dynamo.html#how-to-install) for RNA velocity analysis.

## Download Data

5. [Download the van Hateren IML dataset](http://bethgelab.org/datasets/vanhateren) and unzip it into DATA_PATH/Images.

6. [Download the PBMC scRNAseq dataset](https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_NGSC3_DI_PBMC/Parent_NGSC3_DI_PBMC_filtered_feature_bc_matrix.tar.gz) and unzip it into DATA_PATH/PBMC.

7. [Download the gastrulation scRNAseq dataset](https://content.cruk.cam.ac.uk/jmlab/atlas_data.tar.gz) and unzip it into DATA_PATH/Gastrulation. Rename genes.tsv to features.tsv, raw_counts.mtx to matrix.mtx, and meta.tab to meta.dat.

8. [Download the brain scRNAseq dataset](https://cf.10xgenomics.com/samples/cell-exp/1.3.0/1M_neurons/1M_neurons_filtered_gene_bc_matrices_h5.h5) into DATA_PATH/Brain. [Retrieve](https://cf.10xgenomics.com/samples/cell-exp/1.3.0/1M_neurons/1M_neurons_analysis.tar.gz) the pre-computed PCA/t-SNE projections. Move analysis/pca/20_components/projections.csv to DATA_PATH/Brain/PCA.csv and analysis/tsne/2_components/projections.csv to DATA_PATH/Brain/tsne.csv. [Download](https://github.com/tinglabs/scAIDE/tree/master/Predicted%20labels/Neural%20dataset) files with cell annotations (cell_barcodes.zip, neural_1m_cluster_labels_assignment.RData and neural_64_cluster_prediction.csv) into DATA_PATH/Brain.

## Analysis for Paper

9. Temporarily add PAPER_PATH to MATLAB's search path:

	```MATLAB
	>> addpath(genpath(PAPER_PATH));
	```

10. Compute Laplace-Beltrami operator and eigenvalues for intrinsic approach:

    ```MATLAB
	>> generate_intrinsic_results(RESULTS_PATH);
    ```

11. Simulate toy models for extrinsic approach:

    ```MATLAB
	>> generate_extrinsic_results(RESULTS_PATH);
	>> generate_confounder_results(RESULTS_PATH);
    ```

12. Set up some data structures for computing Klein bottle embeddings:

    ```MATLAB
	>> setup_DCT_basis(RESULTS_PATH);
	>> setup_Carlsson_KB(RESULTS_PATH);	
    ```

13. Preprocess image patch dataset (original and augmented) and map points to Carlsson Klein bottle embedding:

    ```MATLAB
	>> preprocess_Hateren(DATA_PATH, RESULTS_PATH, false);
	>> preprocess_Hateren(DATA_PATH, RESULTS_PATH, true);	
	>> map_to_KB(RESULTS_PATH);
    ```

14. Find optimal Klein bottle embedding for image patch dataset:

	```MATLAB
	>> fit_optimal_KB(RESULTS_PATH);
	```

15. Compute curvature of image datasets and Klein bottle embeddings:

	```MATLAB
	>> generate_image_results(RESULTS_PATH);
	```

16. Subsample transcripts in PBMC dataset:

	```MATLAB
	>> subsample_transcripts('DATA_PATH/PBMC');
	```

17. Preprocess PBMC dataset and downsampled versions:

	```R
	> setwd(PAPER_PATH/SC)
    > source('preprocess_PBMC.R')
	> preprocess_PBMC(DATA_PATH/PBMC)
	> preprocess_PBMC(DATA_PATH/PBMC/SSTranscripts_Fact1)
	> preprocess_PBMC(DATA_PATH/PBMC/SSTranscripts_Fact2)
	```

18. Extract PCs from gastrulation dataset and get total transcripts per cell:

	```R
	> setwd(PAPER_PATH/SC)
    > source('preprocess_gastrulation.R')
	> preprocess_gastrulation(DATA_PATH/Gastrulation)
	```
	```MATLAB
	>> counts = load_10x_mtx(DATA_PATH/Gastrulation);
	>> dlmwrite('DATA_PATH/Gastrulation/TotCounts.csv', full(sum(counts,2)));	
    ```

19. Extract cell annotations for brain dataset and get total transcripts per cell:

	```R
	> setwd(DATA_PATH)
	> load('Brain/neural_1m_cluster_labels_assignment.RData')
	> write.table(neural_1m_cluster_labels_assignment, 'Brain/ClusterPhenotype.csv', sep=',', col.names=F)
	```
    ```MATLAB	
	>> counts = load_10x_h5_matrix(DATA_PATH/Brain/1M_neurons_filtered_gene_bc_matrices_h5.h5);
	>> dlmwrite('DATA_PATH/Brain/TotCounts.csv', full(sum(counts,2)));
    ```

20. Preprocess the dentate gyrus dataset and infer RNA velocity by running the following IPython notebook:

	```Python
	>>> PAPER_PATH/SC/DentateGyrus.ipynb
	```

21. Compute curvature for scRNAseq datasets:

    ```MATLAB
	>> generate_SC_results(DATA_PATH, RESULTS_PATH);
    ```

22. Perform statistical analysis for scRNAseq datasets:

	```MATLAB
	>> SC_err_analysis(RESULTS_PATH, 'PBMC');
	>> SC_err_analysis(RESULTS_PATH, 'Gastrulation');
	>> SC_err_analysis(RESULTS_PATH, 'Brain');
	>> SC_err_analysis(RESULTS_PATH, 'DentateGyrus');
	```

23. Fit linear models for scRNAseq datasets:

	```MATLAB
	>> fit_curvature_to_gene_exp(DATA_PATH, RESULTS_PATH, 'PBMC');
	>> fit_curvature_to_gene_exp(DATA_PATH, RESULTS_PATH, 'DentateGyrus');
	```

## Generate Outputs

24. Prepare main and supplementary figures:

    ```MATLAB
	>> make_intro_subplots(RESULTS_PATH);				% Figure 1
	>> make_intrinsic_subplots(RESULTS_PATH);			% Figure S1
	>> make_extrinsic_subplots(RESULTS_PATH);			% Figure 2
	>> make_confounder_subplots(RESULTS_PATH);			% Figure S2
	>> make_image_subplots(RESULTS_PATH);				% Figures 3, S3
	>> make_SC_subplots(DATA_PATH, RESULTS_PATH);			% Figures 4-5, S4-S8
    ```

25. Compute all the statistics and numerical quantities found in the manuscript text:

    ```MATLAB
	>> generate_manuscript_stats(DATA_PATH, RESULTS_PATH);
    ```

#### Prepared By: Duluxan Sritharan
