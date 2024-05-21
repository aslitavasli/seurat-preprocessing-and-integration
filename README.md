# Integrative Analysis and Visualization of Single-Cell RNA sequencing (scRNA-seq) Data of Tonsil Samples via Seurat Objects


### By Amoli Agarwal and Asli Tavasli



 **Note:**
Due to the extensive size of data, we are unable to upload it on GitHub. For inquiries or access to raw data,
email at: asli.tavasli.27@dartmouth.edu or amoli.agarwal.27@dartmouth.edu


This is the GitHub repo for our project titled "Integrative Analysis and Visualization of Single-Cell RNA sequencing
 (scRNA-seq) Data of Tonsil Samples via Seurat Objects" for the Women In Science Project.
 
 B and T lymphocytes are crucial to the acquired immune response because they are the only cells capable of specifically recognizing and responding to antigenic epitope, a part of an antigen recognized by the immune system. The reactivity of these lymphocytes are determined by their heritable antigenic receptors. 

Seurat objects facilitate the identification and characterization of different subsets of B and T cells, including naive, memory, and effector cells, by enabling precise clustering based on gene expression profiles. This detailed categorization is essential for understanding the roles these cells play in immune responses and how they change in different contexts, such as infection, vaccination, or autoimmune conditions. By integrating multiple datasets, we can compare B and T cell responses across different conditions or diseases, identifying conserved and unique features of the immune response. This comparative analysis is crucial for developing targeted therapies and vaccines.

Data: We used data from tonsils, which are rich in lymphoid tissue, containing both germinal centers and t-cell zones which provide a diverse immune repertoire to analyze. Data was obtained from a study on single-cell analysis of human B cell maturation.

To visualize the high-dimensional clustering data, we used a Uniform Manifold Approximation and Projection (UMAP), a dimensionality reduction technique. We first clustered data of three separate datasets, then integrated the datasets together to boost statistical power and facilitate accurate comparative analysis. 

In our integrated UMAP, we set our resolution such that seven distinct clusters were formed, separated by identity. Cluster marker analysis allows us to identify proteins that differentiate clusters from one another, leading to potential new findings in proteins found in specific cell types. 

