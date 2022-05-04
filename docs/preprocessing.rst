.. _preprocessing:

Preprocessing tutorial
======================

The pre-processing module of CanSig consists of multiple steps.
In the first step, we apply standard quality control measures.
The second step consists in CNV inference to remove cells annotated as malignant that do not
present CNVs and cells annotated as non-malignant that present CNVs.
In the optional, third step, known gene signatures can be scored.

.. todo:: Add graphic.

Loading the data
----------------------
To illustrate the usage of the pre-processing module of CanSig, we use three samples from
[Zhang2021]_.

.. code-block:: python
    from cansig.tutorial import load_data
    adatas, gene_order, scoring_dict = load_data()



The function load_data loads the following objects: adatas, a list of three AnnData
objects; gene_order, a DataFrame that contains
the gene position in the chromosome, and scoring_dict, a dictionary
containing names of known signatures as keys and their associated genes in a list a values.

.. important:: CanSig always expects raw counts (not log library size normalized) as input.

Each cell in the AnnData objects has to have a cell type annotation stored in
`adata.obs[celltype_column]`. The gene order file is expected to be a DataFrame that contains
as index the gene names and three columns: chromosome, start and end. The chromosome
needs to be formatted as "chr<number of the chromosome>".

.. code-block:: python
    print(gene_order.head(2))
                chromosome   start     end
    MIR1302-9.3        chr1   29554   31109
    FAM87B             chr1  752751  755214

Lastly, the key of the scoring_dict indicates into which column of
`adata.obs` the score will be saved and the value should be a list of genes used for scoring.

.. code-block:: python
    print(scoring_dict)
    {'Epi1': ['KRT16', 'CSTB', 'S100A2', 'S100A9', 'SPRR1B', ...], ...}


For scoring, we use scanpy's `score_genes` function with log library size normalized
counts as input. We apply the scoring function per batch and not on the combined
dataset. See `score_genes <https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html>`_
for more details.

Preprocessing
-------------
Before we can start with the preprocessing, we need to define which cell types are
considered as malignant or undetermined. The remaining cell types will be considered as
non-malignant. Furthermore, we need to define which cell types are used as reference for
InferCNV. These are typically non-malignant cells representing a variety of different
cell types. For more details, see `InferCNV <https://github.com/broadinstitute/inferCNV/wiki>`_.

.. code-block:: python
    malignant_celltypes = ["Epi"]
    undetermined_celltypes = ["Undetermined"]
    reference_celltypes = ["Pericytes",
                           "Bcells",
                           "FRC",
                           "Endothelial",
                           "Fibroblasts",
                           "Myeloid",
                           "Tcells"]

Now, we can run preprocessing by importing the function and defining the thresholds
used for quality control.

.. code-block:: python
    from cansig import preprocessing

    adata = preprocessing(adatas,
                         batch_id_column = 'batch_id',
                         celltype_column = 'cell_type',
                         malignant_celltypes=malignant_celltypes,
                         undetermined_celltypes=undetermined_celltypes,
                         reference_celltypes=reference_celltypes,
                         min_counts=1_500,
                         max_counts=50_000,
                         min_genes=700,
                         threshold_pct_mt_counts=30.,
                         gene_order=gene_order,
                         scoring_dict=scoring_dict,
                         figure_dir=None)


.. Note:: Instead of calling the function with the AnnData objects in memory, we can also
    provide a list of paths to .h5ad files. This can save memory if many
    samples are preprocessed. If the data is loaded from memory we have to define
    a column that contains the batch_id. If the data is loaded from disc and the
    batch_id_column is not already in `adata.obs` it will be set to the filename.

Cells with less than `min_counts` counts or more than `max_counts` counts will be filtered out.
Furthermore, cells with fewer than `min_genes` genes expressed or with a higher
percentage count in mitochondrial genes than `threshold_pct_mt_counts` will also be removed.
All plots generated during preprocessing will be stored in `figure_dir`. For more detail,
see preprocessing docs (TODO: add link here).

Outputs
--------
The function `preprocessing` returns a single AnnData object containing all the high quality cells
from the inputted samples.

.. note:: Since the goal of CanSig is to discover shared signatures, we do an inner join
    for the genes. This means only genes present in all samples will be kept in the
    final AnnData (This behavior can be changed by setting `join` to "outer".).


For each cell the following annotations are added in `adata.obs`:

- `n_counts`: The library size of the cell.
- `log_counts`: `log(n_counts)`.
- `n_genes`: The number of genes expressed in the cell.
- `pct_zero_genes`: `n_genes` divided by the number of all genes.
- `pct_counts_mt`: The counts corresponding to mitochondrial RNA divided by `total_counts`.
- `malignant_annotation`: Boolean indicating if the cell is considered malignant based on it cell type.
- `malignant_cnvs`: Boolean indicating if the cell is considered malignant based on its inferred CNV profile.
- `malignant`: Boolean indicating if the cell is considered malignant based on its celltype and CNV profile.

For more details on the malignant/non-malignant status annotation, see  the Methods section
of our `paper <https://www.biorxiv.org/content/10.1101/2022.04.14.488324v1>`_.

.. todo:: Do we want to add cell cycle scores? Problem: When different gene names are used?

.. important:: Rare malignant cells might be difficult to annotate. Therefore, we consider
    cells, that show CNVs but are annotated as undetermined, as malignant. However, cells
    that are annotated as non-malignant but show CNVs will not be considered as
    malignant cells.

In additions to the above annotations, a column for each element in the `scoring_dict` is
added `adata.obs`. For this tutorial, these are the known signature from [Zhang2021]_,

- `Mucosal`: The mucosal immunity-like (Mucosal) program was characterized by the expression of genes associated with innate immune response (e.g., S100P) and mucosal defensive mechanisms including mucosal chemokine (e.g., CXCL17) and mucus production (e.g., AGR2 and MUC20)
- `Stress`: The stress responses (Stress) program consisted of immediate early genes (e.g., EGR1, JUN, and FOS) that are activated in response to widespread cellular stimuli and displayed upregulation of TNFα signaling, UV response, p53, and apoptosis pathways
- `AP`: The antigen presentation (AP) program had increased expression of major histocompatibility complex (MHC) class II molecules (e.g., CD74, HLA-DPA1, and HLA-DRA/B1/B5) that are involved in initiating adaptive antitumor immune responses
- `Cycling`: The cell cycle (Cycling) program was characterized by high expression of genes involved in cell proliferation (e.g., CENPW, CKS1B, and BIRC5) and presented activation of the E2F targets, G2M checkpoint and MYC targets pathways, suggesting tumor cell proliferation
- `Epi1`: The Epi1 program was characterized by the expression of stress keratins (KRT6, KRT16, and KRT17) that are associated with keratinocyte hyperproliferation and therefore may play a role in enhancing tumorigenesis and tumor growth
- `Epi2`: The Epi2 program had the overexpressed genes related to the terminal differentiation such as envelope proteins (SPRR1A/1B) and calprotectin (S100A8/9), apical surface, the PI3K/AKT/mTOR signaling, the complement, and p53 pathways
- `Mes`: he mesenchymal cell-like properties (Mes) program consisted of genes such as VIM and SPARC and showed activation of epithelial-mesenchymal transition (EMT) and angiogenesis pathways.
- `Oxd`: Finally, the oxidative stress or detoxification (Oxd) program was characterized by the expression of multiple peroxidases and reductases (e.g., GPX2 and AKR1C1) involved in the defense against oxidative damage.

Scoring known signatures is an optional step and CanSig can function without. However,
by using known signatures one can assess the quality of the low dimensional representation
found in the next step by the model, in addition to using convergence metrics. We therefore
recommend to try using known signatures for the cancer type studied for better
interpretability of the results.


Furthermore, the CNV profile of each cell is stored in `adata.obsm["X_cnv"]`. Each row
in `adata.obsm["X_cnv"]` corresponds to a cell and each column represents a gene. The genes
are sorted by their position in the chromosome.

.. note:: Since the CNVs inferred by InferCNV are highly correlated we only store the CNVs
    for every 10th gene to save memory. The number of genes skipped is controlled by the
    `step` parameter in `preprocessing`.

In addition to the AnnData object, `preprocessing` also generates
plots for each sample to assess the quality of the data and the split into malignant and
non-malignant cells. The plots are stored in <figure_dir>/<batch_id>. The first plot is
created during the quality control step and gives insights into which cells are being
filtered out. This plot is saved to quality_control.png.

.. todo:: Add image for quality control
    Figure caption:  (A) Histogram of count depth per cell. (B) Histogram of number
    of genes detected per cell. (C) Count depth distribution. (D) Number of genes versus
    the count depth coloured by the fraction of mitochondrial reads. Mitochondrial read
    fractions are only high in particularly low count cells with few detected genes.
    Source: [Luecken2019]_

The next plot is generated after inferring CNVs. It shows the chromosome heatmap
separated into malignant and non-malignant and the malignant cells are further divided
into reference and non-reference Cells. The non-malignant cells should not show CNVs.
This plot is saved to chromosome_heatmap.png

.. todo:: Add image of the chromosome heatmap showing separation of malignant and
    non-malignant cells.

.. todo:: umap for each score + umap for malignant/non-malignant cells in CNV space.

.. note:: For faster pre-processing plotting can be turned off by setting plot to False.

.. todo:: Are there other useful plots that we want to add here?

References
----------

.. [Zhang2021] Zhang, X., Peng, L., Luo, Y. et al. Dissecting esophageal squamous-cell carcinoma ecosystem by single-cell transcriptomic analysis. Nat Commun 12, 5291 (2021). https://doi.org/10.1038/s41467-021-25539-x

.. [Luecken2019] Luecken, M. D., Theis, F. J. Current best practices in single-cell RNA-seq analysis: a tutorial. Molecular systems biology, 15(6), e8746 (2019). https://doi.org/10.15252/msb.20188746
