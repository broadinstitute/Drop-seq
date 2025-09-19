"""
kmeans_effect_clustering
========================

This module provides utilities for performing K-means clustering on gene expression
effect matrices and generating visual outputs for analysis. It integrates with
AnnData and Scanpy for preprocessing, clustering, and visualization, and supports
end-to-end workflows for clustering genes across a range of k values.

Main functionality
------------------
- Input parsing:
  Convert tab-delimited matrices into AnnData objects for downstream processing.

- Cell type ordering:
  Use PCA and dendrogram construction to reorder cell types and visualize
  similarity across cell types with heatmaps.

- K-selection:
  Run K-means across a range of k values, compute silhouette scores and within-cluster
  sum of squares (WCSS), and generate diagnostic plots to guide cluster number choice.

- Gene clustering:
  Cluster genes at a specified k, generate clustered heatmaps, and export cluster labels.

- Output management:
  Create standardized file paths, save figures, merge label files into wide matrices,
  and optionally export combined multi-page PDFs.

Workflow
--------
The high-level `run_gene_clustering_workflow` function orchestrates the complete
process:
1. Parse the input matrix.
2. Order cell types and plot a dendrogram heatmap.
3. Perform k-selection analysis and export metrics.
4. Cluster genes for one or multiple k values.
5. Merge labels into a single matrix.
6. Optionally save all figures into a PDF.

Dependencies
------------
- Python standard library: os
- Data handling: pandas
- Visualization: matplotlib, scanpy
- Data structure: anndata
- Machine learning: scikit-learn (KMeans, silhouette_score)

Intended use
------------
This module is designed for command-line workflows and can be wrapped with
a CLI for reproducible clustering pipelines.
"""

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
from typing import Optional, Tuple
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

from matplotlib.backends.backend_pdf import PdfPages
from typing import Optional, Dict, List
import argparse
import logging

import matplotlib as mpl
mpl.rcParams["figure.max_open_warning"] = 0

def validate_save_path(save):
    """
    Validates that a save path is a string and points to a writable directory.

    Parameters
    ----------
    save : str or None
        The file path to validate. If None, no validation is performed.

    Raises
    ------
    ValueError
        If `save` is not None or str, or if the directory does not exist.
    PermissionError
        If the directory exists but is not writable.
    """
    if save is not None:
        if not isinstance(save, str):
            raise ValueError(f"`save` must be None or str, got {type(save)}")
        save_dir = os.path.dirname(save) or "."
        if not os.path.isdir(save_dir):
            raise ValueError(f"Directory does not exist or is not writable: {save_dir}")
        if not os.access(save_dir, os.W_OK):
            raise PermissionError(f"Directory is not writable: {save_dir}")


def parse_matrix_input(in_file: str) -> ad.AnnData:
    """
    Parses a tab-delimited input matrix file into an AnnData object.

    Parameters
    ----------
    in_file : str
        Path to the input file containing gene scores.

    Returns
    -------
    ad.AnnData
        AnnData object with cells as observations and genes as variables.
    """
    gene_scores_merged = pd.read_csv(in_file, sep='\t').transpose()
    # make an anndata object out of dataframe to take advantage of scanpy functions
    adata = ad.AnnData(gene_scores_merged)
    adata.obs['celltype'] = adata.obs_names.astype('category') # adding cell type label to obs
    return adata


def order_celltypes(
    adata: ad.AnnData,
    num_pcs: int = 10,
    figsize: Tuple[float, float] = (10, 5),
    save: Optional[str] = None,
) -> tuple[ad.AnnData, plt.Figure]:
    """
    Runs PCA, constructs a dendrogram, plots a heatmap, and reorders cell types.

    Parameters
    ----------
    adata : ad.AnnData
        Input AnnData object with cell types in `.obs`.
    num_pcs : int, optional
        Number of principal components to use. Default is 10.
    figsize : tuple, optional
        Size of the heatmap figure. Default is (10, 5).
    save : str or None, optional
        File path to save the heatmap. Default is None.

    Returns
    -------
    ad.AnnData
        Reordered AnnData object with cells arranged by dendrogram order.
    matplotlib.figure.Figure
        Heatmap figure object.

    Raises
    ------
    ValueError
        If the save path is invalid.
    """

    # validate save parameter
    validate_save_path(save)

    sc.pp.pca(adata, n_comps=num_pcs)
    sc.tl.dendrogram(adata, groupby="celltype", n_pcs=num_pcs)

    out = sc.pl.heatmap(
        adata,
        var_names=adata.var_names,
        groupby="celltype",
        dendrogram=True,
        cmap="seismic",
        vcenter=0,
        figsize=figsize,
        show=False,
        show_gene_labels=False,
    )
    fig = out["heatmap_ax"].figure

    if save is not None:
        fig.savefig(save, dpi=300, bbox_inches="tight")

    celltype_order = adata.uns["dendrogram_celltype"]["categories_ordered"]
    adata = adata[celltype_order, :].copy()

    return adata, fig

def k_selection_plot(
    adata: ad.AnnData,
    k_range: tuple[int, int] = (2, 30),
    ax: Optional[plt.Axes] = None,
    figsize: Tuple[float, float] = (10, 5),
    save: Optional[str] = None,
) -> Tuple[plt.Axes, pd.DataFrame]:
    """
    Runs K-means clustering for a range of k values, computes silhouette scores and WCSS,
    and plots them on dual axes.

    Parameters
    ----------
    adata : ad.AnnData
        Input AnnData object with features to cluster.
    k_range : tuple of int, optional
        Range of k values (inclusive). Default is (2, 30).
    ax : matplotlib.axes.Axes, optional
        Axis object to plot on. If None, a new one is created.
    figsize : tuple, optional
        Figure size if a new figure is created. Default is (10, 5).
    save : str or None, optional
        Path to save the plot. Default is None.

    Returns
    -------
    matplotlib.axes.Axes
        Axis containing the plot.
    pandas.DataFrame
        DataFrame with k, silhouette_score, and wcss values.

    Raises
    ------
    ValueError
        If the save path is invalid.
    """

    # validate save parameter
    validate_save_path(save)

    if ax is None:
        with plt.ioff():
            fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    adata_tmp = adata.copy().transpose()
    ks = list(range(k_range[0], k_range[1] + 1))

    wcss = []
    silhouette_scores = []
    for k in ks:
        kmeans = KMeans(n_clusters=k, random_state=12, n_init=10).fit(adata_tmp.X)
        labels = kmeans.labels_
        silhouette_scores.append(silhouette_score(adata_tmp.X, labels))
        wcss.append(kmeans.inertia_)

    df = pd.DataFrame({
        "k": ks,
        "silhouette_score": silhouette_scores,
        "wcss": wcss,
    })

    # Use pandas plotting with secondary_y
    df.plot(
        x="k",
        y="silhouette_score",
        ax=ax,
        color="black",
        marker="o",
        linestyle="--",
        label="Silhouette Score",
    )
    df.plot(
        x="k",
        y="wcss",
        ax=ax,
        secondary_y=True,
        color="blue",
        marker="o",
        linestyle="--",
        label="WCSS",
    )

    # Customize labels
    ax.set_xlabel("Number of Clusters")
    ax.set_ylabel("Silhouette Score", color="black")
    ax.right_ax.set_ylabel("Within-Cluster Sum of Squares (WCSS)", color="blue")

    if save is not None:
        ax.figure.savefig(save, dpi=300)

    return ax, df

def plot_gene_clustering(adata: ad.AnnData,
                         k:int = 20,
                         figsize: Tuple[float, float] = (15, 9),
                         save: Optional[str] = None,
                         out_cluster_labels_filename: Optional[str] = None,):
    """
    Performs K-means clustering of genes, generates a heatmap, and optionally saves labels.

    Parameters
    ----------
    adata : ad.AnnData
        Input AnnData object with genes as features.
    k : int, optional
        Number of clusters. Default is 20.
    figsize : tuple, optional
        Figure size for the heatmap. Default is (15, 9).
    save : str or None, optional
        Path to save the heatmap. Default is None.
    out_cluster_labels_filename : str or None, optional
        Path to save cluster labels as a CSV. Default is None.

    Returns
    -------
    matplotlib.figure.Figure
        Figure object containing the gene clustering heatmap.

    Raises
    ------
    ValueError
        If the save path is invalid.
    """
    # validate save parameter
    validate_save_path(save)

    adata_tmp = adata.copy().transpose()

    kmeans = KMeans(n_clusters=k, random_state=12, n_init=10).fit(adata_tmp.X)
    adata_tmp.obs['gene_clusters'] = kmeans.labels_ # retrieve the labels and add them as a metadata column in our AnnData object
    adata_tmp.obs['gene_clusters'] = adata_tmp.obs['gene_clusters'].astype("category")

    sc.tl.dendrogram(adata_tmp, groupby="gene_clusters",
                     cor_method="spearman", linkage_method="centroid",
                     use_rep="X")

    with plt.ioff():
        out = sc.pl.heatmap(
            adata_tmp,
            var_names=adata_tmp.var_names,
            groupby="gene_clusters",
            dendrogram=True,
            cmap="seismic",
            vcenter=0,
            swap_axes=True,
            figsize=figsize,
            show=False,
            show_gene_labels=True,
        )

    fig = out["heatmap_ax"].figure

    fig.text(
        0.03, 0.9, f"K={k}",
        ha="left", va="top",
        fontsize=20, fontweight="bold"
    )

    if save is not None:
        fig.savefig(save, dpi=300, bbox_inches="tight")

    if out_cluster_labels_filename is not None:
        adata_tmp.obs['gene_clusters'].to_csv(out_cluster_labels_filename)

    return fig

def make_base_paths(out_dir: str, prefix: str) -> dict[str, str]:
    """
    Generates output file paths for clustering results independent of k.

    Parameters
    ----------
    out_dir : str
        Directory for output files.
    prefix : str
        Prefix for output file names.

    Returns
    -------
    dict of str
        Dictionary of standard output paths.
    """
    return {
        "cell_ordering_heatmap": f'{out_dir}/{prefix}.celltype_order_heatmap.png',
        "k_selection_metrics": f'{out_dir}/{prefix}.k_selection_metrics.txt',
        "k_selection_plot": f'{out_dir}/{prefix}.k_selection_plot.png',
        "labels_matrix": os.path.join(out_dir, f"{prefix}.gene_cluster_labels_matrix.txt"),
        "out_pdf": f'{out_dir}/{prefix}.pdf',
    }

def make_cluster_paths (out_dir, prefix, k):
    """
    Generates output file paths for clustering results specific to a given k.

    Parameters
    ----------
    out_dir : str
        Directory for output files.
    prefix : str
        Prefix for output file names.
    k : int
        Number of clusters.

    Returns
    -------
    dict of str
        Dictionary with paths for heatmap and cluster label outputs.
    """
    return {
        "heatmap": os.path.join(out_dir, f"{prefix}.KmeansClusteredGenes_k{k}.png"),
        "labels": os.path.join(out_dir, f"{prefix}.KmeansClusteredGenes_k{k}_labels.txt"),
    }

def _read_labels(path: str) -> pd.DataFrame:
    """
    Reads a label file with two columns (gene, cluster).

    .. note::
        Internal use only. This function is used by `_merge_label_files`
        and `run_gene_clustering_workflow`.

    Parameters
    ----------
    path : str
        Path to the CSV file containing labels.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns ["gene", "cluster"].
    """
    df=pd.read_csv(path, sep=",")
    df.columns = ["gene", "cluster"]
    return df

def _merge_label_files(label_paths: Dict[int, str], out_path: str, delete: bool = True) -> str:
    """
    Merges per-k label files into a single wide matrix.

    .. note::
        Internal use only. Called from `run_gene_clustering_workflow`
        to consolidate clustering results.

    Parameters
    ----------
    label_paths : dict of int to str
        Dictionary mapping k values to label file paths.
    out_path : str
        Output path for the merged labels file.
    delete : bool, optional
        If True, deletes the original label files. Default is True.

    Returns
    -------
    str
        Path to the merged labels file.
    """
    merged: Optional[pd.DataFrame] = None
    for K in sorted(label_paths.keys()):
        df = _read_labels(label_paths[K]).set_index("gene")
        df = df.rename(columns={"cluster": f"k{K}"})
        merged = df if merged is None else merged.join(df, how="outer")
    merged = merged.sort_index()
    merged.to_csv(out_path, sep="\t", index=True)
    if delete:
        for p in label_paths.values():
            try:
                os.remove(p)
            except OSError:
                pass
    return out_path


def run_gene_clustering_workflow(
    in_file: str,
    out_dir: str,
    prefix: str,
    k_range: Tuple[int, int] = (2, 30),
    num_pcs: int = 10,
    clust_k: Optional[int] = None,
    gene_cluster_figsize: Tuple[float, float] = (18, 12),
    cell_order_heatmap_size: Tuple[float, float] = (15, 10),
    k_selection_figsize: Tuple[float, float] = (8, 8),
    labels_matrix: Optional[str] = None,
    out_pdf: Optional[str] = None,
    print_pngs: bool = True,
    delete_label_files: bool = True,
) -> None:
    """
    End-to-end workflow:
      - parse input
      - order cell types and plot heatmap
      - K-selection plot + metrics TSV
      - gene clustering heatmap(s) + labels
      - merged labels file
      - optional multi-page PDF

    Parameters
    ----------
    in_file : str
        Path to the input matrix file.
    out_dir : str
        Directory for saving outputs.
    prefix : str
        Prefix for output files.
    k_range : tuple of int, optional
        Range of k values for clustering. Default is (2, 30).
    num_pcs : int, optional
        Number of principal components to use. Default is 10.
    clust_k : int or None, optional
        If set, run clustering for only this k. Default is None.
    gene_cluster_figsize : tuple, optional
        Figure size for gene clustering heatmaps. Default is (18, 12).
    cell_order_heatmap_size : tuple, optional
        Figure size for cell ordering heatmap. Default is (15, 10).
    k_selection_figsize : tuple, optional
        Figure size for k-selection plot. Default is (8, 8).
    labels_matrix : str or None, optional
        Path to save merged labels file. Default is None.
    out_pdf : str or None, optional
        Path to save combined PDF. Default is None.
    print_pngs : bool, optional
        If True, saves individual PNG files. Default is True.
    delete_label_files : bool, optional
        If True, deletes intermediate label files. Default is True.

    Returns
    -------
        None

    Raises
    ------
    ValueError
        If output paths are invalid.
    """

    if labels_matrix is not None:
        validate_save_path(labels_matrix)
    if out_pdf is not None:
        validate_save_path(out_pdf)

    adata = parse_matrix_input(in_file)
    logging.info("Input matrix shape: %d genes × %d conditions",
                 adata.n_vars, adata.n_obs)

    base_paths = make_base_paths(out_dir, prefix)  # for metrics/plots

    figs: List[plt.Figure] = []
    label_paths: Dict[int, str] = {}

    # 1) Celltype ordering
    adata, fig1 = order_celltypes(
        adata,
        num_pcs=num_pcs,
        figsize=cell_order_heatmap_size,
        save=(base_paths["cell_ordering_heatmap"] if print_pngs else None),
    )
    figs.append(fig1)

    # 2) K-selection
    logging.info("Starting K-selection")
    ax2, df = k_selection_plot(adata, k_range=k_range, figsize=k_selection_figsize)
    fig2 = ax2.figure
    figs.append(fig2)
    if print_pngs:
        fig2.savefig(base_paths["k_selection_plot"], dpi=300, bbox_inches="tight")
    df.to_csv(base_paths["k_selection_metrics"], sep="\t", index=False)

    # 3) Gene clustering
    logging.info("Generating K-means clustering plots")
    def _run_gene_clust(K: int):
        clust_paths = make_cluster_paths(out_dir, prefix, K)
        fig = plot_gene_clustering(
            adata,
            k=K,
            figsize=gene_cluster_figsize,
            save=(clust_paths["heatmap"] if print_pngs else None),
            out_cluster_labels_filename=clust_paths["labels"],
        )
        figs.append(fig)
        label_paths[K] = clust_paths["labels"]

    if clust_k is None:
        for K in range(k_range[0], k_range[1] + 1):
            _run_gene_clust(K)
    else:
        _run_gene_clust(clust_k)

    # 4) Merge labels into matrix
    if labels_matrix is not None:
        _merge_label_files(label_paths, labels_matrix, delete=delete_label_files)

    # 5) Optional combined PDF
    if out_pdf is not None:
        with PdfPages(out_pdf) as pdf:
            for f in figs:
                pdf.savefig(f, bbox_inches="tight")

    for f in figs:
        plt.close(f)

    logging.info("Workflow Finished")

######################
# Command-line interface
######################

def build_arg_parser() -> argparse.ArgumentParser:
    """
    Build the command-line argument parser for the clustering workflow.

    Returns
    -------
    argparse.ArgumentParser
        Configured parser for this module’s CLI.
    """
    p = argparse.ArgumentParser(
        prog="kmeans_effect_clustering",
        description="Run gene clustering workflow from a tab-delimited effect matrix."
    )

    # required
    p.add_argument("--in-file", required=True, help="Path to tab-delimited matrix.")
    p.add_argument("--out-dir", required=True, help="Directory for outputs.")
    p.add_argument("--prefix", required=True, help="Filename prefix for outputs.")

    # reusable parser for pairs
    def _pair(s, cast):
        try:
            a, b = s.split(",")
            return cast(a), cast(b)
        except Exception:
            raise argparse.ArgumentTypeError(f"Expected 'x,y' with {cast.__name__} values.")

    # optional
    p.add_argument("--k-range", default=(2, 30),
                   type=lambda s: _pair(s, int),
                   help="Inclusive k range as 'low,high'. Default: %(default)s")
    p.add_argument("--num-pcs", default=10, type=int,
                   help="Number of PCs for dendrogram. Default: %(default)s")
    p.add_argument("--clust-k", default=None, type=int,
                   help="If set, run clustering only for this k. Default: %(default)s")
    p.add_argument("--gene-cluster-figsize", default=(18, 12),
                   type=lambda s: _pair(s, float),
                   help="Heatmap figsize 'w,h'. Default: %(default)s")
    p.add_argument("--cell-order-heatmap-size", default=(15, 10),
                   type=lambda s: _pair(s, float),
                   help="Cell ordering heatmap figsize 'w,h'. Default: %(default)s")
    p.add_argument("--k-selection-figsize", default=(8, 8),
                   type=lambda s: _pair(s, float),
                   help="K-selection plot figsize 'w,h'. Default: %(default)s")
    p.add_argument("--labels-matrix", default=None,
                   help="Output path for merged labels matrix TSV. Default: %(default)s")
    p.add_argument("--out-pdf", default=None,
                   help="Optional multi-page PDF path. Default: %(default)s")
    p.add_argument("--pngs", action="store_true",
                   help="Write per-figure PNGs (off by default).")
    p.add_argument("--log-level", default="INFO",
                   choices=["CRITICAL","ERROR","WARNING","INFO","DEBUG"],
                   help="Logging level. Default: %(default)s")
    return p

def main(argv: Optional[list[str]] = None) -> None:
    """
    Entry point for the CLI. Parses arguments and runs the workflow.
    """
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(message)s"
    )

    run_gene_clustering_workflow(
        in_file=args.in_file,
        out_dir=args.out_dir,
        prefix=args.prefix,
        k_range=args.k_range,
        num_pcs=args.num_pcs,
        clust_k=args.clust_k,
        gene_cluster_figsize=args.gene_cluster_figsize,
        cell_order_heatmap_size=args.cell_order_heatmap_size,
        k_selection_figsize=args.k_selection_figsize,
        labels_matrix=args.labels_matrix,
        out_pdf=args.out_pdf,
        print_pngs=args.pngs,
    )

if __name__ == "__main__":
    sys.exit(main())
