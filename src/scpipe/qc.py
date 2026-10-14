"""Quality control functions for scRNA-seq data."""

import numpy as np
from scipy.stats import median_abs_deviation
import scanpy as sc

def compute_qc_metrics(adata, mt_prefix = "MT-"):
    adata.var['MT-'] = adata.var_names.str.startswith(mt_prefix) 
    sc.pp.calculate_qc_metrics(adata, qc_vars=[mt_prefix], inplace= True )

def is_outlier(adata,metric , n_mads,upper_only=False):

    """
    Flag cells whose value for a QC metric is far from the median.

    Parameters: 
    adata : AnnData
        Annotated data matrix with QC metrics in adata.obs.
    metric : str
        Column of adata.obs to test (e.g. "log1p_total_counts").
    n_mads : float
        Number of median absolute deviations defining an outlier.
    upper_only : bool
        If True, only flag values above the median.
    
    Returns:
    pandas.Series of bool
        True for outlier cells.
    """

    MAD = median_abs_deviation(adata.obs[metric])
    median = np.median(adata.obs[metric])
    thresh_pos   = median + n_mads * MAD
    
    if upper_only == True:
        mask = adata.obs[metric] > thresh_pos
        
        
    else:
        thresh_under = median - n_mads * MAD
        mask = (adata.obs[metric] > thresh_pos) | (adata.obs[metric] < thresh_under)
        
    return mask

def flag_low_quality_cells(adata,n_mads,n_mads_mt): 
    adata.obs['outlier_n_genes_by_counts'] = is_outlier(adata, "log1p_n_genes_by_counts", n_mads)
    adata.obs['outlier_total_counts']      = is_outlier(adata, "log1p_total_counts", n_mads)
    adata.obs['outlier_pct_counts_MT-']    = is_outlier(adata, "pct_counts_MT-", n_mads_mt, upper_only=True)    

    adata.obs["qc_fail"] = (adata.obs['outlier_n_genes_by_counts']
                        | adata.obs['outlier_total_counts']
                        | adata.obs['outlier_pct_counts_MT-'])
    return adata.obs["qc_fail"]


def filter_cells_and_genes(adata, qc_fail, min_cells=3):
    """Remove cells that fail QC, then genes detected in too few cells.

    Parameters
    ----------
    adata : AnnData
    qc_fail : pandas.Series of bool
        True for cells to remove (output of flag_low_quality_cells).
    min_cells : int
        Minimum number of cells in which a gene must be detected.

    Returns
    -------
    AnnData
        A filtered copy of adata.
    """

    adata_clean = adata[~adata.obs["qc_fail"]].copy() 
    return adata_clean 
    



