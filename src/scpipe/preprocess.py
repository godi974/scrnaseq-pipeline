"""
Preprocessing:doublet removal, normalization and feature selection
"""
import scanpy as sc 




def remove_doublets(adata, random_state=0):
    
    """
    Detect doublets with Scrublet and remove them.

    Parameters
    ----------
    adata : AnnData
        QC-filtered data with raw counts in adata.X.

    returns
    -------
    Copy of adata without predicted doublets 
    (doublet_score and predicted_doublet kept in .obs).


    """
        
    sc.pp.scrublet(adata, random_state=random_state)
    return adata[~adata.obs["predicted_doublet"]].copy()



def normalisation_and_log_transformation(adata,target_sum=None):
    """
    Store raw counts, normalize total counts per cell and log-transform.

    Parameters
    ----------
    adata : AnnData
        QC-filtered + doublets filtered data with raw counts in adata.X.

    """
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata,target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.layers["log1p_norm"] = adata.X.copy()


def feature_selection(adata,layer="counts", flavor = "seurat_v3",n_top_genes = 2000):
    """
    select highly variable genes 

    Parameters
    ----------
    adata : AnnData

    n_top_genes : int
        Number of highly variable genes to keep.
    flavor : str
        Method passed to scanpy.pp.highly_variable_genes.

    Returns
    -------
    Modifies adata in place: adds adata.var['highly_variable']
    


    """
    sc.pp.highly_variable_genes(adata,
                                layer=layer,
                                flavor=flavor,
                                n_top_genes=n_top_genes )