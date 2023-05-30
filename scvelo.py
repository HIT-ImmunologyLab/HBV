#!/usr/bin
# -*-coding:utf-8-*-


import os
import anndata
import cellrank as cr
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
from scipy import io


def mkdir(path):
    cmd = f"mkdir -p {path}"
    if not os.path.exists(path):
        os.system(cmd)


def load_adata(file_dir):
    dir_ = file_dir

    sparse_path = os.path.join(dir_, 'counts.mtx')
    X = io.mmread(sparse_path)
    adata = anndata.AnnData(X=X.transpose().tocsr())

    meta_path = os.path.join(dir_, 'metadata.csv')
    cell_meta = pd.read_csv(meta_path)

    gene_names_path = os.path.join(dir_, 'gene_names.csv')
    with open(gene_names_path, 'r') as f:
        gene_names = f.read().splitlines()

    # set anndata observations and index obs by barcodes, var by gene names
    adata.obs = cell_meta
    adata.obs.index = adata.obs['barcode']
    adata.var.index = gene_names


    # set pca and umap
    adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

    return (adata)


def plot_scvelo(file_dir, grade):
    adata = load_adata(file_dir)
    if grade == "IA":
        ldataA01 = scv.read('A01.loom', cache=True)
        ldataA02 = scv.read('A02.loom', cache=True)
        ldataA04 = scv.read('A04.loom', cache=False)
        ldataA06 = scv.read('A06.loom', cache=False)
        ldataA07 = scv.read('A07.loom', cache=False)


        ##############
        barcodes = [bc.split(':')[1] for bc in ldataA01.obs.index.tolist()]
        ldataA01.obs.index = barcodes
        barcodes = [bc.split(':')[1] for bc in ldataA02.obs.index.tolist()]
        ldataA02.obs.index = barcodes
        barcodes = [bc.split(':')[1] for bc in ldataA04.obs.index.tolist()]
        ldataA04.obs.index = barcodes
        barcodes = [bc.split(':')[1] for bc in ldataA07.obs.index.tolist()]
        ldataA07.obs.index = barcodes
        barcodes = [bc.split(':')[1] for bc in ldataA06.obs.index.tolist()]
        ldataA06.obs.index = barcodes


        #############

        ldataA01.var_names_make_unique()
        adata_sub = adata[adata.obs['orig.ident'] == "A01_LIVER_JGB", :]
        adataA01 = scv.utils.merge(adata_sub, ldataA01)

        ldataA02.var_names_make_unique()
        adata_sub = adata[adata.obs['orig.ident'] == "A02_LIVER_ZLY", :]
        adataA02 = scv.utils.merge(adata_sub, ldataA02)

        ldataA04.var_names_make_unique()
        adata_sub = adata[adata.obs['orig.ident'] == "A04_LIVER_GY", :]
        adataA04 = scv.utils.merge(adata_sub, ldataA04)

        ldataA07.var_names_make_unique()
        adata_sub = adata[adata.obs['orig.ident'] == "A07_XHL_LIVER", :]
        adataA07 = scv.utils.merge(adata_sub, ldataA07)

        ldataA06.var_names_make_unique()
        adata_sub = adata[adata.obs['orig.ident'] == "A06_SXH_LIVER", :]
        adataA06 = scv.utils.merge(adata_sub, ldataA06)

       

        adata = adataA01.concatenate(
            [adataA02, adataA04, adataA06, adataA07])
        sc.pl.umap(adata, color=['donor_name'], frameon=False, save=True)



    elif grade == "IT":
        ldataT01 = scv.read('T01.loom', cache=True)
        ldataT02 = scv.read('T02.loom', cache=False)
        ldataT04 = scv.read('T04.loom', cache=False)
        ldataT05 = scv.read('T05.loom', cache=False)
        ldataT06 = scv.read('T06.loom', cache=False)

        barcodes = [bc.split(':')[1] for bc in ldataT01.obs.index.tolist()]
        barcodes = [bc[0:len(bc) - 1] + '_2' for bc in barcodes]
        ldataT01.obs.index = barcodes
        barcodes = [bc.split(':')[1] for bc in ldataT02.obs.index.tolist()]
        barcodes = [bc[0:len(bc) - 1] + '_7' for bc in barcodes]
        ldataT02.obs.index = barcodes

        barcodes = [bc.split(':')[1] for bc in ldataT04.obs.index.tolist()]
        barcodes = [bc[0:len(bc) - 1] + '_8' for bc in barcodes]
        ldataT04.obs.index = barcodes

        barcodes = [bc.split(':')[1] for bc in ldataT05.obs.index.tolist()]
        barcodes = [bc[0:len(bc) - 1] + '_10' for bc in barcodes]
        ldataT05.obs.index = barcodes

        barcodes = [bc.split(':')[1] for bc in ldataT06.obs.index.tolist()]
        barcodes = [bc[0:len(bc) - 1] + '_10' for bc in barcodes]
        ldataT06.obs.index = barcodes
        ldataT01.var_names_make_unique()
        ldataT02.var_names_make_unique()
        ldataT04.var_names_make_unique()
        ldataT05.var_names_make_unique()
        ldataT06.var_names_make_unique()

        adata_sub = adata[adata.obs['orig.ident'] == "T01_LIVER_WJP", :]
        adataT01 = scv.utils.merge(adata_sub, ldataT01)
        
        adata_sub = adata[adata.obs['orig.ident'] == "T02_LIVER_HYJ", :]
        adataT02 = scv.utils.merge(adata_sub, ldataT02)

        adata_sub = adata[adata.obs['orig.ident'] == "T04_LIVER_LHM", :]
        adataT04 = scv.utils.merge(adata_sub, ldataT04)

        adata_sub = adata[adata.obs['orig.ident'] == "T06_CC_LIVER", :]
        adataT06 = scv.utils.merge(adata_sub, ldataT06)

        adata_sub = adata[adata.obs['orig.ident'] == "T05_LIVER_LXY", :]
        adataT05 = scv.utils.merge(adata_sub, ldataT05)

       

        adata = adataT01.concatenate(
            [adataT02, adataT04, adataT05, adataT06])
        sc.pl.umap(adata, color=['donor_name'], frameon=False, save=True)
    elif grade == "IC":
        ldataC01 = scv.read('C01.loom', cache=False)
        ldataC03 = scv.read('C03.loom', cache=False)
        ldataC04 = scv.read('C04.loom', cache=False)
        ldataC05 = scv.read('C05.loom', cache=False)
        ldataC06 = scv.read('C06.loom', cache=False)

        barcodes = [bc.split(':')[1] for bc in ldataC01.obs.index.tolist()]
        barcodes = [bc[0:len(bc) - 1] + '_4' for bc in barcodes]
        ldataC01.obs.index = barcodes

        barcodes = [bc.split(':')[1] for bc in ldataC03.obs.index.tolist()]
        barcodes = [bc[0:len(bc) - 1] + '_9' for bc in barcodes]
        ldataC03.obs.index = barcodes

        barcodes = [bc.split(':')[1] for bc in ldataC04.obs.index.tolist()]
        barcodes = [bc[0:len(bc) - 1] + '_11' for bc in barcodes]
        ldataC04.obs.index = barcodes

        barcodes = [bc.split(':')[1] for bc in ldataC05.obs.index.tolist()]
        barcodes = [bc[0:len(bc) - 1] + '_12' for bc in barcodes]
        ldataC05.obs.index = barcodes

        barcodes = [bc.split(':')[1] for bc in ldataC06.obs.index.tolist()]
        barcodes = [bc[0:len(bc) - 1] + '_13' for bc in barcodes]
        ldataC06.obs.index = barcodes

        ldataC01.var_names_make_unique()
        ldataC03.var_names_make_unique()
        ldataC04.var_names_make_unique()
        ldataC05.var_names_make_unique()
        ldataC06.var_names_make_unique()

        adata_sub = adata[adata.obs['orig.ident'] == "C01_LIVER_WXM", :]
        adataC01 = scv.utils.merge(adata_sub, ldataC01)

        adata_sub = adata[adata.obs['orig.ident'] == "C04_LIVER_WYM", :]
        adataC04 = scv.utils.merge(adata_sub, ldataC04)

        adata_sub = adata[adata.obs['orig.ident'] == "C03_LIVER_GJT", :]
        adataC03 = scv.utils.merge(adata_sub, ldataC03)

        adata_sub = adata[adata.obs['orig.ident'] == "C05_LIVER_MFB", :]
        adataC05 = scv.utils.merge(adata_sub, ldataC05)

        adata_sub = adata[adata.obs['orig.ident'] == "C06_LIVER_YSC", :]
        adataC06 = scv.utils.merge(adata_sub, ldataC06)

        

        adata = adataC01.concatenate(
            [adataC03, adataC04, adataC05, adataC06])
        sc.pl.umap(adata, color=['donor_name'], frameon=False, save=True)
    elif grade == "AHB":

        ldataAHB01 = scv.read('AHB01.loom', cache=False)
        ldataAHB02 = scv.read('AHB02.loom', cache=False)
        ldataAHB03 = scv.read('AHB03.loom', cache=False)

        barcodes = [bc.split(':')[1] for bc in ldataAHB01.obs.index.tolist()]
        barcodes = [bc[0:len(bc) - 1] + '_1' for bc in barcodes]
        ldataAHB01.obs.index = barcodes

        barcodes = [bc.split(':')[1] for bc in ldataAHB02.obs.index.tolist()]
        barcodes = [bc[0:len(bc) - 1] + '_3' for bc in barcodes]
        ldataAHB02.obs.index = barcodes

        barcodes = [bc.split(':')[1] for bc in ldataAHB03.obs.index.tolist()]
        barcodes = [bc[0:len(bc) - 1] + '_3' for bc in barcodes]
        ldataAHB03.obs.index = barcodes

        ldataAHB01.var_names_make_unique()
        ldataAHB02.var_names_make_unique()
        ldataAHB03.var_names_make_unique()

        adata_sub = adata[adata.obs['orig.ident'] == "AHB01_CY_LIVER", :]
        adataAHB01 = scv.utils.merge(adata_sub, ldataAHB01)

        adata_sub = adata[adata.obs['orig.ident'] == "AHB02_TKL_LIVER", :]
        adataAHB02 = scv.utils.merge(adata_sub, ldataAHB02)

        adata_sub = adata[adata.obs['orig.ident'] == "AHB03_SLP_LIVER", :]
        adataAHB03 = scv.utils.merge(adata_sub, ldataAHB03)

        adata = adataAHB01.concatenate([adataAHB02, adataAHB03])

        sc.pl.umap(adata, color=['donor_name'], frameon=False, save=True)
    elif grade == "HC":
        ldataCon01 = scv.read('Con01.loom', cache=True)
        ldataCon02 = scv.read('Con02.loom', cache=True)
        ldataCon03 = scv.read('Con03.loom', cache=True)

        barcodes = [bc.split(':')[1] for bc in ldataCon01.obs.index.tolist()]
        barcodes = [bc[0:len(bc) - 1] + '_1' for bc in barcodes]
        ldataCon01.obs.index = barcodes

        barcodes = [bc.split(':')[1] for bc in ldataCon02.obs.index.tolist()]
        barcodes = [bc[0:len(bc) - 1] + '_3' for bc in barcodes]
        ldataCon02.obs.index = barcodes

        barcodes = [bc.split(':')[1] for bc in ldataCon03.obs.index.tolist()]
        barcodes = [bc[0:len(bc) - 1] + '_5' for bc in barcodes]
        ldataCon03.obs.index = barcodes

        ldataCon01.var_names_make_unique()
        ldataCon02.var_names_make_unique()
        ldataCon03.var_names_make_unique()

        adata_sub = adata[adata.obs['orig.ident'] == "HC01_LIVER", :]
        adataCon01 = scv.utils.merge(adata_sub, ldataCon01)

        adata_sub = adata[adata.obs['orig.ident'] == "HC02_LIVER", :]
        adataCon02 = scv.utils.merge(adata_sub, ldataCon02)

        adata_sub = adata[adata.obs['orig.ident'] == "HC03_LIVER", :]
        adataCon03 = scv.utils.merge(adata_sub, ldataCon03)

        adata = adataCon01.concatenate([adataCon02, adataCon03])
        sc.pl.umap(adata, color=['donor_name'], frameon=False, save=True)
        #sc.pl.tsne(adata, color=['donor_name'], frameon=False, save=True)
    sc.pl.umap(adata, color='new_cluster', frameon=False, legend_loc='on data', title='', save='celltype_plot.pdf',
               legend_fontsize=5)
    scv.settings.verbosity = 3
    scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
    cr.settings.verbosity = 2
    # pre-process
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata,n_pcs=20, n_neighbors=20)
    # compute velocity
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode='dynamical')
    scv.tl.velocity_graph(adata)

    scv.pl.velocity_embedding_grid(adata, basis='umap', color='new_cluster', save='embedding_grid.pdf', title='',
                                   scale=0.25,size=23)
    scv.pl.velocity_embedding_stream(adata, basis='umap', color='new_cluster', save='embedding_stream.svg',
                                     legend_fontsize=5,size=23,
                                     title='')


if __name__ == '__main__':
    grades = ['AHB','IA', 'IT', 'IC','HC']
    file_dir = './HBV/velocyto'

    for grade in grades:
        out_dir = f"./HBV/velocyto/result/{grade}"
        mkdir(out_dir)
        os.chdir(out_dir)
        plot_scvelo(file_dir, grade)


