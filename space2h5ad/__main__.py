#!/usr/bin/env python3

import argparse as arp
import os.path as osp

import anndata as ad
import h5py

import numpy as np
import pandas as pd

import scipy.sparse as sp_sparse

from PIL import Image

import json
 

def format_spaceranger_output(data_dir : str,
                              use_hgnc : bool = False,
                              use_raw : bool = False,
                              ) -> ad.AnnData:

    if use_raw:
        count_type = 'raw'
    else:
        count_type = 'filtered'

    pths = dict(data = osp.join(data_dir,
                                count_type + "_feature_bc_matrix.h5"),
                spots = osp.join(data_dir,
                                 "spatial",
                                 "tissue_positions_list.csv"),
                img = osp.join(data_dir,
                               "spatial",
                               "tissue_hires_image.png"),
                scf = osp.join(data_dir,
                               "spatial",
                               "scalefactors_json.json"),
                )

    file_status = [osp.exists(p) for p in pths.values()]

    if not all(file_status):
        print("[ERROR] : Data is incomplete.")
        for s,f in zip(file_status,pths.values):
            if not s:
                print("\tMissing : {}".format(f))
        print("Exiting")
        sys.exit(-1)

    else:
        print("[INFO] : Loading Data")

    with h5py.File(pths["data"],'r') as f:
        raw_data = f['matrix']
        matrix = sp_sparse.csc.csc_matrix((np.asarray(raw_data.get('data')),
                                        np.asarray(raw_data.get("indices")),
                                        np.asarray(raw_data.get("indptr")),
                                        ))

        barcodes =  pd.Index(np.asarray(raw_data.get('barcodes')).astype(np.str))
        name_id = np.asarray(raw_data['features']['id'] ).astype(np.str)
        name_hgnc = np.asarray( raw_data['features']['name']  ).astype(np.str)


    with open(pths['scf'],'r+') as jopen:
        scf = json.load(jopen)

    img_high = np.asarray(Image.open(pths['img']))

    scf.update({"image_hires":img_high})
    _ = scf.pop("tissue_lowres_scalef")

    spt = pd.read_csv(pths['spots'],
                    sep = ',',
                    header = None,
                    index_col = None)

    print("[INFO] : Formatting data")

    spt.columns = ['barcodes',
                'under_tissue',
                '_x',
                '_y',
                'x',
                'y']

    spt.index = spt['barcodes']


    cnt = pd.DataFrame(matrix.todense().T,
                    index = barcodes,
                    columns = name_id,
                    )

    barcodes = spt.index.intersection(cnt.index)

    spt = spt.loc[barcodes,:]
    cnt = cnt.loc[barcodes,:]

    var = dict(name = name_hgnc,
               id = name_id,
               n_counts = cnt.sum(axis = 0),
               )



    var = pd.DataFrame(var)


    if use_hgnc:

        _,idx = np.unique(var['name'].values,
                            return_index = True)
        idx = np.sort(idx)
        genes = var['name'].values[idx]

        var = var.iloc[:,idx]
        cnt = cnt.iloc[:,idx]

    var.index = cnt.columns

    adata = ad.AnnData(cnt.values,
                    obs = spt,
                    var = var,
                    uns = scf,
                    )

    return adata

def main():

    prs = arp.ArgumentParser()

    aa = prs.add_argument

    aa('-dd',"--data_dir",type = str,required = True)
    aa('-o',"--output",type = str, default = None)
    aa('-gn',"--gene_names",type = bool, default = False)
    aa('-ur',"--use_raw",type = bool, default = False,action = 'store_true')


    args = prs.parse_args()

    adata = format_spaceranger_output(args.data_dir,
                                      args.gene_names,
                                      args.use_raw,
                                      )

    if args.output is None:
        out_pth = osp.join(args.data_dir,"feature_matrix.h5ad")
    elif not args.output.split('.')[-1] == 'h5ad':
        out_pth = args.output + '.h5ad'
    else:
        out_pth = args.output

    print("[INFO] : Saving Data to : {}".format(out_pth))

    adata.write(out_pth)

if __name__ == '__main__':
    main()

