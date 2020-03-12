# Automatic generation of h5ad files from spaceranger output

```.h5ad``` is a format for storing annotaded data, released in conjunction with the publication of Scanpy. It has gained a lot of traction and is highly suitable for storing Visium data; allowing one to store coordinates, count data, annotations, gene names, images and scaling factors in one single file. There's no established convention for how the Visum data should be stored, thus the structure described in **Structure** is used. Once installed, this package allows for easy conversion from the spaceranger output to the h5ad format, where the only required input is the spaceranger output directory.


## Structure 

* X : The raw count matrix [n\_spots x n\_genes].
* obs 
  * barcodes : 10x barcodes identifiers
  * under_tissue : Binary indiator if spot is under tissue or not (1 = is under,
    0 = is outside)
  * _x : array x-coordinates
  * _y : array y-coordinates
  * x : pixel x-coordinates (use these for visualization)
  * y : pixel y-coordinates (use these for visualization)
  * sample : string indicating which sample a spot belongs to
* var
  * name : HGNC gene symbols (duplicates may occur)
  * id : ENSEMBL annotations 
  * n_counts : total number of observed UMI's (transcripts) of a given gene
* uns
  * spot\_diameter\_fullres : diameter of spots for the full resolution image
  * tissue\_hires\_scalef : scaling factor, transforms full resolution pixel
    coordinates (obs.x and obs.y) to coordinates compatible with the hires image.
  * fiducial\_diameter\_fullres : diameter of fiducials for the full resolution image
  * image\_hires : HE-image 

## Install

A ```setup.py``` file is provided for easy installation simply run

```sh
./setup.py install
```

to install the package. Depending on your system/rights you might have to use the ```--user``` flag, exchanging the above command for:
```sh
./setup.py install --user
```




## Running

To generate a ```.h5ad``` file run

```sh
space2h5ad -dd SPACERANGER_OUTPUT_DIR

```

this will generate a file ```feature_matrix.h5ad``` in the folder ```SPACERANGER_OUTPUT_DIR```. 

You can specify another output file by adding the argument ```-o```, for example 
```sh
space2h5ad -dd SPACERANGER_OUTPUT_DIR -o /tmp/visium-data-sample-1.h5ad

```

by default the filtered data will be used, to change this (including all spots) add the flag ```--use_raw```, for example:

```sh
space2h5ad -dd SPACERANGER_OUTPUT_DIR --use_raw

```

by default ENSEMBL ids will be used, as gene names (index of vars), this can however be changed by adding the flag ```--gene_names```. For example:

```sh
space2h5ad -dd SPACERANGER_OUTPUT_DIR --gene_names

```

Since gene symbols may map to the same ENSEMBL id, there's not a one-to-one relationship between the two. To circumvent this the first instance of a gene name is kept whilst the others are discarded. This feature will likely be refined in the future.
