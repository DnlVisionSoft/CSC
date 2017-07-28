# CSC
The Cell Segmentation and Counting based on Gray Lvel Clustering (CSC) package.  

--------------------------------------------------------------------------

New technological advancements in automated microscopy gave rise to large volumes of data, which  have made human-based analysis infeasible and have heightened the need of automatic systems for high-throughput microscopy applications.
In particular in the field of fluorescence microscopy, automatic tools for image analysis make an essential contribution to increase the statistical power of the cell analysis process. 
This is a difficult task due to both diversification of staining patterns and local variability of the images. In order to cope with this challenges, this method divides the whole image in overlapping patches, which then undergo a gray level clustering. Thus, adaptive thresholding on clustered gray levels is applied to single patches to extract the foreground, while a merging process is implemented to recompose the foreground of overlapping patches into a final binary mask. The foreground represents the input of the cell counting stage. Since separating clustered cells is crucial for an accurate cell counting, centres of cells are detected by a two-stage process that combines the distance transform and curvature analysis. Labeling the detected centres, a  partition of the image in single cells is obtained. The method has been tested on several publicly available image datasets with respect to both segmentation and cell counting. 
This package offers a powerful tool for automatic cells segmentation and counting.

The foreground detection algorithm as well as the cell counting process is driven by a few number of parameters, whose purpose is detailed in the following. The image pre-processing is an optional step regulated by a boolean parameter, that is set to "true", when image correction is required and is "false", otherwise. 
The foreground detection is based on a sliding window, whose behavior is determined by two parameters that are the window size "n" and the sliding step "beta". The patch size "n" strongly depends on the resolution and homogeneity of the input image. Indeed, the more homogeneous are the objects into the image, the larger the value that can be assigned to "n". The parameter "beta" determines the degree of overlap of different patches. An appropriate value for "beta" can be selected according to the same heuristic adopted for "n".

In the binarization process, a key role is played by the parameter "alpha", as it regulates the adapthive thresholding of the quantized patches. In particular, we notice that higher values of "alpha" must be set when pre-processing is applied to the image. This is motivated by the fact that pre-processing produces significant changes in the image contrast by spreading the original grey levels on a larger range of values. On the contrary, when pre-processing is omitted, the quantized patches are generally characterized by low contrast that induces smaller values for the parameter "alpha".
The parameter "lambda" is also involved in the binarization process and regulates how much the thresholding of the current quantized patch is influenced by threshold values adopted for the preceding ones. The value of the parameter epsilon depends on the size of the smallest cell into the images. 
The only parameter involved in the counting process is "delta", which drives the incremental clustering of seeds produced by the distance transform. The value of this parameter strongly depends on the size of cells. Indeed, the smaller the cells, the lower the value of "delta" to be set, to avoid that different cells into a cluster will be merged.

--------------------------------------------------------------------------

Matlab Implementation by Daniel Riccio, June 06, 2017. 
Copyright (C) 2017 Daniel Riccio (dnl.riccio@gmail.com)
