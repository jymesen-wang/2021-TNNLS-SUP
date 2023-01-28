# SUP

To use the function "SUP" or "GOSUP", please follow the input/output format:

 [P,F] = SUP(X,S,d1,k,lambda)
 [P,F] = GOSUP(X,S,d1,k,lambda)

X   ： d*n samples matrix, d and n are dimension and number of samples, respectively
S   :  n*n the similarity matrix between samples
d1  :  subspace dimension, d1 < d
k   ： int， d1 < k < d
lambda: regularization parameter
P   ： d*d1 row-sparse projection matrix
F   ： n*d1 samples after noise removal in subspace

Please make sure that the documents Eu2_distance.m and ClusteringMeasure.m are in the same folder as SUP.m and GOSUP.m

Use the codes, please cite
Wang J, Wang L, Nie F, et al. Joint Feature Selection and Extraction With Sparse Unsupervised Projection[J]. IEEE Transactions on Neural Networks and Learning Systems, 2021, doi: 10.1109/TNNLS.2021.3111714.

If you have any questions, please connect wanglinjun@mail.nwpu.edu.cn
