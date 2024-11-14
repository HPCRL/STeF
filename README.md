# STeF: Sparse Tensor Decomposition Framework

STeF (Sparse Tensor Decomposition Framework) is a C++ library that provides efficient implementations of sparse tensor decomposition algorithms, with a focus on Canonical Polyadic Decomposition (CPD). This framework is designed to accelerate data analytics and machine learning tasks that involve sparse tensors.

## Overview

Sparse tensor decomposition, such as Canonical Polyadic Decomposition (CPD), plays a crucial role in data analytics and machine learning. The computational complexity of CPD is dominated by Matricized Tensor Times Khatri Rao Product (MTTKRP) operations. Sparse CPD requires performing a collection of MTTKRP operations, with common sub-computations across them. Several existing approaches aim to factorize and reuse these common sub-expressions to minimize computational overhead. However, previous work primarily focused on reducing the number of high-level operators.

In this code, we present an innovative design space exploration for sparse CPD. Our exploration covers various factors, including whether partial MTTKRP results should be saved, different mode permutations, and modeling the total volume of data movement to and from memory. Additionally, we propose a fine-grained load balancing method that supports higher levels of parallelization.

This repository houses the C++ implementation of our proposed techniques and optimizations discussed in the paper. By leveraging this code, you can explore different strategies to enhance the performance and efficiency of sparse CPD. Our implementations include features such as saving partial MTTKRP results, mode permutations, and minimizing data movement to and from memory. Furthermore, we introduce a fine-grained load-balancing method that enables efficient parallelization.

## Paper

For more detailed information, please refer to our research paper, titled "Sparse Tensor Decomposition: Exploring Design Space and Load Balancing for Canonical Polyadic Decomposition". The full paper can be found [here](https://ieeexplore.ieee.org/abstract/document/9820702).

## Getting Started

To utilize the STeF framework, follow the instructions below:

1. Clone this repository:

```git clone https://github.com/your-username/STeF.git```


2. Build the C++ code using make.

``` make ```

## How to Use 

Currently, a clean-up is in progress. You can use the version described as `STeF` in "Sparsity-Aware Tensor Decomposition" paper as follows. Changing the cache size parameter will only affect the model data size movement and it should be in Bytes. 
```
./bin/STeF.exe <tensor name> <number of ranks (optional, defaults to 32)> <cache size (optional, defaults to 1MB)>
```
