# HLC2
A highly efficient cross-matching framework for large astronomical catalogues in hybrid computing platforms

This is the main source code of HLC2, which can perform cross-matching between two astronomical catalogues. HLC2 is a high performance command-line cross-matching tool running on the Linux platform, which is implemented in C and C++.


- [HLC2](#hcgrid)
  * [More About HLC2](#more-about-hcgrid)
    + [Implementation](#implementation)
    + [Features](#features)
  * [Installation](#installation)
    + [Dependencies](#dependencies)
    + [Build from source](#build-from-source)
  * [Getting Started](#getting-started)
    + [Minimal example](#minimal-example)
  * [Community Contribution and Advice](#community-contribution-and-advice)



## More About HLC2


### Implementation

The specific steps of cross-matching for two astronomical catalogues are shown in the following figure. 

<P align="center"><img src=figs/architecture_and_workflow_of_HLC2.png hidth="50%" width="75%"></img></p>

First the **Extraction module** extracts celestial position information (mainly RA, DEC) from the input data and filters out useless information. The **First-level partition module** divides the extracted location information into HEALPix data blocks, while our quad-direction strategy is implemented in the **Boundary-solved module** to reduce the loss of accuracy without adding too much redundant data. Then through **Second-level partition module**, calculation blocks can be obtained for storing and subsequent parallel accessing catalogue records on GPU. **Source reading module** is designed to retrieve the current CPU and GPU computing status to dynamically adjust splitting strategy. On the GPU, the inter-catalogue parallelization is adopted to calculate the radius distances on **Kernel module** with I/O optimizations on the **Compression module**. Finally, the matching results transferred to CPU will be exported the final products and be visualized. This function is under development. 


### Features

- Supports end-to-end efficent cross-matching of catalogues.
- Scales well in CPU-GPU heterogeneous platforms.

## Installation

### Dependencies


- CUDA Toolkit

All of these packages can be found in "Dependencies" directory or get from follow address:


- CUDA: https://developer.nvidia.com/cuda-toolkit-archive

## Build from source

1. Change the direction to 'HLC2' directory
2. Update the library file paths in the Makefile according to the paths of installed dependencies.
3. make
It will generate an executable file：Crossmatch

### Minimal example

In the terminal window, after successful compilation, you can do the following thing:

1. Type "**./HCGrid -h**" to get the detail parameter guide.

3. Create the target map:

```shell
$ python Creat_target_file.py -p /home/summit/Project/HCGrid/data/ -t target -n 1 -b 300
```

**Note:** You need to set the relevant parameters of target_map according to the coverage sky area and beam width of the sampled data. For details, please refer to "Creat_target_file.py " file.

4. Do the gridding:

```shell
$ ./HCGrid --fits_path /home/summit/Project/HCGrid/data/ --input_file input --target_file target --output_file output --fits_id 1 --beam_size 300 --register_num 64 --sp_num 64 --order_arg 1
```

## Community Contribution and Advice

HLC2 is being further improved, if you have any question or ideas, please don't skimp on your suggestions and welcome make a pull request. Moreover, you can contact us through the follow address.

- zyj0928@tju.edu.cn

