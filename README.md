# HLC2
A highly efficient cross-matching framework for large astronomical catalogues in hybrid computing platforms

This is the main source code of HLC2, which can perform cross-matching between two astronomical catalogues. HLC2 is a high performance command-line cross-matching tool running on the Linux platform, which is implemented in C and C++.


- [HCGrid](#hcgrid)
  * [More About HCGrid](#more-about-hcgrid)
    + [Implementation](#implementation)
    + [Features](#features)
  * [Installation](#installation)
    + [Dependencies](#dependencies)
    + [Build from source](#build-from-source)
  * [Getting Started](#getting-started)
    + [Minimal example](#minimal-example)
  * [Community Contribution and Advice](#community-contribution-and-advice)



## More About HCGrid 


### Implementation

The specific steps of cross-matching for two astronomical catalogues are shown in the following figure. 

<P align="center"><img src=figs/architecture_and_workflow_of_HLC2.png hidth="50%" width="78%"></img></p>

First the **Extraction module** extracts celestial position information (mainly RA, DEC) from the input data and filters out useless information. The **First-level partition module** divides the extracted location information into HEALPix data blocks, while our quad-direction strategy is implemented in the **Boundary-solved module** to reduce the loss of accuracy without adding too much redundant data. Then through **Second-level partition module**, calculation blocks can be obtained for storing and subsequent parallel accessing catalogue records on GPU. **Source reading module** is designed to retrieve the current CPU and GPU computing status to dynamically adjust splitting strategy. On the GPU, the inter-catalogue parallelization is adopted to calculate the radius distances on **Kernel module** with I/O optimizations on the **Compression module**. Finally, the matching results transferred to CPU will be exported the final products and be visualized. This function is under development. 


### Features

- Supports WCS projection system as target.
- Scales well in CPU-GPU heterogeneous platforms.

## Installation

### Dependencies


- CUDA Toolkit

All of these packages can be found in "Dependencies" directory or get from follow address:


- CUDA: https://developer.nvidia.com/cuda-toolkit-archive

### Build from source

1. Change the direction to "HCGrid" folder
2. Update the library file paths in the Makefile according to the paths of installed dependencies, e.g. CUDA, fits, wcslib, etc.
3. make 

## Getting Started

1. Defined and create a target grid map according to specific scenario, in Creat_target_file.py, such as:

``` python
# define target FITS/WCS header
header = {
	'NAXIS': 3,
	'NAXIS1': dnaxis1,
	'NAXIS2': dnaxis2,
	'NAXIS3': 1,  
	'CTYPE1': 'RA---SIN',
	'CTYPE2': 'DEC--SIN',
	'CUNIT1': 'deg',
	'CUNIT2': 'deg',
	'CDELT1': -pixsize,
	'CDELT2': pixsize,
	'CRPIX1': dnaxis1 / 2.,
	'CRPIX2': dnaxis2 / 2.,
	'CRVAL1': mapcenter[0],
	'CRVAL2': mapcenter[1],
	}
```

2. Set the related kernel parameters in HCGrid.cpp

```c++
/*Set kernel*/
kernel_type = GAUSS1D;
kernelsize_fwhm = 300./3600.;
kernelsize_sigma = 0.2;
kernel_params[0] = kernelsize_sigma;
sphere_radius = 3.*kernelsize_sigma;
hpx_max_resolution=kernelsize_sigma/2;
_prepare_grid_kernel(
	kernel_type, 
	kernel_params, 
	sphere_radius, 
	hpx_max_resolution
	);
```

3. make

### Minimal example

In the terminal window, after successful compilation, you can do the following thing:

1. Type "**./HCGrid -h**" to get the detail parameter guide.
2. ./HCGrid [options]. The options include the following parameter:

| Parameter | Description |
| :----------| :-----------------------------------|
| fits_path  | Absolute path of FITS file             |
| input_file | Name of unsorted input FITS file       |
| target_file| Name of target FITS file               |
| output_file| Name of output FITS file               |
| sorted_file| Name of sorted input FITS file         |
| fits_id    | ID of FITS file                        |
| beam_size  | Beam size of FITS file                 |
|register_num| total number of registers for each thread block of the GPU |
|  sp_num    |the number of SPs in each SM of the GPU.|
| ord_arg    | Select the pre_order function      |
| block_num  | The number of thread in each block     |
| coarsening_factor| The value of coarsening factor   |

3. Create the target map:

```shell
$ python Creat_target_file.py -p /home/summit/Project/HCGrid/data/ -t target -n 1 -b 300
```

**Note:** You need to set the relevant parameters of target_map according to the coverage sky area and beam width of the sampled data. For details, please refer to "Creat_target_file.py " file.

4. Do the gridding:

```shell
$ ./HCGrid --fits_path /home/summit/Project/HCGrid/data/ --input_file input --target_file target --output_file output --fits_id 1 --beam_size 300 --register_num 64 --sp_num 64 --order_arg 1
```

or

```shell
$ ./HCGrid --fits_path /home/summit/HCGrid/data/ --input_file input --target_file target --output_file output --fits_id 1 --beam_size 300 --order_arg 1 --block_num 64
```

The former further specifies the relevant hardware parameters, please refer to our article for details.



## Community Contribution and Advice

HLC2 is being further improved, if you have any question or ideas, please don't skimp on your suggestions and welcome make a pull request. Moreover, you can contact us through the follow address.

- zyj0928@tju.edu.cn

