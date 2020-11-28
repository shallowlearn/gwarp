# gwarp++

Software to create true ortho photos of satellite images using rational polynomial coefficients (RPCs) and Digital Surface Models (DSMs).

While the popular and powerful [gdalwarp](https://gdal.org/programs/gdalwarp.html) tool works well with a Digital Elevation Model (DEM), it does not account for the heights of elevated structures such as buildings and trees in a DSM.

On the other hand, gwarp++ detects the occluded portions of the ground and creates a true ortho image. The occluded portions are marked with a no-data value. 

If you use gwarp++ in your research, **we would appreciate your citing our [paper](https://arxiv.org/abs/2008.10271v2).**

```bash
@misc{comandur2020semantic,
      title={Semantic Labeling of Large-Area Geographic Regions Using Multi-View and Multi-Date Satellite Images and Noisy OSM Training Labels},
      author={Bharath Comandur and Avinash C. Kak},
      year={2020},
      eprint={2008.10271},
      archivePrefix={arXiv},
      primaryClass={cs.CV}
}
```

For additional details regarding gwarp++, please see Sections 2.5 and 5.1.2 of the dissertation [here](https://hammer.figshare.com/articles/thesis/Semantic_Labeling_of_Large_Geographic_Areas_Using_Multi-Date_and_Multi-View_Satellite_Images_and_Noisy_OpenStreetMap_Labels/12739556).

## Requirements:

Install gdal (with C++ and python bindings). It is recommended to use an anaconda environment for this.

Also install libboost program options as follows:
```bash
$ sudo apt-get install -y libboost-program-options-dev
```
### Compile gwarp C++
```bash
env=<Name of anaconda environment>

$ conda activate $env

$ bash compile.sh
```
### Run orthorectification

There are two ways to do this. We can either use "orthorectify.py" to orthorectify multiple small images in parallel. Otherwise we can use "gwarp++.py" to orthorectify a single large image by dividing it into smaller pieces and orthorectifying the pieces in parallel.

#### Orthorectify multiple images in parallel
```python
$ python orthorectify.py [-h] -folder_img FOLDER_IMG
  	 		      [-folder_rpc FOLDER_RPC]
                       	      -dsm_file DSM_FILE
			      -dem_file DEM_FILE
			      -folder_output FOLDER_OUTPUT
			      -cache_dir CACHE_DIR
			      [-parallel]
                       	      [-num_workers NUM_WORKERS]
                       	      [-dsm_mask_file DSM_MASK_FILE]
                       	      [-image_interp IMAGE_INTERP]
			      [-output_res OUTPUT_RES]
                       	      [-nodata NODATA]
```

##### Required Arguments

``` -folder_img FOLDER_IMG ``` - Folder containing images to orthorectify. Images should be in GTiff file format. The tifs can have the RPCs embedded in their metadata.

``` -dsm_file DSM_FILE ``` - DSM file with heights respect to WGS84 ellipsoid

``` -dem_file DEM_FILE ``` - DEM file indicating height of ground with respect to WGS84 ellipsoid

``` -folder_output FOLDER_OUTPUT ``` - Output folder

``` -cache_dir CACHE_DIR ``` - A cache/temporary directory. It should not exist. The program creates this directory and deletes it at the end. Please make sure it is not in any important folder such as /root or / .

##### Optional Arguments

``` -folder_rpc FOLDER_RPC ``` - Folder containing RPCs in RPB file format. There should be a one-to-one correspondence between the name of the tif files and the RPB files, i.e., an image A.tif should have a RPB file A.RPB. If this folder is specified, then only the images that have corresponding RPB files will be processed. Moreover this RPB file can be used to override any RPCs that are embedded in the GTiff metadata.

``` -parallel ``` - To process large tifs in parallel. The software chops up the image into smaller pieces. Warning - Consumes more memory.

``` -num_workers NUM_WORKERS ``` - Number of parallel processes to launch. The more workers, the more memory is used.

``` -dsm_mask_file DSM_MASK_FILE ``` - An optional file indicating points to ignore in the DSM. For example, water bodies might have noisy height values in a DSM.

``` -image_interp IMAGE_INTERP ``` - Interpolation method for the image. Can be either 'near' or 'bilinear'.

``` -output_res OUTPUT_RES ``` - Output ground sampling distance (GSD) in meters. Default is 0.5 m per pixel. Please note that the final output GTiff has its projection in the WGS84 coordinate system. The sampling distances in longitude and latitude are estimated using the specified output_res.

``` -nodata NODATA ``` - A value to indicate the occluded points in the final true ortho image. Default is -9999.

#### Quickly orthorectify a large image
```bash
Coming Soon
```

## Acknowledgments

```bash
This work was supported by the Intelligence Advanced
Research Projects Activity (IARPA) via Department of
Interior / Interior Business Center (DOI/IBC) contract
number D17PC00280. The U.S. Government is authorized to
reproduce and distribute reprints for Governmental purposes
notwithstanding any copyright annotation
thereon. Disclaimer: The views and conclusions contained
herein are those of the authors and should not be
interpreted as necessarily representing the official
policies or endorsements, either expressed or implied, of
IARPA, DOI/IBC, or the U.S. Government.
```
