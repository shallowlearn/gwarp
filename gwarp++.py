'''
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

Author: Bharath Comandur, cjrbharath@gmail.com
Date: 11/24/2020
'''

from __future__ import print_function
import argparse 
import os, sys, glob
import numpy as np 
import gdal, osr
import pyproj
import pprint as pp 

import subprocess
import shutil

from scipy import interpolate

import time

from multiprocessing import Pool

import copy

import traceback

#from skimage import morphology

from rpc_parser import load_rpb
    
import pickle,cv2

codePath = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

# GDAL data type to the corresponding Numpy data type.
DICT_GDAL_TO_NP = {gdal.GDT_Byte : np.dtype('uint8'),
                   gdal.GDT_UInt16 : np.dtype('uint16'),
                   gdal.GDT_UInt32 : np.dtype('uint32'),
                   gdal.GDT_Int16 : np.dtype('int16'),
                   gdal.GDT_Int32 : np.dtype('int32'),
                   gdal.GDT_Float32 : np.dtype('float32'),
                   gdal.GDT_Float64 : np.dtype('float64')
}
###############################################################################

def extract_subraster_mem_efficient(rasterds, rasterband, bbox, fill_value = 0, raiseError = True):
    '''
    extract subraster with empty padding.
    bottom right corner of bbox is included
    This is primarily useful for gdal when the image is too large to be loaded fully
    '''
    
    raster_h = rasterds.RasterYSize
    raster_w = rasterds.RasterXSize
    rasterdtype = rasterds.GetRasterBand(1).DataType
    rasterdtype = DICT_GDAL_TO_NP[rasterdtype]
    
    top_left_col_row = bbox[0]
    bottom_right_col_row = bbox[1]
    
    if bbox[0][0] > raster_w - 1 or bbox[0][1] > raster_h - 1 or bbox[1][0] < 0 or bbox[1][1] < 0:
        if raiseError:
            raise ValueError("bbox outside raster")
        else:
            print ("bbox outside raster, returning empty array")
            height = bottom_right_col_row[1] - top_left_col_row[1] + 1
            width = bottom_right_col_row[0] - top_left_col_row[0] + 1
            empty_array = fill_value*np.ones([height,width])
            return empty_array.astype(rasterdtype)
    
    ### Extremely important to use int/np.int and not np.int64. Refer this bug - https://github.com/conda-forge/gdal-feedstock/issues/167
    ### gdal ReadAsArray will throw strange error if it is np.int64 claiming that it is double
    ### By default np.minimum and np.maximum will return np.int64 even if the input to it are ints
    left_zero_pad = int(np.abs(np.minimum(0, top_left_col_row[0])))
    
    top_zero_pad = int(np.abs(np.minimum(0, top_left_col_row[1])))
    
    right_zero_pad = int(np.abs(np.minimum(0, raster_w - bottom_right_col_row[0] - 1)))
    
    bottom_zero_pad = int(np.abs(np.minimum(0, raster_h - bottom_right_col_row[1] - 1 )))
    
    top_left_col = np.int(np.maximum(0, top_left_col_row[0]))
    
    top_left_row =  int(np.maximum(0, top_left_col_row[1]))
    
    bottom_right_col = int(np.minimum(raster_w, bottom_right_col_row[0] + 1))
    
    bottom_right_row = int(np.minimum(raster_h, bottom_right_col_row[1] + 1))
    
    return cv2.copyMakeBorder(rasterband.ReadAsArray(top_left_col, top_left_row, bottom_right_col - top_left_col, bottom_right_row - top_left_row), 
                              top_zero_pad, bottom_zero_pad, left_zero_pad, right_zero_pad, cv2.BORDER_CONSTANT, value = fill_value)


def verify_if_mask_file_is_needed(start_lon, start_lat, end_lon, end_lat, spacing_lon, spacing_lat, mask_file):
    
    verify_flag = False
    mask_file_ds = gdal.Open(mask_file)
    mask_file_gt = mask_file_ds.GetGeoTransform()
    padding = 200        
    bbox = [[ start_lon - (padding*spacing_lon), start_lat - (padding*spacing_lat) ],
            [ end_lon + (padding*spacing_lon), end_lat + (padding*spacing_lat) ] ]
    bbox_pixels = convert_projected_coords_to_pixel_coords(np.array(bbox), mask_file_gt)  
    bbox_pixels = np.floor(bbox_pixels).astype(np.int)
    
    """
    bbox_pixels[0,0] = np.maximum(0, bbox_pixels[0,0])
    bbox_pixels[0,1] = np.maximum(0, bbox_pixels[0,1])
    
    bbox_pixels[1,0] = np.minimum(mask_file_ds.RasterXSize - 1, bbox_pixels[1,0])
    bbox_pixels[1,1] = np.minimum(mask_file_ds.RasterYSize - 1, bbox_pixels[1,1])    
    
    bbox_pixels = np.floor(bbox_pixels).astype(np.int)
    """
    
    bbox_vals = extract_subraster_mem_efficient(mask_file_ds, mask_file_ds.GetRasterBand(1), bbox_pixels, fill_value = 0, raiseError = False)
                
    if np.sum(bbox_vals) > 0:
        verify_flag = True
    
    mask_file_ds = None
    
    return verify_flag
    

def create_memory_raster_obj_without_data(source_raster_gdal_obj, name_string = ''):
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromWkt(source_raster_gdal_obj.GetProjectionRef())
    projection = outRasterSRS.ExportToWkt()
    
    geotransform = source_raster_gdal_obj.GetGeoTransform()
    
    outputDataType = source_raster_gdal_obj.GetRasterBand(1).DataType
    
    width = source_raster_gdal_obj.RasterXSize
    height = source_raster_gdal_obj.RasterYSize
    nbands = source_raster_gdal_obj.RasterCount
    
    # create the output tif
    driver = gdal.GetDriverByName ('MEM')
    dstGdalObj = driver.Create (name_string, width, height, nbands, outputDataType)
    
    dstGdalObj.SetGeoTransform (geotransform)
    
    """
    outputNoDataValue = source_raster_gdal_obj.GetRasterBand(1).GetNoDataValue()
    for i,outBandArray in enumerate(outBandsAsArray):
        dstGdalObj.GetRasterBand(i + 1).WriteArray (outBandArray)
        
        if outputNoDataValue is not None:
            dstGdalObj.GetRasterBand(i + 1).SetNoDataValue(outputNoDataValue)
            
        dstGdalObj.GetRasterBand(i + 1).FlushCache()
    """
    
    dstGdalObj.SetProjection (projection)
        
    # copy all meta data:
    dstGdalObj.SetMetadata (source_raster_gdal_obj.GetMetadata(""), "")
    for domain in ["RPC", "IMAGE_STRUCTURE", "SUBDATASETS", "GEOLOCATION"]:
        dstGdalObj.SetMetadata (source_raster_gdal_obj.GetMetadata(domain), domain)
    dstGdalObj.SetGCPs (source_raster_gdal_obj.GetGCPs(), source_raster_gdal_obj.GetGCPProjection())
    
    return dstGdalObj  

def update_rpcs_from_rpb(raster_object, new_rpc_file):
    '''
    Function to update rpc metadata of a gdal object 
    using a new rpc file 
    '''
    rpc_metadata = raster_object.GetMetadata('RPC')   
    if len(rpc_metadata.keys()) == 0:
        #print(" No rpc metadata found. Using standard dict\n")
        
        rpc_metadata = { 'HEIGHT_OFF': None,
                         'SAMP_OFF': None,
                         'LINE_NUM_COEFF' : None,
                         'LONG_OFF' : None,
                         #'MIN_LAT' : None,
                         #'MAX_LONG' : None,
                         'LINE_SCALE' : None,
                         'SAMP_NUM_COEFF' : None,
                         'LONG_SCALE' : None,
                         'SAMP_DEN_COEFF' : None,
                         #'MIN_LONG' : None,
                         'SAMP_SCALE' : None,
                         #'MAX_LAT' : None,
                         'LAT_SCALE' : None,
                         'LAT_OFF' : None,
                         'LINE_OFF' : None,
                         'LINE_DEN_COEFF' : None,
                         'HEIGHT_SCALE' : None}

    
    check_if_exists(new_rpc_file)
    
    new_rpcs = load_rpb(new_rpc_file)
    
    for key in rpc_metadata:
        # RAPDR rpc class does not bother with these following keys. So it is ok if they are missing
        # in new_rpc_file
        if key in ['MIN_LAT', 'MAX_LONG', 'MIN_LONG', 'MAX_LAT']:
            if key not in new_rpcs:
                continue
                    
        x = new_rpcs[key]
        if isinstance(x, np.ndarray):
            x = [str(i) for i in x]
            x = ' '.join(x)
        else:
            x = str(x)
        rpc_metadata[key] = x

    for key in rpc_metadata:
        assert rpc_metadata[key] is not None
          
    for domain in ["RPC"]:
        raster_object.SetMetadata(rpc_metadata, domain)     
    
    new_rpcs = None
    return raster_object

def get_epsg_code_from_lon_lat(lon, lat):
    
    utm_zone = np.int(np.floor( 31.0 + 1.0*lon/6.0 ))
    
    if lat >= 0:
        epsg_utm_code = "EPSG:" + str(32600 + utm_zone)
    else:
        epsg_utm_code = "EPSG:" + str(32700 + utm_zone)
    
    return epsg_utm_code
        
def get_epsg_code_of_raster(input_gtiff):
    
    srs_ds = gdal.Open(input_gtiff)
    
    mid_pt_pix = [[(srs_ds.RasterXSize/2) + 0.5, (srs_ds.RasterYSize/2) + 0.5]]
    
    mid_pt_lon_lat = convert_pixel_coords_to_projected_coords(mid_pt_pix, srs_ds.GetGeoTransform())
    
    return get_epsg_code_from_lon_lat(mid_pt_lon_lat[0,0], mid_pt_lon_lat[0,1])


def bilinear_interpolation(input_array_coords, mask_rows, input_array, no_data_value = None, pixelis = 'area'):
    
    # Bilinear interpolation for image
    x = input_array_coords[mask_rows,0]
    y = input_array_coords[mask_rows,1]

    """
    # In this case pixel centers represent pixel, so
    # weights must be measured relative to these
    
    * Turns out the near and iinterpolated image look shifted
    * if I have an extra halfshift here
    """
    """
    if pixelis == "area":
        
        halfshift = 0.5
        x = x - halfshift
        y = y - halfshift
    """
    x_min = 0
    y_min = 0
    x0 = np.floor(x).astype(int)
    y0 = np.floor(y).astype(int)
    x1 = x0 + 1
    y1 = y0 + 1
    w_x = x - x0
    w_y = y - y0
    
    x0[x0 < 0] = x1[x0 < 0]
    y0[y0 < 0] = y1[y0 < 0]

    x1[x1 >= input_array.shape[2] ] = x0[x1 >= input_array.shape[2]]
    y1[y1 >= input_array.shape[1] ] = y0[y1 >= input_array.shape[1]]
            
    img00 = input_array[:, y0-y_min, x0-x_min]
    img01 = input_array[:, y0-y_min, x1-x_min]
    img10 = input_array[:, y1-y_min, x0-x_min]
    img11 = input_array[:, y1-y_min, x1-x_min]
    
    # If we have no data values
    if no_data_value is not None:
        
        no_data_mask = np.ones_like(img00, dtype = np.bool)
        # Find points where all four coords used for bilinear interpolation have no data
        no_data_mask = np.logical_and(no_data_mask, (img00 == no_data_value))
        no_data_mask = np.logical_and(no_data_mask, (img01 == no_data_value))
        no_data_mask = np.logical_and(no_data_mask, (img10 == no_data_value))
        no_data_mask = np.logical_and(no_data_mask, (img11 == no_data_value))
        
        # No data is implicitl set to 0 while interpolating
        output = (img00 != no_data_value)*img00*(1-w_x)*(1-w_y) +\
                 (img01 != no_data_value)*img01*w_x *(1-w_y) +\
                 (img10 != no_data_value)*img10*(1-w_x)*w_y  +\
                 (img11 != no_data_value)*img11*w_x *w_y
        
        # Set points where all four surrounding values are no data         
        output[no_data_mask] = no_data_value
        
    else:
        output = img00*(1-w_x)*(1-w_y) +\
                 img01*w_x *(1-w_y) +\
                 img10*(1-w_x)*w_y  +\
                 img11*w_x *w_y


    return output
    
def nearest_interpolation(input_array_coords, mask_rows, input_array):
    
    # Nearest interpolation for image
    x = input_array_coords[mask_rows,0]
    y = input_array_coords[mask_rows,1]

    x0 = np.round(x).astype(int)
    y0 = np.round(y).astype(int)
    
    output = input_array[:, y0, x0]
        
    
    return output, y0, x0

def writeGTiffImage (projection, geotransform, outBandsAsArray, outputRasterName, outputDataType, outputNoDataValue = None, compress_flag = False):
    
    """ create an output geotiff file using explicitly provided metadata params """
    
    if os.path.isfile(outputRasterName):
        os.remove(outputRasterName)  
        #command = "rm -rf " + outputRasterName;
        #os.system(command);
    
    
    (height, width) = np.shape(outBandsAsArray[0])
    nbands = len(outBandsAsArray)
    # create the output tif
    driver = gdal.GetDriverByName ('GTiff')
    if compress_flag:
        dstGdalObj = driver.Create (outputRasterName, width, height, nbands, outputDataType, options=["TILED=YES", "COMPRESS=LZW", "BIGTIFF=YES"])
    else:
        dstGdalObj = driver.Create (outputRasterName, width, height, nbands, outputDataType, options=["TILED=YES", "BIGTIFF=YES"])
    
    dstGdalObj.SetGeoTransform (geotransform)
    
    for i,outBandArray in enumerate(outBandsAsArray):
        dstGdalObj.GetRasterBand(i + 1).WriteArray (outBandArray)
        
        if outputNoDataValue is not None:
            dstGdalObj.GetRasterBand(i + 1).SetNoDataValue(outputNoDataValue)
            
        dstGdalObj.GetRasterBand(i + 1).FlushCache()
    
    dstGdalObj.SetProjection (projection)
    
    dstGdalObj = None      


def convert_projected_coords_to_pixel_coords(projected_x_y_array, geotransform):
        
    pixel_col_row_array = projected_x_y_array - np.array((geotransform[0], geotransform[3]))
    
    pixel_col_row_array[:,0] = pixel_col_row_array[:,0]/geotransform[1]
    
    pixel_col_row_array[:,1] = pixel_col_row_array[:,1]/geotransform[5]

    return pixel_col_row_array
       

def convert_pixel_coords_to_projected_coords(pixel_col_row_list, geotransform):
    
    geotr_mat = np.reshape(geotransform, (2,3))
    
    n_pixels = len(pixel_col_row_list) 
    
    return np.transpose(np.dot(geotr_mat, np.concatenate([ np.ones((1,n_pixels)), np.transpose(np.array(pixel_col_row_list)) ])))
        
#  aa

def check_if_exists(input_file_or_folder):
    if os.path.isfile(input_file_or_folder) or os.path.isdir(input_file_or_folder):
        return True
    else:
        raise ValueError("\n\nERROR: " + input_file_or_folder + " does not exist. Stopping \n")
        return False


def get_ortho_grid_worker(input_list):
    """
    function that implements gwarp++ for a block. A block could also cover an entire raster if it is small enough.
    By block we refer to a grid in lon lat coordinates
    """
        
    # Encapsulate within a try catch statement to detect blocks that fail
    try:
        #print(input_list)        
        outputRasterName, start_lon, start_lat, end_lon, end_lat, spacing_lon, spacing_lat, \
                raster_file, dem_file, dem_interpolate_method, image_interpolate_method, dst_nodata, block_idx, \
                save_mapping_flag, border_block_flag, rpc_file, compress_flag, dem_mask_file, dem_low_res_file = input_list
    
        # verify if we need mask file
        if dem_mask_file is not None and str(dem_mask_file) != "None":
            verify_flag = verify_if_mask_file_is_needed(start_lon, start_lat, end_lon, end_lat, spacing_lon, spacing_lat, dem_mask_file)
            if not verify_flag:
                print("\nDEM MASK NOT REQUIRED\n")
                dem_mask_file = None    
    
        command = codePath + "/c++/gwarp++ --output %s --start_lon %s --start_lat %s --end_lon %s --end_lat %s --spacing_lon %s --spacing_lat %s " + \
                  " --raster_file %s --dem_file %s --dem_interpolate_method %s --image_interpolate_method %s  --dst_nodata %s --block_idx %s " + \
                  " --save_mapping_flag %s --border_block_flag %s --rpc_file %s --compress_flag %s --dem_mask_file %s --dtm_file %s"
        # note block_idx is set to 0 for worker cpp function
        command = (command % (outputRasterName, start_lon, start_lat, end_lon, end_lat, spacing_lon, spacing_lat, \
                              raster_file, dem_file, dem_interpolate_method, image_interpolate_method, dst_nodata, block_idx, \
                              save_mapping_flag, border_block_flag, rpc_file, compress_flag, str(dem_mask_file), dem_low_res_file  ))
        
        print(command)
    
        os.system(command)
        
    except Exception as e:
        
        print("ID is " + str(block_idx))
        traceback.print_exc()
            
    return 
        
def get_ortho_grid(outputRasterName, utm_code, raster_file, dem_file, dem_low_res_file, dem_interpolate_method, image_interpolate_method, parallel_flag,
                   output_resolution = None, n_workers = 12, gdal_merge = None, dst_nodata = -9999, save_mapping_flag = False, rpc_file = None,
                   compress_flag = False, dem_mask_file = None):
    halfshift = 0.5
    """
    function that prepares blocks and grids for orthorectification
    """
    
    if parallel_flag and gdal_merge is None:
        raise Exception('When using parallel, need to supply gdal_merge.py path')

    # Check if input raster has rpc metadata
    raster_object = gdal.Open(raster_file)
    
    # Create rpc based transformer
    if rpc_file is not None:
        check_if_exists(rpc_file)
        raster_obj_rpc = create_memory_raster_obj_without_data(raster_object) 
        raster_obj_rpc = update_rpcs_from_rpb(raster_obj_rpc, rpc_file)
    else:
        raster_obj_rpc = raster_object
    
    rpc_metadata = raster_obj_rpc.GetMetadata('RPC')
    if len(rpc_metadata.keys()) == 0:
        raise ValueError(" ERROR: No rpc metadata found\n")

    rpc_dem_missing_value = rpc_metadata['HEIGHT_OFF']
    
    
    # default error is 0 pixel
    # Specify dem missing value as 0. Otherwise due to a bug in gdal, if the dem has Nans, then
    # the code will fail with errors. This is true for standard gdalwarp as well
    rpc_transformer_options_list = ['METHOD=RPC','RPC_PIXEL_ERROR_THRESHOLD=0.001','RPC_DEM=%s' % dem_file, 
                                    'RPC_DEMINTERPOLATION=%s' % dem_interpolate_method,
                                    'RPC_DEM_MISSING_VALUE='+rpc_dem_missing_value]
    
    # rpc_transformer_options_list = ['METHOD=RPC','RPC_PIXEL_ERROR_THRESHOLD=0.001']
    # Create rpc based transformer
    rpc_transformer = gdal.Transformer(raster_obj_rpc, None, rpc_transformer_options_list)
    
    # Get raster width and height
    raster_w = raster_object.RasterXSize
    raster_h = raster_object.RasterYSize 
    
    # Corners of raster in pixel space
    raster_corners_col_row = [(halfshift,halfshift), (raster_w + halfshift, halfshift), (halfshift, raster_h + halfshift), (raster_w + halfshift, raster_h + halfshift)]
    
    # Corners of raster in orthorectified space
    raster_corners_lon_lat, flags  = rpc_transformer.TransformPoints(0, raster_corners_col_row)
    
    # Check if all four corners were correctly found. Most common reason for failure is missing DEM value.
    # In this case there will be a 0 in the flags
    rpc_error = 0.001
    
    while (0 in flags):
        
        if rpc_error >= 0.1:
            break
            
        rpc_error = rpc_error + 0.005
        rpc_transformer_options_list = ['METHOD=RPC','RPC_PIXEL_ERROR_THRESHOLD=' + str(rpc_error),'RPC_DEM=%s' % dem_file, 
                                        'RPC_DEMINTERPOLATION=%s' % dem_interpolate_method,
                                        'RPC_DEM_MISSING_VALUE='+rpc_dem_missing_value]            
        # Create rpc based transformer
        rpc_transformer = gdal.Transformer(raster_obj_rpc, None, rpc_transformer_options_list)
        
        # Corners of raster in orthorectified space
        raster_corners_lon_lat, flags  = rpc_transformer.TransformPoints(0, raster_corners_col_row)
        
    if 0 in flags:
        rpc_error = -0.004
        print("Trying low res dem\n")
        
    while (0 in flags):
                
        if rpc_error >= 0.1:
            break
            
        rpc_error = rpc_error + 0.005
        
        rpc_transformer_options_list = ['METHOD=RPC','RPC_PIXEL_ERROR_THRESHOLD=' + str(rpc_error),'RPC_DEM=%s' % dem_low_res_file, 
                                        'RPC_DEMINTERPOLATION=%s' % dem_interpolate_method,
                                        'RPC_DEM_MISSING_VALUE='+rpc_dem_missing_value]            
        # Create rpc based transformer
        rpc_transformer = gdal.Transformer(raster_obj_rpc, None, rpc_transformer_options_list)
        
        # Corners of raster in orthorectified space
        raster_corners_lon_lat, flags  = rpc_transformer.TransformPoints(0, raster_corners_col_row)
        
    if 0 in flags:   
        
        rpc_transformer_options_list = ['METHOD=RPC','RPC_PIXEL_ERROR_THRESHOLD=0.001']
        
        print("NOT USING DEM FOR CORNERS\n")
        
        # Create rpc based transformer
        rpc_transformer = gdal.Transformer(raster_obj_rpc, None, rpc_transformer_options_list)
        
        # Corners of raster in orthorectified space
        raster_corners_lon_lat, flags  = rpc_transformer.TransformPoints(0, raster_corners_col_row)
    
    else:
        print ("Final rpc error ", rpc_error)
        
    if 0 in flags:
        print("\nRaster corners in pixel are ", raster_corners_col_row)
        print("\nRaster corners in lon lat are ", raster_corners_lon_lat)
        raise ValueError("\n\nERROR: Unable to get height from dem at bounds.\n")
        
    raster_object = None
    raster_obj_rpc = None
    
    raster_corners_lon = [raster_corners_lon_lat[i][0] for i in range(4)]
    raster_corners_lat = [raster_corners_lon_lat[i][1] for i in range(4)]
    
    min_lon = np.min(raster_corners_lon)
    max_lon = np.max(raster_corners_lon)
    
    min_lat = np.min(raster_corners_lat)
    max_lat = np.max(raster_corners_lat)
    
    raster_corners_lon = [min_lon, max_lon, min_lon, max_lon]
    raster_corners_lat = [max_lat, max_lat, min_lat, min_lat]
    
    # If the output resolution is not specified in metres, then it is set to 0.5*the resolution of the epsg4326 dem
    # i.e in degrees
    if output_resolution is None:
        dem_object = gdal.Open(dem_file)
        dem_geotransform  = dem_object.GetGeoTransform()
        spacing_lon,  spacing_lat = dem_geotransform[1], dem_geotransform[5]
        dem_object = None 
        scale = 0.5
        spacing_lon = scale*spacing_lon
        spacing_lat = scale*spacing_lat
    else:
        # If the output resolution is specified in metres, then depending upon the location, the resolution in degrees can vary.
        # Hence we first project the corner of the raster into utm, find the next point. Project it back into lon lat. 
        # The difference  gives us the resolution in degrees
        output_resolution_x = output_resolution
        # y coordinate is flipped. Hence the negative sign.
        output_resolution_y = -output_resolution
        epsg4326_proj = pyproj.Proj(init = 'EPSG:4326')
        utm_proj = pyproj.Proj(init = utm_code)
        # Find raster corners in projected coordinates
        raster_corners_proj_col, raster_corners_proj_row = pyproj.transform(epsg4326_proj, utm_proj, raster_corners_lon, raster_corners_lat)
        # Find the next point by stepping output resolution metres along x and y axis
        next_point_grid_lon, next_point_grid_lat = pyproj.transform(utm_proj, epsg4326_proj,
                                                   raster_corners_proj_col[0] + output_resolution_x,
                                                   raster_corners_proj_row[0] + output_resolution_y)
        
        
        spacing_lon = 1.0*(next_point_grid_lon - raster_corners_lon[0])
        spacing_lat = 1.0*(next_point_grid_lat - raster_corners_lat[0])
        
    print("\nCorners of the output raster are ")
    pp.pprint(list(zip(raster_corners_lon, raster_corners_lat)))
    print("\nOutput resolution in degrees is %s, %s\n" % (spacing_lon, spacing_lat))
    
    if parallel_flag:
        # We split the lon lat grid into blocks. The number of points in each block is 2000x2000. In practice
        # we add padding to the block to account for height effects
        n_pts_block = 2000
        
        # Starting lon lat for each block
        blocks_lon = np.arange(raster_corners_lon[0], raster_corners_lon[-1], spacing_lon*n_pts_block)
        blocks_lat = np.arange(raster_corners_lat[0], raster_corners_lat[-1], spacing_lat*n_pts_block)

        # Check distance between last pt in arange and actual end
        dist_lon_end = (raster_corners_lon[-1] - blocks_lon[-1])/(1.0*spacing_lon*n_pts_block)
        dist_lat_end = (raster_corners_lat[-1] - blocks_lat[-1])/(1.0*spacing_lat*n_pts_block)

        # If the last sample pt is close to actual end, just replace it.
        # If it is not, add the actual end to the list
        # By close, we mean if its distance to actual end < block_size/4
        if dist_lon_end < n_pts_block/4.0:
            blocks_lon[-1] = raster_corners_lon[-1]
        else:
            blocks_lon.append(raster_corners_lon[-1])
            
        if dist_lat_end < n_pts_block/4.0:
            blocks_lat[-1] = raster_corners_lat[-1]
        else:
            blocks_lat.append(raster_corners_lat[-1])
            
        input_mp_list = []
        
        block_idx = 0
        
        outputRasterFolder, outputRasterName_only = os.path.split(outputRasterName)
        if outputRasterFolder == '':
            outputRasterFolder = '.'
        
        outputRasterFolder_workers = os.path.join(outputRasterFolder, 'worker',  os.path.splitext(outputRasterName_only)[0])
        if not os.path.isdir(outputRasterFolder_workers):
            os.makedirs(outputRasterFolder_workers)
        
        # For loop to get the start and end coordinates for each block (without padding)
        for lon_index,lon in enumerate(blocks_lon):
            
            for lat_index,lat in enumerate(blocks_lat):
            
                if (lon_index == blocks_lon.size - 1) or (lat_index == blocks_lat.size - 1):
                    continue
                
                border_block_flag = False
                
                if (lon_index == blocks_lon.size -2) or (lat_index == blocks_lat.size - 2):
                    border_block_flag = True
                                                
                start_lon = lon
                start_lat = lat
                                
                end_lon = blocks_lon[lon_index + 1]
                end_lat = blocks_lat[lat_index + 1]
                    
                # name the output raster based on block index                
                outputRasterName_worker = outputRasterName_only.replace('.tif', '_' + str(block_idx) + '.tif')
                outputRasterName_worker = os.path.join(outputRasterFolder_workers, outputRasterName_worker)
                blk_list = [outputRasterName_worker, start_lon, start_lat, end_lon, end_lat, spacing_lon, spacing_lat, raster_file,
                            dem_file, dem_interpolate_method, image_interpolate_method, dst_nodata, block_idx, save_mapping_flag,
                            border_block_flag, rpc_file, compress_flag, dem_mask_file, dem_low_res_file]
                input_mp_list.append(blk_list)
                
                block_idx += 1
        
        print("\nTotal number of blocks is %s" % block_idx)
        
        # Launch pool
        pool_workers = Pool(n_workers)
        
        pool_workers.map(get_ortho_grid_worker, input_mp_list)
        
        #with open('debug.pickle','wb') as f:
        #    pickle.dump(input_mp_list, f)
        pool_workers.close()
        
        tmpFileList = os.path.join(outputRasterFolder_workers, "tmpFileList.txt");
        f = open(tmpFileList, 'w')
        
        if f is None:
            print('Error opening tmp file {0} for writing'.format(tmpFileList))
            return None
    
        tileCount = 0
                
        for file in os.listdir(outputRasterFolder_workers):
            if file.endswith(".tif"):
                # skip coordmap tifs
                if not file.endswith("_coord_map.tif"):
                    f.write(os.path.join(outputRasterFolder_workers, file) + '\n')
                tileCount = tileCount + 1
                
        f.close()
        
        my_env = os.environ.copy()
        
        outputRasterVrt = outputRasterName.replace('.tif', '.vrt')
        command = "%s/gdalbuildvrt %s -tr %s %s -srcnodata %i -input_file_list %s"%(gdal_merge, outputRasterVrt, np.abs(spacing_lon), np.abs(spacing_lat), dst_nodata, tmpFileList)    
        p = subprocess.Popen(command,  shell=True, env=my_env)
        retval = p.wait()
        
        if compress_flag:
            #command = "%s -o %s -co TILED=YES -co BIGTIFF=YES -co COMPRESS=LZW -a_nodata %i --optfile %s"%(gdal_merge, outputRasterName, dst_nodata, tmpFileList)
            command = "%s/gdal_translate  -co TILED=YES -co BIGTIFF=YES -co COMPRESS=LZW -a_nodata %i %s %s"%(gdal_merge, dst_nodata, outputRasterVrt, outputRasterName)
        else:
            #command = "%s -o %s -co TILED=YES -co BIGTIFF=YES -a_nodata %i --optfile %s"%(gdal_merge, outputRasterName, dst_nodata, tmpFileList)
            command = "%s/gdal_translate -co TILED=YES -co BIGTIFF=YES -a_nodata %i %s %s"%(gdal_merge, dst_nodata, outputRasterVrt, outputRasterName)
            
        p = subprocess.Popen(command,  shell=True, env=my_env)
        retval = p.wait()
        
        # save pixel coords map between ortho image and raw image
        if save_mapping_flag:
            
            # First merge all the  worker 2 band rasters with the same metadata as the orthorectified image
            # First band will contain the column index of the corresponding raw image pixel
            # Second band will contain the row index of the corresponding raw image pixel
            outputCoordMapRasterName = outputRasterName.replace('.tif', '_coord_map.tif')
            
            tmpFileList_coordmap = os.path.join(outputRasterFolder_workers, "tmpFileList_coordmap.txt");
            f = open(tmpFileList_coordmap, 'w')
            
            if f is None:
                print('Error opening tmp file {0} for writing'.format(tmpFileList_coordmap))
                return None
        
            tileCount = 0
                    
            for file in os.listdir(outputRasterFolder_workers):
                if file.endswith("_coord_map.tif"):
                    f.write(os.path.join(outputRasterFolder_workers, file) + '\n')
                    tileCount = tileCount + 1
                    
            f.close()
            
            outputCoordMapRasterVrt = outputCoordMapRasterName.replace('.tif', '.vrt')
            command = "%s/gdalbuildvrt %s -tr %s %s -srcnodata %i -input_file_list %s"%(gdal_merge, outputCoordMapRasterVrt, np.abs(spacing_lon), np.abs(spacing_lat), dst_nodata, tmpFileList_coordmap)    
            p = subprocess.Popen(command,  shell=True, env=my_env)
            retval = p.wait()
                        
            if compress_flag:
                """
                command = "%s -o %s -co TILED=YES -co BIGTIFF=YES -co COMPRESS=LZW --optfile %s -a_nodata %i -init %i "%(gdal_merge, 
                                                                                                                         outputCoordMapRasterName, 
                                                                                                                         tmpFileList_coordmap, 
                                                                                                                         dst_nodata, dst_nodata)
                """
                command = "%s/gdal_translate  -co TILED=YES -co BIGTIFF=YES -co COMPRESS=LZW -a_nodata %i %s %s"%(gdal_merge, dst_nodata, outputCoordMapRasterVrt, outputCoordMapRasterName)
            else:
                """
                #command = "%s -o %s -co TILED=YES -co BIGTIFF=YES --optfile %s -a_nodata %i -init %i "%(gdal_merge, 
                                                                                                        outputCoordMapRasterName, 
                                                                                                        tmpFileList_coordmap, 
                                                                                                        dst_nodata, dst_nodata)
                """
                command = "%s/gdal_translate -co TILED=YES -co BIGTIFF=YES -a_nodata %i %s %s"%(gdal_merge, dst_nodata, outputCoordMapRasterVrt, outputCoordMapRasterName)
                
            
            q = subprocess.Popen(command,  shell=True, env=my_env)
            retval_q = q.wait()
        
        shutil.rmtree(outputRasterFolder_workers)
        #command = 'rm -rf ' + outputRasterFolder_workers
        #os.system(command)
    
    else:
        # Ensure that the image is small enough. Otherwise we will run into memory issues.
        start_lon = raster_corners_lon[0]
        start_lat = raster_corners_lat[0]
        
        end_lon = raster_corners_lon[-1]
        end_lat = raster_corners_lat[-1]
                    
        input_mp_list = [outputRasterName, start_lon, start_lat, end_lon, end_lat, spacing_lon, spacing_lat, raster_file,
                         dem_file, dem_interpolate_method, image_interpolate_method, dst_nodata, 0, save_mapping_flag, True,
                         rpc_file, compress_flag, dem_mask_file, dem_low_res_file]
        get_ortho_grid_worker(input_mp_list)
    
    
    # Convert coordmap tif to npy file for easy use. Only for small images. NOT FOR BIG IMAGES. 
    # Will take about 30GB of memory for huge images as locations can be 32 bit or 64 bit
    # unlike images which can be 16 bit.  
    """
    if save_mapping_flag:

        outputCoordMapRasterName = outputRasterName.replace('.tif', '_coord_map.tif')
        if not os.path.isfile(outputCoordMapRasterName):
            print("ERROR:Coord map file not found.Stopping\n")
        else:
            coordmap_object = gdal.Open(outputCoordMapRasterName)
            coordmap_raster_col = coordmap_object.GetRasterBand(1).ReadAsArray().flatten()
            coordmap_raster_row = coordmap_object.GetRasterBand(2).ReadAsArray().flatten()
    
            # Append the ortho image width and height as the last column
            '''
            output_coord_map_file_name = outputCoordMapRasterName.replace('.tif', '.csv')
            np.savetxt(output_coord_map_file_name, 
                       np.vstack((np.transpose([coordmap_raster_col,
                                              coordmap_raster_row]),
                                  np.array([coordmap_object.RasterXSize,
                                            coordmap_object.RasterYSize])[None,:] )), delimiter = ',', fmt="%d")
        
            '''
            output_coord_map_file_name = outputCoordMapRasterName.replace('.tif', '.npy')
            np.save(output_coord_map_file_name, 
                    np.vstack((np.transpose([coordmap_raster_col,
                                              coordmap_raster_row]),
                               np.array([coordmap_object.RasterXSize,
                                         coordmap_object.RasterYSize])[None,:] )))
            
            # Delete the coordmap raster        
            #os.remove(outputCoordMapRasterName)
    """

def check_input(args):

    check_if_exists(args.input_raster)
    
    if args.RPC_DEM is not None:
        check_if_exists(args.RPC_DEM)   
        
    if args.dem_interp not in ['near', 'bilinear', 'cubic']:
        raise ValueError("\n\nERROR:DEM interpolation should be one of near, bilinear or cubic. You have provided %s\n" % args.dem_interp)

    if args.image_interp not in ['near', 'bilinear']:
        raise ValueError("\n\nERROR:Image interpolation should be one of near or bilinear. You have provided %s\n" % args.image_interp)
    
    
def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    
    parser = argparse.ArgumentParser(prog = 'gwarp++', description = 'gwarp to orthorectify images using high resolution dsms')
    
    parser.add_argument('input_raster', help = 'ntf or tif to orthorectify')
    
    parser.add_argument('-RPC_DEM', help = 'optional dem or dsm')
    
    parser.add_argument('-low_res_dem', help = 'low res dem wrt to wgs84 ellipsoid')
    
    parser.add_argument('-utm', help = 'flag to specify if dem is in utm. If so a epsg4326 dem is created from it', action='store_true')
    
    parser.add_argument('-dem_interp', help = 'dem interpolation - one of near, bilinear or cubic', default = 'bilinear')

    parser.add_argument('-image_interp', help = 'raster interpolation - one of near or bilinear', default = 'bilinear')
    
    parser.add_argument('-output_res', type=float , help = 'resolution of orthorectified image in metres', default = None)
    
    parser.add_argument('-parallel', help = 'flag for parallel processing', action='store_true')
    
    parser.add_argument('-n_workers', type=int , help = 'number of workers for parallel processing', default = 12)
    
    parser.add_argument('-output', type=str , help = 'optional name of output raster', default = None)

    parser.add_argument('-gdal_merge', type=str , help = 'path to gdal_merge.py', default = None)

    parser.add_argument('-dst_nodata', type=int, help='no data value for destination', default=-9999)

    parser.add_argument('-save_mapping_flag', help = 'flag to save pixel mapping between raw and orthorectified image', action='store_true')
    
    parser.add_argument('-rpc_file', help = 'optional rpc file with updated bias values', type = str, default = None)
    
    parser.add_argument('-compress_flag', help = 'flag to turn on compression for output image', action='store_true')
    
    parser.add_argument('-dem_mask_file', help = 'uint8 bit mask file for dem indicating points to skip', type = str, default = None)
    
    args = parser.parse_args(argv)
    
    check_input(args)
    
    if args.parallel and args.gdal_merge is None:
        raise Exception('When using parallel, need to supply gdal_merge.py path')
    '''
    create temporary folder 
    
    temp_folder = '/tmp/gwarp/'
    if os.path.isdir(temp_folder):
        shutil.rmtree(temp_folder)
    
    os.makedirs(temp_folder)
    if not os.path.isdir(temp_folder):
        raise ValueError("\n\nERROR: Could not create %s \n" % temp_folder)
    '''
    '''
    Output raster name
    '''
    if args.output is None :
        
        # We add _ortho at the end of the name of the input tif. We create a folder called orthorectified in the same folder as the input image
        # and save the output in that folder
        extension = os.path.splitext(args.input_raster)[-1]
        outputRasterName = args.input_raster.replace(extension, '_ortho.tif')
        
        outputRasterFolder = os.path.join(os.path.split(outputRasterName)[0], 'orthorectified')
        if not os.path.isdir(outputRasterFolder):
            os.makedirs(outputRasterFolder)
            
        outputRasterName = os.path.join(outputRasterFolder, os.path.split(outputRasterName)[-1])
        
    else:
        # Use specified name for output
        outputRasterName = os.path.abspath(args.output)
        # Create output folder if needed
        outputRasterFolder = os.path.split(outputRasterName)[0]
        if not os.path.isdir(outputRasterFolder):
            os.makedirs(outputRasterFolder)
    
    dem_epsg4326 = None
    # If DEM is none, then use standard gdalwarp
    if args.RPC_DEM is None:    
        command =   "set GDAL_CACHEMAX=3000;set GDAL_SWATH_SIZE=3073741824;set GDAL_NUM_THREADS=2; " \
                    "gdalwarp -to RPC_DEM_MISSING_VALUE=0 -co COMPRESS=PACKBITS -co INTERLEAVE=BAND -co BIGTIFF=YES " \
                    "-co TILED=YES -rpc -et 0 -overwrite -wo WRITE_FLUSH=YES " \
                    "-wo OPTIMIZE_SIZE=YES -multi -r cubic -dstnodata %i -wm 3000 -of GTIFF -t_srs EPSG:4326 " % args.dst_nodata + args.input_raster + " " + outputRasterName

        os.system(command)
    
    else:   
        
        dem_object = gdal.Open(args.RPC_DEM)
        dem_srs= osr.SpatialReference(wkt=dem_object.GetProjectionRef())  
        dem_projcs = dem_srs.GetAttrValue('PROJCS')
            
        if args.utm: 
            ''' 
            convert dem or dsm from utm to epsg4326 i.e. longitude, latitude 
            ''' 
            # Get dem projcs. i.e. utm code
            
            dem_srs_epsg = dem_srs.GetAttrValue('AUTHORITY', 0) + ':' + dem_srs.GetAttrValue('AUTHORITY', 1)
            
            # Name of dem epsg4326 tif
            dem_ext = os.path.splitext(args.RPC_DEM)[-1]
            dem_epsg4326 = args.RPC_DEM.replace(dem_ext, '_epsg4326.tif')
            
            if dem_projcs is not None and 'UTM' in dem_projcs:    
                if not os.path.isfile(dem_epsg4326):
                    # warp dem from utm to epsg 4326 if not done already
                    print("Converting dem from " + dem_srs_epsg + " to EPSG:4326" )
                    
                    dem_nodataval = dem_object.GetRasterBand(1).GetNoDataValue()
                    
                    if dem_nodataval is None:
                        dem_nodataval = np.nan 
                    
                    dem_nodataval = str(dem_nodataval)
                    
                    command = 'gdalwarp -t_srs EPSG:4326 -r near -srcnodata ' + dem_nodataval +  ' -dstnodata ' + dem_nodataval  + ' ' + args.RPC_DEM + ' ' + dem_epsg4326
                    os.system(command)
            else:
                raise ValueError("\nAre you sure if DEM is in UTM projection? Stopping")
            
        else:
            if dem_projcs is not None and 'UTM' in dem_projcs:
                raise ValueError("UTM Dem provided. Use the -utm flag\n")
            
            dem_srs_epsg = get_epsg_code_of_raster(args.RPC_DEM)
            dem_epsg4326 = args.RPC_DEM

        dem_object = None 
        
        get_ortho_grid(outputRasterName, dem_srs_epsg, args.input_raster, 
                       dem_epsg4326, args.low_res_dem, args.dem_interp, args.image_interp, args.parallel,
                       args.output_res, args.n_workers, args.gdal_merge, args.dst_nodata, args.save_mapping_flag,
                       args.rpc_file, args.compress_flag, args.dem_mask_file)
        
        
if __name__ == "__main__":
    sys.exit(main())
