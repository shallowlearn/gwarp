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

import importlib
gwarp = importlib.import_module('gwarp++')
import glob
import os, sys, shutil
from utilities import remove_trailing_slash,remove_folder_ext, \
                        create_directory, check_if_exists, remove_folder_from_name, get_folder
                                
from multiprocessing import Pool

#from pprint import pprint

import argparse

def call_ortho(input_list):
    
    #pp.pprint(input_list)
    """
    If the chip and RPB file are in different folders,
    then c++ might not find the RPC metadata. So create
    tmp folder with symbolic links
    """
    tmp_dir = input_list[0]
    rpc_name = input_list[-2]
    img_name = input_list[-1]
    output_name = input_list[-4]
    just_name = remove_folder_ext(img_name)
    tmp_dir = tmp_dir + '/' + just_name + '/'
    create_directory(tmp_dir)
    create_directory(tmp_dir + '/orthorectified')
        
    new_img_name = tmp_dir + '/' + remove_folder_from_name(img_name)
    if rpc_name:
        new_rpc_name = tmp_dir + '/' + remove_folder_from_name(rpc_name)
    else:
        new_rpc_name = None
    
    # Save output locally first
    new_output_name = tmp_dir + '/orthorectified/' + remove_folder_from_name(output_name)
    # create symbolic links
    os.symlink(os.path.abspath(img_name), new_img_name)
    print("\nSymbolinking ", img_name, new_img_name)
    
    if rpc_name:
        os.symlink(os.path.abspath(rpc_name), new_rpc_name)
        print("\nSymbolinking ", rpc_name, new_rpc_name)
    
    # Remove cache dir,rpc and img from input list
    copy_list = [i for i in input_list[1:-2]]
    if new_rpc_name:
        copy_list.extend([new_rpc_name])
        opidx = -4
    else:
        copy_list = copy_list[:-1]
        opidx = -2
    copy_list.extend([new_img_name])
    copy_list[opidx] = new_output_name
    
    #pprint(copy_list)
    # call orthorect
    gwarp.main(copy_list)
    
    # move temp output back to original
    shutil.move(new_output_name, output_name)
    
    # if coordmap was saved move it also to original
    new_output_coordmap_name = new_output_name.replace('.tif', '_coord_map.tif')
    if os.path.isfile(new_output_coordmap_name):
        shutil.move(new_output_coordmap_name, output_name.replace('.tif', '_coord_map.tif'))
    return

def batch_orthorectify(folder_img, folder_rpc, dsm_file, dem_file,
                       folder_output, CACHE_DIR, parallel_flag = False, num_workers = 1, dsm_mask_file = None,
                       image_interp = 'bilinear', output_res=0.5, dst_nodata = -9999):
    
    """
    """

    ortho_cache_dir = CACHE_DIR + "/orthorectify/"
    if os.path.isdir(CACHE_DIR):
        print("ERROR: Please delete CACHE_DIR and rerun")
        sys.exit(1)
        # shutil.rmtree(CACHE_DIR)
    create_directory(ortho_cache_dir)

    # Only orthorecitfy chips for which rpc file is found
    # We require rpc files to be in RPB format
    folder_img = remove_trailing_slash(folder_img)
    check_if_exists(folder_img)
    if folder_rpc:
        folder_rpc = remove_trailing_slash(folder_rpc)
        check_if_exists(folder_rpc)
    check_if_exists(dsm_file)
    check_if_exists(dem_file)
    if dsm_mask_file is not None:
        check_if_exists(dsm_mask_file)
    
    folder_output = remove_trailing_slash(folder_output)    
    create_directory(folder_output)
    
    initial_files_list = glob.glob(folder_img + "/*.tif")
    
    
    rpc_files_list = []
    img_files_list = []
    
    for img_file in initial_files_list:
        just_name = remove_folder_ext(img_file)
        if folder_rpc:
            rpc_file = os.path.join(folder_rpc, just_name + '.RPB')
            if os.path.isfile(rpc_file):
                img_files_list.append(img_file)
                rpc_files_list.append(rpc_file)
            else:
                print("WARNING:RPC file %s does not exist for img %s. Skipping" % (rpc_file, img_file))            
        else:
            img_files_list.append(img_file)
            rpc_files_list.append(None)

            
    # If parallel do orthorectify for img,rpc combo in batches
    input_mp_list = []
    dem_interp = 'near'
    image_interp = image_interp
    output_res = str(output_res)
    dst_nodata = str(dst_nodata)

    for img_file, rpc_file in zip(img_files_list, rpc_files_list):
        input_mp = []
        input_mp.extend([ortho_cache_dir])
        input_mp.extend(['-RPC_DEM', dsm_file])
        input_mp.extend(['-low_res_dem', dem_file])        
        input_mp.extend(['-dem_interp', dem_interp])
        input_mp.extend(['-image_interp', image_interp])
        input_mp.extend(['-output_res', output_res])
        input_mp.extend(['-dst_nodata', dst_nodata])
        input_mp.extend(['-save_mapping_flag'])
        input_mp.extend(['-dem_mask_file', dsm_mask_file])
        output_file = os.path.join(folder_output, os.path.split(img_file)[-1])
        input_mp.extend(['-output', output_file])
        #input_mp.extend(['-utm'])
        #input_mp.extend(['-compress_flag'])
        #input_mp.extend(['-parallel'])
        #input_mp.extend(['-gdal_merge', GDAL_MERGE])
        #input_mp.extend(['-n_workers' N_WORKERS])
        input_mp.extend(['-rpc_file', rpc_file])
        input_mp.extend([img_file])
        input_mp_list.append(input_mp)

    if parallel_flag:
        assert num_workers > 1
        mypool = Pool(num_workers)
        mypool.map(call_ortho, input_mp_list)
        mypool.close()
    else:
        for input_mp in input_mp_list:
            call_ortho(input_mp)            
    
    # remove cache dir
    shutil.rmtree(CACHE_DIR)
     
    return

 
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description = "a module to do orthorectification")

    parser.add_argument('-folder_img', '--folder_img', type = str, help = 'folder_img', required = True)
    parser.add_argument('-folder_rpc', '--folder_rpc', type = str, help = 'folder_rpc', default = None)
    parser.add_argument('-dsm_file', '--dsm_file', type = str, help = 'dsm_file', required = True)
    parser.add_argument('-dem_file', '--dem_file', type = str, help = 'dem_file', required = True)
    parser.add_argument('-folder_output', '--folder_output', type = str, help = 'output folder', required = True)
    parser.add_argument('-cache_dir', '--cache_dir', type = str, help = 'cache_dir', required = True)
    parser.add_argument('-parallel', '--parallel', action="store_true", help = 'parallel')
    parser.add_argument('-num_workers', '--num_workers', default = 1, help = 'num_workers', type = int)
    parser.add_argument('-dsm_mask_file', '--dsm_mask_file', default = None, type = str, help = 'dsm_mask_file')
    parser.add_argument('-image_interp', '--image_interp', default = 'bilinear', type = str, help = 'image_interp')
    parser.add_argument('-output_res', '--output_res', default = 0.5, type = float, help = 'output GSD in metres')
    parser.add_argument('-nodata', '--nodata', default = -9999, type = int, help = 'no data value to indicate occlusion in output')

    args = parser.parse_args()

    batch_orthorectify(args.folder_img, args.folder_rpc, args.dsm_file, args.dem_file,
                       args.folder_output, args.cache_dir, parallel_flag = args.parallel,
                       num_workers = args.num_workers, dsm_mask_file = args.dsm_mask_file,
                       image_interp = args.image_interp, output_res = args.output_res,
                       dst_nodata = args.nodata)
