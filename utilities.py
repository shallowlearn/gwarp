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

import sys,os

                                                                                                                                                      
def remove_folder_from_name(input_name):                                                                                                              
    return os.path.split(os.path.abspath(input_name))[-1]                                                                                             
                                                                                                                                                      
def get_folder(input_name):                                                                                                                           
    return os.path.split(os.path.abspath(input_name))[0]                                                                                              
                                                                                                                                                      
def get_ext(input_name):                                                                                                                              
    return os.path.splitext(input_name)[-1]                                                                                                           
                                                                                                                                                      
def remove_folder_ext(input_name):                                                                                                                    
    return os.path.splitext(os.path.split(os.path.abspath(input_name))[-1])[0]                                                                        
                                                                                    
def create_directory(folder_name):                                                                                                                    
                                                                                                                                                      
    if not os.path.isdir(folder_name):                                                                                                                
        # Sometimes we get race conditions. To stop script from throwing error and exiting                                                            
        # I am actually surprised that I did not hit this till now                                                                                    
        try:                                                                                                                                          
            os.makedirs(folder_name)                                                                                                                  
        except OSError as e:                                                                                                                              
            if e.errno != errno.EEXIST:                                                                                                               
                print("\nERROR: Could not create %s" % folder_name)                                                                                   
                raise                                                                                                                                 
                                                                                                                                                      
def remove_trailing_slash(input_folder):                                                                                                              
    if input_folder.endswith('/'):                                                                                                                    
        return input_folder[:-1]                                                                                                                      
    else:                                                                                                                                             
        return input_folder                                                                                                                           
                                                                                                                                                      
def check_if_exists(input_file_or_folder):                                                                                                            
    if os.path.isfile(input_file_or_folder) or os.path.isdir(input_file_or_folder):                                                                   
        return True                                                                                                                                   
    else:                                                                                                                                             
        raise ValueError("\n \n ERROR: " + input_file_or_folder + " does not exist. Stopping \n")                                                     
        return False