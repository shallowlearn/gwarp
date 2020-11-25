/*
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
*/

#include "gdal_priv.h" // for common GDAL functions
#include "gdal_alg.h" // for GDALCreateRPCTransformer and GDALRPCTransform function
#include "cpl_conv.h" // for CPLMalloc()
#include <iostream>	// for cout, cerr
#include <string>	// for string datatype
#include <cmath> //NaN
#include <boost/program_options.hpp> // for nice argument parsing
#include <vector>
#include <cstdint> //uint16_t

void convert_pixel_coords_to_projected_coords(double col, double row, std::vector<double> geotransform, std::vector <double>&  proj_x_y)
{

	proj_x_y[0] = geotransform[0] + (geotransform[1]*col) + (geotransform[2]*row);

	proj_x_y[1] = geotransform[3] + (geotransform[4]*col) + (geotransform[5]*row);


}

void convert_projected_coords_to_pixel_coords(double proj_x, double proj_y, std::vector<double>  geotransform, double* col, double* row)
{

	col[0]  = (proj_x - geotransform[0])/geotransform[1];

	row[0]  = (proj_y - geotransform[3])/geotransform[5];

}

template <class T>
bool createOrthoImage(int start_row_final, int end_row_final, int start_col_final, int end_col_final,
					  int ngrid_cols, std::vector<bool>& valid, std::vector<int> image_grid_col, std::vector<int> image_grid_row,
					  int image_pre_ortho_height, int image_pre_ortho_width, std::vector<double> height, std::vector<float> lookup,
					  double height_thresh_filter_metres, int nbands, int outputWidth, int outputHeight,
					  std::string image_interpolate_method, GDALDataset  *poRasterDataset,
					  std::vector<T>& outBandsAsArray, GDALDataType outputDataType,
					  std::vector<float> image_grid_col_float, std::vector<float> image_grid_row_float, float halfshift,
					  std::vector<float>& output_raster_coord_map, bool save_mapping_flag)
{
	//#pragma omp parallel for // Causes segmentation fault for large array
	int x,y, twod_idx, clip_twod_idx, agnostic_idx;
	float w_x, w_y, x_f, y_f;
	T img00, img01, img10, img11;
	int x0, x1, y0, y1;
	int lookup_idx;
	int readflag;
	double imgval;

	for (int j = start_row_final; j < end_row_final; j++)
	{
		for (int i = start_col_final; i < end_col_final; i++)
		{
			twod_idx = (j*ngrid_cols) + i;
			if (valid[twod_idx])
			{
				agnostic_idx = ((j - start_row_final)*outputWidth) + (i - start_col_final);
				if (save_mapping_flag)
				{
					output_raster_coord_map[agnostic_idx] =  image_grid_col_float[twod_idx];
					output_raster_coord_map[agnostic_idx + (outputWidth*outputHeight)] =  image_grid_row_float[twod_idx];
				}

				x = image_grid_col[twod_idx];
				y = image_grid_row[twod_idx];
				lookup_idx = (x*image_pre_ortho_height) + y;
				if ( (lookup[lookup_idx] - height[twod_idx]) > height_thresh_filter_metres )
				{
					valid[twod_idx] = 0;
				}
				else
				{
					if (save_mapping_flag)
						output_raster_coord_map[agnostic_idx + 2*(outputWidth*outputHeight)] =  1;

					if (image_interpolate_method == "near")
					{
						for (int bandIdx = 1; bandIdx <= nbands; bandIdx++)
						{
							clip_twod_idx = ((bandIdx -1)*outputWidth*outputHeight) + ((j - start_row_final)*outputWidth) + (i - start_col_final);
							readflag = poRasterDataset->GetRasterBand(bandIdx)->RasterIO(GF_Read, x, y, 1, 1, &outBandsAsArray[clip_twod_idx], 1, 1, outputDataType, 0, 0);
						}
					}
					else if (image_interpolate_method == "bilinear")
					{

						/*
						 * Turns out the near and iinterpolated image look shifted
						 * if I have an extra halfshift here
						 */
						x_f = image_grid_col_float[twod_idx]; //- halfshift;
						y_f = image_grid_row_float[twod_idx]; //- halfshift;

						x0 = floor(x_f);
						y0 = floor(y_f);

						x1 = x0 + 1;
						y1 = y0 + 1;

						w_x = x_f - (float)(x0);
						w_y = y_f - (float)(y0);

						if (x0 < 0)
							x0 = x1;
						if (y0 < 0)
							y0 = y1;
						if (x1 >= image_pre_ortho_width)
							x1 = x0;
						if (y1 >= image_pre_ortho_height)
							y1 = y0;

						for (int bandIdx = 1; bandIdx <= nbands; bandIdx++)
						{
							clip_twod_idx = ((bandIdx -1)*outputWidth*outputHeight) + ((j - start_row_final)*outputWidth) + (i - start_col_final);
							readflag = poRasterDataset->GetRasterBand(bandIdx)->RasterIO(GF_Read, x0, y0, 1, 1, &img00, 1, 1, outputDataType, 0, 0);
							readflag = poRasterDataset->GetRasterBand(bandIdx)->RasterIO(GF_Read, x1, y0, 1, 1, &img01, 1, 1, outputDataType, 0, 0);
							readflag = poRasterDataset->GetRasterBand(bandIdx)->RasterIO(GF_Read, x0, y1, 1, 1, &img10, 1, 1, outputDataType, 0, 0);
							readflag = poRasterDataset->GetRasterBand(bandIdx)->RasterIO(GF_Read, x1, y1, 1, 1, &img11, 1, 1, outputDataType, 0, 0);

							imgval = ( (((double)(img00))*(1.0 - w_x)*(1.0 - w_y)) + (((double)(img01))*w_x*(1.0 - w_y)) +
									   (((double)(img10))*(1.0 - w_x)*w_y) + (((double)(img11))*w_x*w_y) );

							/*if (outputDataType == GDT_Int16 ||  outputDataType == GDT_UInt16 || outputDataType == GDT_Int32
									|| outputDataType == GDT_Byte || outputDataType == GDT_UInt32)
								outBandsAsArray[clip_twod_idx] = round(imgval);
							else if (outputDataType == GDT_Float32 ||  outputDataType == GDT_Float64)
								outBandsAsArray[clip_twod_idx] = imgval;*/
							outBandsAsArray[clip_twod_idx] = imgval;

						}

					}

				}
			}

		}
	}

	return 1;

}

template <class T>
bool writeGTiffImage(std::string outputRasterName, int outputWidth, int outputHeight, int nbands, GDALDataType outputDataType,
					 std::vector <double> raster_corners_lon, std::vector <double> raster_corners_lat,
					 double *spacing_lon, double *spacing_lat,GDALDataset  *poDemDataset,
					 std::vector<T>& outBandsAsArray,  T dst_nodata, bool *compress_flag,
					 bool tiled_flag = 1)
{
	int writeflag;
	const char *pszFormat = "GTiff";
	GDALDriver *poDriver;
	char **papszMetadata;
	poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
	if( poDriver == NULL )
		exit( 1 );
	papszMetadata = poDriver->GetMetadata();
	/*if( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE ) )
		printf( "Driver %s supports Create() method.\n", pszFormat );
	else
		return 0;*/

	GDALDataset *poDstDS;
	char **ds_papszOptions = NULL;
	if (tiled_flag)
		ds_papszOptions = CSLSetNameValue( ds_papszOptions, "TILED", "YES" );
	if (compress_flag[0])
	{	// In case of compression, gdal wont know if BIGTIFF is needed or not. So force for safety
		ds_papszOptions = CSLSetNameValue( ds_papszOptions, "COMPRESS", "LZW" );
		ds_papszOptions = CSLSetNameValue( ds_papszOptions, "BIGTIFF", "YES" );
	}

	poDstDS = poDriver->Create(outputRasterName.c_str(), outputWidth, outputHeight, nbands, outputDataType, ds_papszOptions);
	double geotransform[6] = {raster_corners_lon[0] - (spacing_lon[0]/2.0), spacing_lon[0], 0, raster_corners_lat[0] - (spacing_lat[0]/2.0), 0, spacing_lat[0]};
	poDstDS->SetGeoTransform( geotransform );

	poDstDS->SetProjection( poDemDataset->GetProjectionRef());
	GDALRasterBand *poDstBand;

	//std::vector<T> outBand (outputWidth*outputHeight, dst_nodata);
	for (int bandIdx = 1; bandIdx <= nbands; bandIdx++)
	{
		poDstBand = poDstDS->GetRasterBand(bandIdx);
		poDstBand->SetNoDataValue(dst_nodata);
		//std::copy (outBandsAsArray.begin() + (bandIdx -1)*outputWidth*outputHeight, outBandsAsArray.begin() + (bandIdx)*outputWidth*outputHeight, outBand.begin());
		//writeflag = poDstBand->RasterIO( GF_Write, 0, 0, outputWidth, outputHeight, &outBand[0], outputWidth, outputHeight, outputDataType, 0, 0);
		writeflag = poDstBand->RasterIO( GF_Write, 0, 0, outputWidth, outputHeight, &outBandsAsArray[(bandIdx -1)*outputWidth*outputHeight], outputWidth, outputHeight, outputDataType, 0, 0);
	}
	// Once we're done, close the dataset

	GDALClose( (GDALDatasetH) poDstDS );

	return 1;
}


std::string replaceFirstOccurrence(std::string& s, const std::string& toReplace, const std::string& replaceWith)
{
    std::size_t pos = s.find(toReplace);
    if (pos == std::string::npos) return s;
    return s.replace(pos, toReplace.length(), replaceWith);
}

bool get_ortho_grid_worker (std::string outputRasterName, double *start_lon, double *end_lon,
						   double *start_lat, double *end_lat, double *spacing_lon, double *spacing_lat,
						   std::string raster_file, std::string dem_file, std::string dem_interpolate_method,
						   std::string image_interpolate_method, double *dst_nodata, int*block_idx, bool *save_mapping_flag,
						   bool *border_block_flag, std::string rpc_file, bool *compress_flag, std::string dem_mask_file,
						   std::string dtm_file)

{

	/* function that implements gwarp++ for block. A block could also cover an entire raster if it is small enough.
	By block we refer to a grid in lon lat coordinates */

	float halfshift = 0.5;
	/* print inputs
	std::cout << std::endl;
	std::cout << outputRasterName << std::endl;
	std::cout << *start_lon << std::endl;
	std::cout << *end_lon << std::endl;
	std::cout << *start_lat << std::endl;
	std::cout << *end_lat << std::endl;
	std::cout << *spacing_lon << std::endl;
	std::cout << *spacing_lat << std::endl;
	std::cout << raster_file << std::endl;
	std::cout << dem_file << std::endl;
	std::cout << dem_interpolate_method << std::endl;
	std::cout << image_interpolate_method << std::endl;
	std::cout << *dst_nodata << std::endl;
	std::cout << *block_idx << std::endl;
	std::cout << *save_mapping_flag << std::endl;
	std::cout << *border_block_flag << std::endl << std::endl; */

	/* We have to account for the effects of occlusion. Since we split the lon lat grid into blocks, we need to account for tall objects
   nearby. So we add a padding of 100 points in each direction.
   # BTODO - Make padding dependent on data. */
	double padding = 200;

	/* A typical storey height is 3 metres. So all points that get projected to the same pixel whose height is within 1.5 metres
	of the maximum height for that pixel are considered valid. This is to reduce the holes in actual building rooftops */
	double height_thresh_filter_metres = 1.5;

	bool valid_flag = 0;

	float no_data_val = -32000;
	// Encapsulate within a try catch statement to detect blocks that fail
	try
	{
		// register all known drivers
		GDALAllRegister();

		// pointers to the raster and dem dataset
		GDALDataset  *poRasterDataset, *poDemDataset, *poMaskDemDataset, *poDtmDataset;

		// open raster dataset
		poRasterDataset = (GDALDataset *) GDALOpen(raster_file.c_str(), GA_ReadOnly);

		if( poRasterDataset == NULL )
		{
			std::cout << "ERROR IN OPENING RASTER DATASET" << std::endl;
			abort();
		}

		/*std::cout << poRasterDataset->GetRasterBand(1)->GetNoDataValue() << std::endl;

		if (std::isnan(no_data_input_value))
		{
			std::cout << "JJD";
			no_data_input_value = no_data_val;
		}
		 */
		// open dem dataset
		poDemDataset = (GDALDataset *) GDALOpen(dem_file.c_str(), GA_ReadOnly);

		if( poDemDataset == NULL )
		{
			std::cout << "ERROR IN OPENING DEM DATASET" << std::endl;
			abort();
		}

		poMaskDemDataset = NULL;
		if (dem_mask_file != "None")
		{
			poMaskDemDataset = (GDALDataset *) GDALOpen(dem_mask_file.c_str(), GA_ReadOnly);

			if( poMaskDemDataset == NULL )
			{
				std::cout << "ERROR IN OPENING MASK FOR DEM DATASET" << std::endl;
				abort();
			}
		}

		poDtmDataset = NULL;
		if (dtm_file != "None")
		{
			poDtmDataset = (GDALDataset *) GDALOpen(dtm_file.c_str(), GA_ReadOnly);

			if( poDtmDataset == NULL )
			{
				std::cout << "ERROR IN OPENING DTM DATASET" << std::endl;
				abort();
			}
		}


		// Height and width of dem
		int dem_width = poDemDataset->GetRasterXSize();
		int dem_height = poDemDataset->GetRasterYSize();

		// The dem can have no data value as negative integers if it is stored as uint16.
		double no_data_height_value = (double)poDemDataset->GetRasterBand(1)->GetNoDataValue();

		/* Declaring it as GDALRPCInfo rpcinfo could blow out the stack memory */
		GDALRPCInfo *rpcinfo = new GDALRPCInfo;
		char *papszOptions = NULL;
		void *rpc_transformer;

		if( !GDALExtractRPCInfo(poRasterDataset->GetMetadata( "RPC" ), rpcinfo) )
		{
			std::cout << " ERROR: No rpc metadata found" << std::endl;
			abort();
		}

		// Create roc transformer. Note gdal's convention is opposite to normal rpc. So by default
		// it maps pixel line sample to lon, lat
		/*GDALCreateRPCTransformer( GDALRPCInfo *psRPCInfo, int bReversed, double dfPixErrThreshold, char **papszOptions )*/
		rpc_transformer = GDALCreateRPCTransformer(rpcinfo, 0, 0, &papszOptions);

		// Raster corners for this block
		std::vector <double> raster_corners_lon(2);
		raster_corners_lon[0] = start_lon[0];
		raster_corners_lon[1] = end_lon[0];

		std::vector <double> raster_corners_lat(2);
		raster_corners_lat[0] = start_lat[0];
		raster_corners_lat[1] = end_lat[0];

		std::vector<double> demGeoTransform(6);
		poDemDataset->GetGeoTransform(&demGeoTransform[0]);


		std::vector <double>  dem_top_left_lon_lat(2);
		convert_pixel_coords_to_projected_coords(0.0, 0.0, demGeoTransform, dem_top_left_lon_lat);

		std::vector <double> dem_bottom_right_lon_lat(2);
		convert_pixel_coords_to_projected_coords((double) dem_width, (double)dem_height, demGeoTransform, dem_bottom_right_lon_lat);

        /* Skip if dsm does not cover block. Otherwise the rpc transformer assumes
         * 0 for NaNs and creates ghosts at the boundaries of the dem and the image.*/

        if (start_lon[0] < dem_top_left_lon_lat[0] || start_lat[0] > dem_top_left_lon_lat[1])
		{
            printf("\nDEM does not cover block (start_lon,start_lat=%lf,%lf dem_corners=%lf,%lf)\n", start_lon[0], start_lat[0], dem_top_left_lon_lat[0], dem_top_left_lon_lat[1]);
            //return 0;
		}

        if (end_lon[0] > dem_bottom_right_lon_lat[0] || end_lat[0] < dem_bottom_right_lon_lat[1])
		{
            printf("\nDEM does not cover block (end_lon,end_lat=%lf,%lf dem_corners=%lf,%lf)\n", end_lon[0], end_lat[0], dem_bottom_right_lon_lat[0], dem_bottom_right_lon_lat[1]);
            //return 0;
		}

        //printf("%.15f,%.15f\n", dem_top_left_lon_lat[0], dem_top_left_lon_lat[1]);

        //printf("%.15f,%.15f\n", dem_bottom_right_lon_lat[0], dem_bottom_right_lon_lat[1]);

        // # Create the grid of points in lon lat coordinates with the padding
        int ngrid_cols = round( ((raster_corners_lon[1] + (padding*spacing_lon[0]) - (raster_corners_lon[0] - (padding*spacing_lon[0])))/spacing_lon[0]) );
        int ngrid_rows = round( ((raster_corners_lat[1] + (padding*spacing_lat[0]) - (raster_corners_lat[0] - (padding*spacing_lat[0])))/spacing_lat[0]) );
        //std::cout << ngrid_cols << " " << ngrid_rows << " " << raster_corners_lat[1]<< " " << raster_corners_lat[0] << " " << spacing_lat[0] << std::endl;

		// Get the indices of the actual grid we save, i.e  minus the padding
        int start_col_final = padding;
        int end_col_final = (start_col_final + (1.0*(raster_corners_lon[1] - raster_corners_lon[0])/spacing_lon[0]) + 1);
		/*if (border_block_flag[0])
			//end_col_final = end_col_final + 1;*/
        int start_row_final = padding;
        int end_row_final = (start_row_final + (1.0*(raster_corners_lat[1] - raster_corners_lat[0])/spacing_lat[0]) + 1);
		/*if (border_block_flag[0])
			end_row_final = end_row_final + 1;*/

		//std::cout << end_col_final << " " << start_col_final << " "  << end_row_final << " " << start_row_final << std::endl << std::endl;
		//std::cout << end_col_final - start_col_final << " "  << end_row_final - start_row_final << std::endl << std::endl;

        //double* ortho_grid_lon = new double[ngrid_cols*ngrid_rows];
        //double* ortho_grid_lat = new double[ngrid_cols*ngrid_rows];

        std::vector <float> image_grid_col_float (ngrid_cols*ngrid_rows, no_data_val);
        std::vector <float> image_grid_row_float (ngrid_cols*ngrid_rows, no_data_val);

        std::vector <int> image_grid_col (ngrid_cols*ngrid_rows, no_data_val);
		std::vector <int> image_grid_row (ngrid_cols*ngrid_rows, no_data_val);

        std::vector <double> height (ngrid_cols*ngrid_rows, no_data_height_value);
        std::vector <bool> valid (ngrid_cols*ngrid_rows, 0);

        double dem_col;
        double dem_row;
        double height_val;
        double lon;
        double lat;


        GDALRasterBand  *poDemBand;
        poDemBand = poDemDataset->GetRasterBand( 1 );

        GDALRasterBand *poMaskDemBand;
        if ( dem_mask_file != "None" )
        {
        	poMaskDemBand = poMaskDemDataset->GetRasterBand( 1 );
        }

        GDALRasterBand *poDtmBand;
		std::vector<double> dtmGeoTransform(6,0);
		double no_data_dtm_value;
        if ( dtm_file != "None" )
        {
        	poDtmBand = poDtmDataset->GetRasterBand( 1 );
    		poDtmDataset->GetGeoTransform(&dtmGeoTransform[0]);
    		no_data_dtm_value = (double)poDtmBand->GetNoDataValue();
        }

        int image_pre_ortho_width = poRasterDataset->GetRasterXSize();
        int image_pre_ortho_height = poRasterDataset->GetRasterYSize();
        int nbands = poRasterDataset->GetRasterCount();
        int lookup_idx;
        int readflag, writeflag, readmaskflag,readdtmflag;

        std::vector<float> lookup (image_pre_ortho_width*image_pre_ortho_height, no_data_height_value);

        //std::cout << no_data_height_value << std::endl;
        //std::cout << lookup[1] << std::endl;
        int twod_idx;
        int wall_twod_idx;

        int num_wallpts = 10; //21;
        double wall_step = 0; //0.5;
        double wall_step_size = 0.5;

        uint8_t dem_mask_value = 0;

        int16_t dtm_height_val;
        double dtm_col;
        double dtm_row;

		// Height and width of dtm
		int dtm_width = poDtmDataset->GetRasterXSize();
		int dtm_height = poDtmDataset->GetRasterYSize();

		//#pragma omp parallel for // Causes segmentation fault for large array
		for (int j = 0; j < ngrid_rows; j++)
		{
			for (int i=0; i < ngrid_cols; i++)
			{
				twod_idx = (j*ngrid_cols) + i;
				lon = raster_corners_lon[0] + (((double)(i - padding))*spacing_lon[0]);

				lat = raster_corners_lat[0] + (((double)(j - padding))*spacing_lat[0]);

				convert_projected_coords_to_pixel_coords(lon, lat, demGeoTransform, &dem_col, &dem_row);

				if ( dem_col < 0 ||  dem_col >= dem_width || dem_row < 0 ||  dem_row >= dem_height)
				{
					continue;
				}
				else
				{
					readflag = poDemBand->RasterIO(GF_Read, floor(dem_col), floor(dem_row), 1, 1, &height[twod_idx], 1, 1, GDT_Float64, 0, 0);

					if ( dem_mask_file != "None" )
					{
						readmaskflag = poMaskDemBand->RasterIO(GF_Read, floor(dem_col), floor(dem_row), 1, 1, &dem_mask_value, 1, 1, GDT_Byte, 0, 0);
						if (dem_mask_value == 255)
						{
							std::cout << "Skipping Water Point" << std::endl;
							dem_mask_value = 0; // Reset dem_mask_value
							continue;
						}
					}

					if (height[twod_idx] != no_data_height_value && ! std::isnan(height[twod_idx]))
					{

						if ( dtm_file != "None ")
						{
							// Load height value from dtm and calculate num_wallpts accordingly
							convert_projected_coords_to_pixel_coords(lon, lat, dtmGeoTransform, &dtm_col, &dtm_row);

							if ( dtm_col < 0 ||  dtm_col >= dtm_width || dtm_row < 0 ||  dtm_row >= dtm_height) // check if we are inside dtm
							{
								std::cout << "WARNING OUTSIDE DTM FOR LON: " << lon << " , LAT: " << lat << std::endl;
							}
							else
							{
								readdtmflag = poDtmBand->RasterIO(GF_Read, floor(dtm_col), floor(dtm_row), 1, 1, &dtm_height_val, 1, 1, GDT_Int16, 0, 0);
								if ((double) dtm_height_val == no_data_dtm_value || std::isnan(dtm_height_val)) // check if we have valid dtm value
								{
									std::cout << "WARNING NO DTM FOUND FOR LON: " << lon << " , LAT: " << lat << std::endl;
									std::cout << "WARNING NO DTM FOUND FOR LON: " << dtm_col << " , LAT: " << dtm_row << std::endl;
									std::cout << "WARNING NO DTM FOUND : " << dtm_height_val <<  std::endl;
								}
								else if (1.0*(double)(dtm_height_val) < height[twod_idx]) // check if dem > dtm. If so update num_wallpts
								{
									num_wallpts = ceil((height[twod_idx] - (1.0*(double)(dtm_height_val)))/wall_step_size);
									if (num_wallpts <= 0)
									{
										num_wallpts = 1;
									}
								}
							}
						}
						//std::cout << num_wallpts << std::endl;
				        std::vector<double> height_walls(num_wallpts, 0);
				        std::vector<double> lon_walls(num_wallpts, 0);
				        std::vector<double> lat_walls(num_wallpts, 0);
				        std::vector<int> success_flag(num_wallpts, 0);

						height_walls[0] = height[twod_idx];
						lon_walls[0] = lon;
						lat_walls[0] = lat;
						success_flag[0] = 0;

						for (int wall_slice_idx = 1; wall_slice_idx < num_wallpts; wall_slice_idx++)
						{
							height_walls[wall_slice_idx] =  height[twod_idx] - wall_step_size*(((double)(wall_slice_idx)) + wall_step); // 1.0*(wall_slice_idx*wall_step);//- 1.0*(wall_slice_idx - 1 + wall_step);
							lon_walls[wall_slice_idx] = lon;
							lat_walls[wall_slice_idx] = lat;
							success_flag[wall_slice_idx] = 0;
						}

						GDALRPCTransform (rpc_transformer, 1, num_wallpts, &lon_walls[0], &lat_walls[0], &height_walls[0], &success_flag[0]);

						for (int wall_slice_idx = 0; wall_slice_idx < num_wallpts; wall_slice_idx++)
						{
							lon = lon_walls[wall_slice_idx];
							lat = lat_walls[wall_slice_idx];
							height_val = height_walls[wall_slice_idx];
							if (success_flag[wall_slice_idx] == 1 && lon >= 0 && lon < image_pre_ortho_width && lat >= 0 && lat < image_pre_ortho_height)
							{
								if (wall_slice_idx == 0)
								{
									valid[twod_idx] = 1;
									image_grid_col_float[twod_idx] = lon - halfshift;// I believe this is causing a bug as 0.49999 is getting cast to 0.5
									image_grid_row_float[twod_idx] = lat - halfshift;

									image_grid_col[twod_idx] = floor(lon);
									image_grid_row[twod_idx] = floor(lat);

									valid_flag = 1;
								}

								lookup_idx = (floor(lon)*image_pre_ortho_height) + floor(lat);

								if ( lookup[lookup_idx] == no_data_height_value  || std::isnan(lookup[lookup_idx]) )
								{
									lookup[lookup_idx] = height_val;
								}
								else
								{
									if ((height_val - lookup[lookup_idx]) > 0)//height_thresh_filter_metres)
									{
										lookup[lookup_idx] = height_val;
									}
								}
							}
							//std::cout << lon << " " << lat << " " <<  height[200] << std::endl;
						}

					}
				}
			}

		}

		if (!valid_flag)
		{
		  printf("\nNo valid points for block #%d\n", block_idx[0]);
		  printf("\nFinished processing block #%d\n", block_idx[0]);
		  return 0;
		}
		int outputWidth = end_col_final - start_col_final;
		int outputHeight = end_row_final - start_row_final;

		GDALDataType outputDataType = poRasterDataset->GetRasterBand(1)->GetRasterDataType();

		//printf("%s",GDALGetDataTypeName(outputDataType));
		int dst_nodata_int = (int) dst_nodata[0];

		if (outputDataType == GDT_UInt16 && dst_nodata_int < 0)
			outputDataType = GDT_Int16;

		if (outputDataType == GDT_Byte)
		{
			if (dst_nodata_int < 0 || dst_nodata_int > 255)
				outputDataType = GDT_Int16;
		}

		std::vector<float> output_raster_coord_map(3*outputWidth*outputHeight, (float)dst_nodata[0]);
		// free memory allocated to coord_map if not needed
		if (! save_mapping_flag[0])
			std::vector<float>().swap(output_raster_coord_map);
		bool tiled_flag = 0;

		if (outputDataType == GDT_Int16)
		{

			std::vector<int16_t> outBandsAsArray(nbands*outputWidth*outputHeight, dst_nodata_int);
			//printf("\nCreating output ortho grid\n");
			bool flag = createOrthoImage(start_row_final, end_row_final, start_col_final, end_col_final,
										 ngrid_cols, valid, image_grid_col, image_grid_row,
										 image_pre_ortho_height, image_pre_ortho_width, height, lookup,
										 height_thresh_filter_metres, nbands, outputWidth, outputHeight, image_interpolate_method,
										 poRasterDataset, outBandsAsArray, outputDataType,
										 image_grid_col_float, image_grid_row_float, halfshift, output_raster_coord_map,
										 save_mapping_flag[0]);

			//printf("\nSaving output image\n");
			bool writeSuccess = writeGTiffImage(outputRasterName, outputWidth, outputHeight, nbands, outputDataType,
												raster_corners_lon, raster_corners_lat,
												spacing_lon, spacing_lat, poDemDataset,
												outBandsAsArray, (int16_t) dst_nodata_int, compress_flag);

			if (save_mapping_flag[0])
			{

				//printf("\nSaving output coordmap\n");
				std::string outputCoordMapRasterName = replaceFirstOccurrence(outputRasterName, ".tif", "_coord_map.tif");
				bool writeSuccess = writeGTiffImage(outputCoordMapRasterName, outputWidth, outputHeight, 3, GDT_Float32,
													raster_corners_lon, raster_corners_lat,
													spacing_lon, spacing_lat, poDemDataset,
													output_raster_coord_map, (float) dst_nodata[0],
													compress_flag, tiled_flag);
			}

		}
		else if (outputDataType == GDT_UInt16)
		{

			std::vector<uint16_t> outBandsAsArray(nbands*outputWidth*outputHeight, dst_nodata_int);
			//printf("\nCreating output ortho grid\n");
			bool flag = createOrthoImage(start_row_final, end_row_final, start_col_final, end_col_final,
										 ngrid_cols, valid, image_grid_col, image_grid_row,
										 image_pre_ortho_height, image_pre_ortho_width, height, lookup,
										 height_thresh_filter_metres, nbands, outputWidth, outputHeight, image_interpolate_method,
										 poRasterDataset, outBandsAsArray, outputDataType,
										 image_grid_col_float, image_grid_row_float, halfshift, output_raster_coord_map,
										 save_mapping_flag[0]);

			//printf("\nSaving output image\n");
			bool writeSuccess = writeGTiffImage(outputRasterName, outputWidth, outputHeight, nbands, outputDataType,
												raster_corners_lon, raster_corners_lat,
												spacing_lon, spacing_lat, poDemDataset,
												outBandsAsArray, (uint16_t) dst_nodata_int, compress_flag);

			if (save_mapping_flag[0])
			{

				//printf("\nSaving output coordmap\n");
				std::string outputCoordMapRasterName = replaceFirstOccurrence(outputRasterName, ".tif", "_coord_map.tif");
				bool writeSuccess = writeGTiffImage(outputCoordMapRasterName, outputWidth, outputHeight, 3, GDT_Float32,
													raster_corners_lon, raster_corners_lat,
													spacing_lon, spacing_lat, poDemDataset,
													output_raster_coord_map, (float) dst_nodata[0],
													compress_flag, tiled_flag);
			}

		}
		else if (outputDataType == GDT_Byte)
		{

			std::vector<uint8_t> outBandsAsArray(nbands*outputWidth*outputHeight, dst_nodata_int);
			//printf("\nCreating output ortho grid\n");
			bool flag = createOrthoImage(start_row_final, end_row_final, start_col_final, end_col_final,
										 ngrid_cols, valid, image_grid_col, image_grid_row,
										 image_pre_ortho_height, image_pre_ortho_width, height, lookup,
										 height_thresh_filter_metres, nbands, outputWidth, outputHeight, image_interpolate_method,
										 poRasterDataset, outBandsAsArray, outputDataType,
										 image_grid_col_float, image_grid_row_float, halfshift, output_raster_coord_map,
										 save_mapping_flag[0]);

			//printf("\nSaving output image\n");
			bool writeSuccess = writeGTiffImage(outputRasterName, outputWidth, outputHeight, nbands, outputDataType,
												raster_corners_lon, raster_corners_lat,
												spacing_lon, spacing_lat, poDemDataset,
												outBandsAsArray, (uint8_t) dst_nodata_int, compress_flag);

			if (save_mapping_flag[0])
			{

				//printf("\nSaving output coordmap\n");
				std::string outputCoordMapRasterName = replaceFirstOccurrence(outputRasterName, ".tif", "_coord_map.tif");
				bool writeSuccess = writeGTiffImage(outputCoordMapRasterName, outputWidth, outputHeight, 3, GDT_Float32,
													raster_corners_lon, raster_corners_lat,
													spacing_lon, spacing_lat, poDemDataset,
													output_raster_coord_map, (float) dst_nodata[0],
													compress_flag, tiled_flag);
			}

		}
		else if (outputDataType == GDT_Float32)
		{

			std::vector<float> outBandsAsArray(nbands*outputWidth*outputHeight, dst_nodata_int);
			//printf("\nCreating output ortho grid\n");
			bool flag = createOrthoImage(start_row_final, end_row_final, start_col_final, end_col_final,
										 ngrid_cols, valid, image_grid_col, image_grid_row,
										 image_pre_ortho_height, image_pre_ortho_width, height, lookup,
										 height_thresh_filter_metres, nbands, outputWidth, outputHeight, image_interpolate_method,
										 poRasterDataset, outBandsAsArray, outputDataType,
										 image_grid_col_float, image_grid_row_float, halfshift, output_raster_coord_map,
										 save_mapping_flag[0]);

			//printf("\nSaving output image\n");
			bool writeSuccess = writeGTiffImage(outputRasterName, outputWidth, outputHeight, nbands, outputDataType,
												raster_corners_lon, raster_corners_lat,
												spacing_lon, spacing_lat, poDemDataset,
												outBandsAsArray, (float) dst_nodata[0], compress_flag);

			if (save_mapping_flag[0])
			{

				//printf("\nSaving output coordmap\n");
				std::string outputCoordMapRasterName = replaceFirstOccurrence(outputRasterName, ".tif", "_coord_map.tif");
				bool writeSuccess = writeGTiffImage(outputCoordMapRasterName, outputWidth, outputHeight, 3, GDT_Float32,
													raster_corners_lon, raster_corners_lat,
													spacing_lon, spacing_lat, poDemDataset,
													output_raster_coord_map, (float) dst_nodata[0],
													compress_flag, tiled_flag);
			}

		}
		else if (outputDataType == GDT_Float64)
		{

			std::vector<double> outBandsAsArray(nbands*outputWidth*outputHeight, dst_nodata_int);
			//printf("\nCreating output ortho grid\n");
			bool flag = createOrthoImage(start_row_final, end_row_final, start_col_final, end_col_final,
										 ngrid_cols, valid, image_grid_col, image_grid_row,
										 image_pre_ortho_height, image_pre_ortho_width, height, lookup,
										 height_thresh_filter_metres, nbands, outputWidth, outputHeight, image_interpolate_method,
										 poRasterDataset, outBandsAsArray, outputDataType,
										 image_grid_col_float, image_grid_row_float, halfshift, output_raster_coord_map,
										 save_mapping_flag[0]);

			//printf("\nSaving output image\n");
			bool writeSuccess = writeGTiffImage(outputRasterName, outputWidth, outputHeight, nbands, outputDataType,
												raster_corners_lon, raster_corners_lat,
												spacing_lon, spacing_lat, poDemDataset,
												outBandsAsArray, (double) dst_nodata[0], compress_flag);

			if (save_mapping_flag[0])
			{

				//printf("\nSaving output coordmap\n");
				std::string outputCoordMapRasterName = replaceFirstOccurrence(outputRasterName, ".tif", "_coord_map.tif");
				bool writeSuccess = writeGTiffImage(outputCoordMapRasterName, outputWidth, outputHeight, 3, GDT_Float32,
													raster_corners_lon, raster_corners_lat,
													spacing_lon, spacing_lat, poDemDataset,
													output_raster_coord_map, (float) dst_nodata[0],
													compress_flag, tiled_flag);
			}

		}

		// close datasets
		GDALClose(poRasterDataset);
		GDALClose(poDemDataset);

		if (dem_mask_file != "None")
		{
			GDALClose(poMaskDemDataset);
		}

		if (dtm_file != "None")
		{
			GDALClose(poDtmDataset);
		}

		// destroy rpc transformer
		GDALDestroyRPCTransformer(rpc_transformer);

		// destroy variables on heap
		delete rpcinfo;
		printf("\nFinished processing block #%d\n", block_idx[0]);
		return 1;
	}
	catch(std::exception &err)
	{
		std::cerr << "Unhandled Exception" << err.what() << ", Exiting " << std::endl;
		return 0;
	}

}


int main(int argc, char **argv)
{
	// Set precision of std::cout
	std::cout.precision(17);

	// parse input args
	std::string outputRasterName;
	double start_lon;
	double end_lon;
	double start_lat;
	double end_lat;
	double spacing_lon;
	double spacing_lat;
	std::string raster_file;
	std::string dem_file;
	std::string dem_interpolate_method;
	std::string image_interpolate_method;
	std::string rpc_file;
	double dst_nodata;;
	int block_idx;
	bool save_mapping_flag;
	bool border_block_flag;
	bool compress_flag;
	std::string dem_mask_file;
	std::string dtm_file;


	try
	{
		namespace po = boost::program_options;
		po::options_description desc("Allowed options");
		desc.add_options()
				("help", "Help for gwarp++.cpp")
				("output", po::value<std::string>(&outputRasterName)->required(), "Output raster file name")
				("start_lon", po::value<double>(&start_lon)->required(), "Start value for longitude")
				("end_lon", po::value<double>(&end_lon)->required(), "End value for longitude")
				("start_lat", po::value<double>(&start_lat)->required(), "Starting value for latitude")
				("end_lat", po::value<double>(&end_lat)->required(), "End value for latitude")
				("spacing_lon", po::value<double>(&spacing_lon)->required(), "Spacing for grid in longitude axis")
				("spacing_lat", po::value<double>(&spacing_lat)->required(), "Spacing for grid in latitude axis")
				("raster_file", po::value<std::string>(&raster_file)->required(),"Input raster file name")
				("dem_file", po::value<std::string>(&dem_file)->required(),"Input dem file name")
				("dem_interpolate_method", po::value<std::string>(&dem_interpolate_method)->required(),"dem interpolation method, One of near, bilinear or cubic is supported.")
				("image_interpolate_method", po::value<std::string>(&image_interpolate_method)->required(),"image interpolation method. Near or Bilinear is supported")
				("dst_nodata", po::value<double>(&dst_nodata)->required(), "No data value in output")
				("block_idx", po::value<int>(&block_idx)->required(), "Index of current block")
				("save_mapping_flag", po::value<bool>(&save_mapping_flag)->required(), "Flag to save mapping")
				("border_block_flag", po::value<bool>(&border_block_flag)->required(), "Flag to indicate whether current block is a border block or not")
				("rpc_file", po::value<std::string>(&rpc_file)->required(),"Input rpc file name")
				("compress_flag", po::value<bool>(&compress_flag)->required(), "Flag to turn on compression")
				("dem_mask_file", po::value<std::string>(&dem_mask_file)->default_value("None"),"Optional mask file for dsm")
				("dtm_file", po::value<std::string>(&dtm_file)->default_value("None"),"Optional input dtm_file file name");

		po::variables_map vm;

		try
		{
			po::store(po::command_line_parser(argc, argv)
					.options(desc)
					.style(
							po::command_line_style::unix_style |
							po::command_line_style::allow_long_disguise
						   )
					.run(), vm);

			po::notify(vm);

		}
		catch(po::error& err)
		{
			std::cerr << "ERROR: " << err.what() << std::endl << std::endl;
			std::cerr << desc << std::endl;
			return 0;
		}

	}
	catch (std::exception& err)
	{
		std::cerr << "Unhandled Exception" << err.what() << ", Exiting " << std::endl;
		return 0;

	}

	bool flag = get_ortho_grid_worker (outputRasterName, &start_lon, &end_lon,
								  	   &start_lat, &end_lat, &spacing_lon, &spacing_lat,
									   raster_file, dem_file, dem_interpolate_method,
									   image_interpolate_method, &dst_nodata, &block_idx, &save_mapping_flag,
									   &border_block_flag, rpc_file, &compress_flag, dem_mask_file, dtm_file);


}

