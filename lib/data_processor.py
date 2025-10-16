"""
Data processing module for spiral sun observations.
Handles FITS file processing, coordinate transformations, and data preparation.
"""

import os
import pandas as pd
from lib.utils import bintable_to_pandas, processData, getFinalProcessedData, process_all_heliocentric_coordinates


class SolarDataProcessor:
    """Processes FITS files and prepares data for visualization."""
    
    def __init__(self, observation_config):
        """
        Initialize the data processor.
        
        Args:
            observation_config: Configuration object containing observation parameters
        """
        self.config = observation_config
        
    def process_fits_file(self, fit_file_path, sun_positions):
        """
        Process FITS file and return processed data for all bands.
        
        Args:
            fit_file_path: Path to the FITS file
            sun_positions: DataFrame with UTC, SunX, SunY columns
            
        Returns:
            dict: Processed data for each frequency band
        """
        print('-------------------------------------------------------------')
        print('Start processing the FITS file:', fit_file_path)
        
        # Convert binary table to pandas DataFrame
        hdu_number = 1  # Number of the extension containing the binary table
        data_df = bintable_to_pandas(fit_file_path, hdu_number)        
  
        print('Processing the data...')
        band_data_dfs = processData(data_df)

        print('Getting final processed data...')
        processed_dfs = getFinalProcessedData(self.config, sun_positions, band_data_dfs)

        print('Processing heliocentric coordinates...')
        band_processed_helio_dfs = process_all_heliocentric_coordinates(processed_dfs, self.config)
        
        return band_processed_helio_dfs
    
    def create_sun_position_dataframe(self, utc, xx1, yy1):
        """
        Create DataFrame with sun position data.
        
        Args:
            utc: Array of UTC times
            xx1, yy1: Sun position coordinates
            
        Returns:
            pd.DataFrame: DataFrame with UTC, SunX, SunY columns
        """
        return pd.DataFrame({'UTC': utc, 'SunX': xx1, 'SunY': yy1})
    
    def prepare_output_directory(self, images_path):
        """
        Create output directory for images.
        
        Args:
            images_path: Base path for images
            
        Returns:
            str: Full path to the created directory
        """
        directory = f'{self.config.year:04d}-{self.config.month:02d}-{self.config.day:02d}T{self.config.hour_start:02d}_{self.config.minute_start:02d}_00'
        full_path = f'{images_path}/{directory}'
        
        # Create the directory if it doesn't exist
        if not os.path.exists(full_path):
            os.makedirs(full_path)
            
        return full_path