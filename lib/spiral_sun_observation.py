"""
Refactored SpiralSunObservation class.
Main orchestrator for spiral solar observations using modular components.
"""

from datetime import datetime
from lib.observation_config import ObservationConfig
from lib.position_generator import SpiralPositionGenerator
from lib.file_generator import PTFFileGenerator
from lib.data_processor import SolarDataProcessor
from lib.solar_plotter import SolarMapPlotter
from lib.antenna import Antenna
from lib.weather import Weather


class SpiralSunObservation:
    """
    Main class for spiral solar observations.
    
    This class orchestrates the entire observation process by coordinating
    specialized modules for position generation, file creation, data processing,
    and visualization.
    """
    
    def __init__(self, weather: Weather, antenna_instance: Antenna, observation_datetime: datetime):
        """
        Initialize the spiral sun observation.
        
        Args:
            weather: Weather instance with atmospheric conditions
            antenna_instance: Antenna instance with location and properties  
            observation_datetime: Start time of the observation
        """
        # Create configuration object that holds all parameters
        self.config = ObservationConfig(weather, antenna_instance, observation_datetime)
        
        # Initialize specialized modules
        self.position_generator = SpiralPositionGenerator(self.config)
        self.file_generator = PTFFileGenerator(self.config)
        self.data_processor = SolarDataProcessor(self.config)
        self.plotter = SolarMapPlotter(self.config)
        
        # Expose commonly used attributes for backwards compatibility
        self.weather = weather
        self.antenna = antenna_instance
        self.year = self.config.year
        self.month = self.config.month
        self.day = self.config.day
        self.hour_start = self.config.hour_start
        self.minute_start = self.config.minute_start
        self.start_time = self.config.start_time
        self.file_name1 = self.config.file_name1
    
    def calculatePositions(self):
        """
        Calculate antenna positions for the spiral observation.
        
        Returns:
            tuple: (az_anten, el_anten, az_sun, el_sun, xx1, yy1, utc)
        """
        return self.position_generator.calculate_positions()
    
    def generateFile(self, path, az_anten, el_anten, utc):
        """
        Generate PTF file with antenna pointing coordinates.
        
        Args:
            path: Directory path where to save the file
            az_anten: Array of antenna azimuth coordinates
            el_anten: Array of antenna elevation coordinates
            utc: Array of UTC times in Julian Date format
            
        Returns:
            bool: True if file was successfully created
        """
        return self.file_generator.generate_ptf_file(path, az_anten, el_anten, utc)
    
    def createImages(self, fit_file_path, images_path):
        """
        Create solar maps using RBF interpolation.
        
        Args:
            fit_file_path: Path to the FITS file containing observation data
            images_path: Base directory for saving images
            
        Returns:
            str: Path where images were saved
        """
        # Calculate positions
        az_anten, el_anten, az_sun, el_sun, xx1, yy1, utc = self.calculatePositions()
        
        # Create sun position dataframe
        sun_position_df = self.data_processor.create_sun_position_dataframe(utc, xx1, yy1)
        
        # Process FITS file
        band_processed_helio_dfs = self.data_processor.process_fits_file(fit_file_path, sun_position_df)
        
        # Prepare output directory
        output_path = self.data_processor.prepare_output_directory(images_path)
        
        # Create and save solar maps
        return self.plotter.create_solar_maps_rbf(band_processed_helio_dfs, output_path)
    
    def createImages_GAUSS(self, fit_file_path, images_path):
        """
        Create solar maps using Gaussian filtering.
        
        Args:
            fit_file_path: Path to the FITS file containing observation data
            images_path: Base directory for saving images
            
        Returns:
            str: Path where images were saved
        """
        # Calculate positions
        az_anten, el_anten, az_sun, el_sun, xx1, yy1, utc = self.calculatePositions()
        
        # Create sun position dataframe
        sun_position_df = self.data_processor.create_sun_position_dataframe(utc, xx1, yy1)
        
        # Process FITS file
        band_processed_helio_dfs = self.data_processor.process_fits_file(fit_file_path, sun_position_df)
              
        # Prepare output directory
        output_path = self.data_processor.prepare_output_directory(images_path)     

        # Create and save solar maps
        return self.plotter.create_solar_maps_gaussian(band_processed_helio_dfs, output_path , rotation=True, P_angle=self.config.P_angle)