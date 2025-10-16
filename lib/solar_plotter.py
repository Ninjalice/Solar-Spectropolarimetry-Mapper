"""
Plotting and visualization module for spiral sun observations.
Handles the creation of solar maps and visualization of observational data.
"""

import numpy as np
import matplotlib.pyplot as plt
import sunpy.map
import astropy.units as u
from astropy.coordinates import SkyCoord
from scipy.interpolate import Rbf
from scipy.ndimage import gaussian_filter


class SolarMapPlotter:
    """Creates solar maps and visualizations from processed observational data."""
    
    def __init__(self, observation_config):
        """
        Initialize the plotter.
        
        Args:
            observation_config: Configuration object containing observation parameters
        """
        self.config = observation_config
        # Get bands from antenna configuration
        self.bands = getattr(self.config.antenna, 'bands', ['4.07GHZ', '6.42GHZ', '8.40GHZ', '9.80GHZ', '11.90GHZ'])
        # Get antenna diameter from antenna configuration
        self.antenna_diameter = getattr(self.config.antenna, 'diameter', 32.0)
        
    def create_solar_maps_rbf(self, band_processed_helio_dfs, output_path):
        """
        Create solar maps using RBF interpolation.
        
        Args:
            band_processed_helio_dfs: Processed data for all bands
            output_path: Directory to save the images
            
        Returns:
            str: Path where images were saved
        """
        print('Creating solar maps with RBF interpolation...')
        
        for band in self.bands:
            sun_x = band_processed_helio_dfs[band]['tx_helio_anten']
            sun_y = band_processed_helio_dfs[band]['ty_helio_anten']
            stoke_i = band_processed_helio_dfs[band][f'STOKE_I_{band}'].values
            stoke_v = band_processed_helio_dfs[band][f'STOKE_V_{band}'].values

            # Create interpolated maps
            grid_step = 10
            interp_power_stoke_i, metadata = self._create_interpolated_map_rbf(
                sun_x, sun_y, stoke_i, grid_step, 2200)
            interp_power_stoke_v, _ = self._create_interpolated_map_rbf(
                sun_x, sun_y, stoke_v, grid_step, 2200)

            # Create and save STOKE I map
            self._plot_and_save_map(
                interp_power_stoke_i, metadata, band, 'STOKE_I', 
                'gist_heat', output_path, grid_step, 2200)
            
            # Create and save STOKE V map
            self._plot_and_save_map(
                interp_power_stoke_v, metadata, band, 'STOKE_V', 
                'nipy_spectral', output_path, grid_step, 2200)
         
        plt.close('all')
        print('Solar maps created successfully!')
        return output_path
    
    def create_solar_maps_gaussian(self, band_processed_helio_dfs, output_path):
        """
        Create solar maps using Gaussian filtering.
        
        Args:
            band_processed_helio_dfs: Processed data for all bands
            output_path: Directory to save the images
            
        Returns:
            str: Path where images were saved
        """
        print('Creating solar maps with Gaussian filtering...')
        
        for band in self.bands:
            sun_x = band_processed_helio_dfs[band]['tx_helio_anten']
            sun_y = band_processed_helio_dfs[band]['ty_helio_anten']
            stoke_i = band_processed_helio_dfs[band][f'STOKE_I_{band}'].values
            stoke_v = band_processed_helio_dfs[band][f'STOKE_V_{band}'].values

            # Create interpolated maps using Gaussian filtering
            grid_step = 2
            size = 1500
            interp_power_stoke_i, metadata = self._create_interpolated_map_gaussian(
                sun_x, sun_y, stoke_i, grid_step, size, band)
            interp_power_stoke_v, _ = self._create_interpolated_map_gaussian(
                sun_x, sun_y, stoke_v, grid_step, size, band)

            # Create and save STOKE I map
            self._plot_and_save_map(
                interp_power_stoke_i, metadata, band, 'STOKE_I', 
                'gist_heat', output_path, grid_step, size)
            
            # Create and save STOKE V map
            self._plot_and_save_map(
                interp_power_stoke_v, metadata, band, 'STOKE_V', 
                'nipy_spectral', output_path, grid_step, size)
         
        plt.close('all')
        print('Solar maps created successfully!')
        return output_path
    
    def _create_interpolated_map_rbf(self, sun_x, sun_y, values, grid_step, size):
        """Create interpolated map using RBF interpolation."""
        # Define the grid covering the helioprojective coordinate space
        tx_min, tx_max = -size, size
        ty_min, ty_max = -size, size

        # Create a grid
        tx, ty = np.meshgrid(np.arange(tx_min, tx_max, grid_step),
                            np.arange(ty_min, ty_max, grid_step))

        # Interpolate power values using RBF
        rbf = Rbf(sun_x, sun_y, values, function='linear')
        interp_power = rbf(tx, ty)

        # Create metadata
        metadata = self._create_map_metadata(grid_step, tx_max, tx_min, ty_max, ty_min)
        
        return interp_power, metadata
    
    def _create_interpolated_map_gaussian(self, sun_x, sun_y, values, grid_step, size, band):
        """Create interpolated map using Gaussian filtering."""
        # Calculate beam size
        beam_size = 1.22 * (3e8 / (float(band[:-3]) * 1e9)) / self.antenna_diameter  # antenna diameter from config
        beam_size_arcsec = (beam_size * (180 / np.pi) * 3600)  # Convert to arcseconds

        # Define grid bounds
        tx_min, tx_max = -size, size
        ty_min, ty_max = -size, size

        # Create power matrix
        power_matrix = np.zeros((size, size))

        # Map coordinates to grid indices
        grid_indices_x = ((sun_x - tx_min) / grid_step).astype(int)
        grid_indices_y = ((sun_y - ty_min) / grid_step).astype(int)

        # Populate power matrix
        for i in range(len(values)):
            x_idx = grid_indices_x.iloc[i]
            y_idx = grid_indices_y.iloc[i]
            if 0 <= x_idx < size and 0 <= y_idx < size:
                power_matrix[y_idx, x_idx] += values[i]

        # Apply Gaussian filter to simulate beam size distribution
        sigma = (beam_size_arcsec / grid_step) / (2 * np.sqrt(2 * np.log(2)))  # Convert FWHM to sigma
        
        # Preserve value range while applying filter
        original_max = np.max(power_matrix)
        original_min = np.min(power_matrix)
        power_matrix = gaussian_filter(power_matrix, sigma=sigma)
        power_matrix = (power_matrix - np.min(power_matrix)) / (np.max(power_matrix) - np.min(power_matrix))
        power_matrix = power_matrix * (original_max - original_min) + original_min

        # Create metadata
        metadata = self._create_map_metadata(grid_step, tx_max, tx_min, ty_max, ty_min)
        
        return power_matrix, metadata
    
    def _create_map_metadata(self, grid_step, tx_max, tx_min, ty_max, ty_min):
        """Create metadata dictionary for sunpy map."""
        return {
            'date-obs': f'{self.config.year}-{self.config.month}-{self.config.day}T{self.config.hour_start}:{self.config.minute_start}:00',
            'crval1': 0,
            'crval2': 0,
            'cdelt1': grid_step,
            'cdelt2': grid_step,
            'cunit1': 'arcsec',
            'cunit2': 'arcsec',
            'ctype1': 'HPLN-TAN',
            'ctype2': 'HPLT-TAN',
            'crpix1': (tx_max - tx_min) / (2 * grid_step),
            'crpix2': (ty_max - ty_min) / (2 * grid_step),
            'waveunit': 'm',
            'wavelnth': 0.0262897 * u.m,
            'obsrvtry': 'Ventspils International Radio Astronomy Center',
            'detector': 'LNSP4',
            'dsun_ref': 149597870691,
            'dsun_obs': 151846026489,
            'rsun': 1573.89688496,
            'rsun_ref': 696000000,                           
            'hglt_obs': 0 * u.deg,
            'hgln_obs': 0 * u.deg,              
        }
    
    def _plot_and_save_map(self, interp_power, metadata, band, stoke_type, colormap, 
                          output_path, grid_step, size):
        """Plot and save a single solar map."""
        # Create sunpy map
        interpolated_map = sunpy.map.Map((interp_power, metadata))
        
        # Configure plotting style
        plt.ioff()
        plt.rcParams['text.color'] = 'white'
        plt.rcParams['axes.labelcolor'] = 'white'
        plt.rcParams['xtick.color'] = 'white'
        plt.rcParams['ytick.color'] = 'white'   
        
        # Create figure
        fig = plt.figure(figsize=(8, 8))
        fig.patch.set_facecolor('black')            
        ax = fig.add_subplot(projection=interpolated_map)
        
        # Plot map
        interpolated_map.plot(axes=ax, cmap=colormap)
        interpolated_map.draw_limb(axes=ax)
        interpolated_map.draw_grid(axes=ax)
        plt.colorbar(label='Power')
        
        # Add title and labels
        plt.title(f'{stoke_type} | {band} {self.config.year}-{self.config.month:02d}-{self.config.day:02d}T{self.config.hour_start:02d}:{self.config.minute_start:02d}')
        plt.xlabel('X (arcsec)')
        plt.ylabel('Y (arcsec)')
        
        # Add beam size indicator
        self._add_beam_size_circle(ax, interpolated_map, band, grid_step, size, stoke_type)
        
        plt.grid(True)
        
        # Save image
        directory = f'{self.config.year:04d}-{self.config.month:02d}-{self.config.day:02d}T{self.config.hour_start:02d}_{self.config.minute_start:02d}_00'
        filename = f'{output_path}/LNSP4-{directory}-{stoke_type}-{band}.jpeg'
        plt.savefig(filename, format='jpeg', dpi=300)
        plt.close()
    
    def _add_beam_size_circle(self, ax, interpolated_map, band, grid_step, size, stoke_type):
        """Add beam size indicator circle to the plot."""
        # Calculate beam size
        beam_size = 1.22 * (3e8 / (float(band[:-3]) * 1e9)) / self.antenna_diameter  # antenna diameter from config
        beam_size_arcsec = (beam_size * (180 / np.pi) * 3600)  # Convert to arcseconds
        
        # Position the circle
        if stoke_type == 'STOKE_I':
            circle_pos = (size-300) * u.arcsec, -(size-300) * u.arcsec
        else:  # STOKE_V
            circle_pos = 900 * u.arcsec, -900 * u.arcsec
            
        pixel_coords = interpolated_map.world_to_pixel(
            SkyCoord(circle_pos[0], circle_pos[1], frame=interpolated_map.coordinate_frame))
        
        # Create and add circle
        circle = plt.Circle((pixel_coords.x.value, pixel_coords.y.value), 
                           (beam_size_arcsec / grid_step)/2, 
                           color='white', fill=False, linestyle='--', 
                           transform=ax.transData)
        
        # Add legend
        legend_circle = plt.Line2D([0], [0], linestyle='--', color='white', label='Beam Size')
        ax.legend(handles=[legend_circle], loc='upper right', fontsize='small', frameon=False)
        ax.add_patch(circle)