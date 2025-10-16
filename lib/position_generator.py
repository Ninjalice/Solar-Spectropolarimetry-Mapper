"""
Position generation module for spiral sun observations.
Handles the calculation of spiral trajectories and coordinate transformations.
"""

import numpy as np
from astropy.time import Time
from astropy.coordinates import AltAz, get_sun
from lib.utils import RT32_SUN_PARA


class SpiralPositionGenerator:
    """Generates spiral scanning positions for solar observations."""
    
    def __init__(self, observation_config):
        """
        Initialize the position generator.
        
        Args:
            observation_config: Configuration object containing observation parameters
        """
        self.config = observation_config
        
    def calculate_spiral_coordinates(self):
        """
        Calculate the basic spiral coordinates (x, y) for one scan.
        
        Returns:
            tuple: (x, y) arrays containing the spiral coordinates in arcminutes
        """
        x = np.zeros(700)
        y = np.zeros(700)

        # Calibration at Sun center
        for i in range(int(self.config.t_cal)):
            x[i] = 0.
            y[i] = 0.

        # First spiral turn
        for i in range(int(self.config.t1)):
            i0 = int(self.config.t_cal) 
            fi = i * 360. / self.config.t1
            r = i * self.config.step1 / self.config.t1
            x[int(i + i0)] = r * np.cos(np.deg2rad(fi))
            y[int(i + i0)] = r * np.sin(np.deg2rad(fi))

        # Second spiral turn
        for i in range(int(self.config.t2)):
            i0 = int(self.config.t_cal) + int(self.config.t1)
            fi = i * 360. / self.config.t2
            r = self.config.step1 + i * (self.config.step2 / self.config.t2)
            x[i + i0] = r * np.cos(np.deg2rad(fi))
            y[i + i0] = r * np.sin(np.deg2rad(fi))

        # Third spiral turn
        for i in range(int(self.config.t3)):
            i0 = int(self.config.t_cal) + int(self.config.t1) + int(self.config.t2)
            fi = i * 360. / self.config.t3
            r = self.config.step1 + self.config.step2 + i * (self.config.step3 / self.config.t3)
            x[i + i0] = r * np.cos(np.deg2rad(fi))
            y[i + i0] = r * np.sin(np.deg2rad(fi))

        # Fourth spiral turn
        for i in range(int(self.config.t4)):
            i0 = int(self.config.t_cal) + int(self.config.t1) + int(self.config.t2) + int(self.config.t3)
            fi = i * 360. / self.config.t4
            r = self.config.step1 + self.config.step3 + self.config.step3 + (i * self.config.step4 / self.config.t4)
            x[i + i0] = r * np.cos(np.deg2rad(fi))
            y[i + i0] = r * np.sin(np.deg2rad(fi))

        # Fifth spiral turn
        for i in range(int(self.config.t5)):
            i0 = int(self.config.t_cal) + int(self.config.t1) + int(self.config.t2) + int(self.config.t3) + int(self.config.t4)
            fi = i * 360. / self.config.t5
            r = self.config.step1 + self.config.step2 + self.config.step3 + self.config.step4 + i * (self.config.step5 / self.config.t5)
            x[i + i0] = r * np.cos(np.deg2rad(fi))
            y[i + i0] = r * np.sin(np.deg2rad(fi))

        # Slew to calibration sky position
        for i in range(int(self.config.t_slew)):
            i0 = int(self.config.t_cal) + int(self.config.t1) + int(self.config.t2) + int(self.config.t3) + int(self.config.t4) + int(self.config.t5)
            y[i + i0] = 0.
            x0 = self.config.step1 + self.config.step2 + self.config.step3 + self.config.step4 + self.config.step5
            x[i + i0] = x0 + i * (self.config.sky - x0) / self.config.t_slew

        # Calibration at sky position
        for i in range(int(self.config.t_cal)):
            i0 = int(self.config.t_cal) + int(self.config.t1) + int(self.config.t2) + int(self.config.t3) + int(self.config.t4) + int(self.config.t5) + int(self.config.t_slew)
            x[i + i0] = self.config.sky
            y[i + i0] = 0.

        # Slew back to Sun center
        for i in range(int(self.config.t_slew)):
            i0 = int(self.config.t_cal) + int(self.config.t1) + int(self.config.t2) + int(self.config.t3) + int(self.config.t4) + int(self.config.t5) + int(self.config.t_slew) + int(self.config.t_cal)
            y[i + i0] = 0.
            x[i + i0] = (int(self.config.t_slew) - i - 1) * self.config.sky / self.config.t_slew

        return x, y
    
    def generate_multiple_scans(self, x, y):
        """
        Generate coordinates for multiple spiral scans rotated around the Sun.
        
        Args:
            x, y: Single scan coordinates
            
        Returns:
            tuple: (xx, yy) arrays containing all scan coordinates
        """
        xx = np.zeros(self.config.num_scan * int(self.config.t_scan))
        yy = np.zeros(self.config.num_scan * int(self.config.t_scan))

        for j in range(self.config.num_scan):
            ff = j * (360.0 / self.config.num_scan)
            for i in range(int(self.config.t_scan)):
                ii = j * int(self.config.t_scan) + i
                xx[ii] = x[i] * np.cos(np.deg2rad(ff)) - y[i] * np.sin(np.deg2rad(ff))
                yy[ii] = x[i] * np.sin(np.deg2rad(ff)) + y[i] * np.cos(np.deg2rad(ff))
                
        return xx, yy
    
    def transform_to_antenna_coordinates(self, xx, yy):
        """
        Transform spiral coordinates to antenna azimuth and elevation.
        
        Args:
            xx, yy: Spiral coordinates in arcminutes
            
        Returns:
            tuple: (az_anten, el_anten, az_sun, el_sun, xx1, yy1, utc)
        """
        # Generate time array
        utc = self.config.utime_start + np.arange(self.config.num_scan * self.config.t_scan) / 3600. / 24   
        
        # Get parallactic angle correction
        q = RT32_SUN_PARA(utc, self.config.antenna.location)

        # Apply parallactic angle rotation
        xx1 = xx * np.cos(np.deg2rad(q)) - yy * np.sin(np.deg2rad(q))
        yy1 = xx * np.sin(np.deg2rad(q)) + yy * np.cos(np.deg2rad(q))

        # Calculate Sun's position at each time
        az_sun = self.config.sun_location.transform_to(
            AltAz(obstime=Time(utc, format='jd'), 
                  location=self.config.antenna.location, 
                  pressure=self.config.weather.pressure,
                  temperature=self.config.weather.temperature, 
                  relative_humidity=self.config.weather.relative_humidity,
                  obswl=self.config.weather.obswl)
        ).az.deg
        
        el_sun = self.config.sun_location.transform_to(
            AltAz(obstime=Time(utc, format='jd'), 
                  location=self.config.antenna.location, 
                  pressure=self.config.weather.pressure,
                  temperature=self.config.weather.temperature, 
                  relative_humidity=self.config.weather.relative_humidity,
                  obswl=self.config.weather.obswl)
        ).alt.deg

        # Calculate antenna pointing coordinates
        az_anten = az_sun + xx1 / np.cos(np.deg2rad(el_sun)) / 60.
        el_anten = el_sun + yy1 / 60.    

        return az_anten, el_anten, az_sun, el_sun, xx1, yy1, utc
    
    def calculate_positions(self):
        """
        Calculate all antenna positions for the spiral observation.
        
        Returns:
            tuple: (az_anten, el_anten, az_sun, el_sun, xx1, yy1, utc)
        """
        # Generate basic spiral coordinates
        x, y = self.calculate_spiral_coordinates()
        
        # Generate multiple rotated scans
        xx, yy = self.generate_multiple_scans(x, y)
        
        # Transform to antenna coordinates
        return self.transform_to_antenna_coordinates(xx, yy)