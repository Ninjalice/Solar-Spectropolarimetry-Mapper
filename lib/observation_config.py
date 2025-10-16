"""
Configuration module for spiral sun observations.
Contains observation parameters and configuration management.
"""

from datetime import datetime
import numpy as np
from astropy.time import Time
from astropy.coordinates import AltAz, get_sun
from lib.utils import time_to_seconds


class ObservationConfig:
    """Configuration class that holds all observation parameters."""
    
    def __init__(self, weather, antenna_instance, observation_datetime: datetime):
        """
        Initialize observation configuration.
        
        Args:
            weather: Weather instance with atmospheric conditions
            antenna_instance: Antenna instance with location and properties
            observation_datetime: Start time of the observation
        """
        print("ObservationConfig", observation_datetime.isoformat())
        
        # Store instances
        self.weather = weather
        self.antenna = antenna_instance
        
        # Extract datetime components
        self.year = observation_datetime.year
        self.month = observation_datetime.month
        self.day = observation_datetime.day
        self.hour_start = observation_datetime.hour
        self.minute_start = observation_datetime.minute

        # Time calculations
        self.start_time = Time(observation_datetime)
        self.t_start_seconds = time_to_seconds(self.start_time.datetime)

        # Sun position at start time
        self.sun_location = get_sun(self.start_time)
        sun_altaz = self.sun_location.transform_to(
            AltAz(obstime=self.start_time, 
                  location=self.antenna.location, 
                  pressure=self.weather.pressure,
                  temperature=self.weather.temperature, 
                  relative_humidity=self.weather.relative_humidity,
                  obswl=self.weather.obswl)
        )
        self.az_start = sun_altaz.az
        self.el_start = sun_altaz.alt

        # File naming
        self.file_out = f"{(self.year - 2000):02d}{self.month:02d}{self.day:02d}_{self.hour_start:02d}{self.minute_start:02d}"
        self.file_name1 = f"sun_scan_{self.file_out}.ptf"

        # Scanning parameters
        self._init_scanning_parameters()
        
        # Time calculations
        self._init_time_parameters()
    
    def _init_scanning_parameters(self):
        """Initialize spiral scanning parameters."""
        self.num_scan = 5
        self.step1 = 6.  # arcmin
        self.step2 = 6.
        self.step3 = 6.
        self.step4 = 10.
        self.step5 = 12.
        self.sky = 16. * 4  # arcmin

        # Timing parameters (seconds)
        self.t_cal = 60.  # calibration time
        self.t1 = 40.     # first spiral
        self.t2 = 100.    # second spiral
        self.t3 = 120.    # third spiral
        self.t4 = 120.    # fourth spiral
        self.t5 = 120.    # fifth spiral
        self.t_slew = 20. # slew time
        
        # Total scan time
        self.t_scan = (self.t_cal + self.t1 + self.t2 + self.t3 + 
                      self.t4 + self.t5 + self.t_slew + self.t_cal + self.t_slew)

        # Arrays for convenience
        self.t_spirals = np.array([self.t1, self.t2, self.t3, self.t4, self.t5])
        self.times = np.array([self.t_cal, self.t1, self.t2, self.t3, self.t4, 
                              self.t5, self.t_slew, self.t_cal, self.t_slew])
        self.labels = np.array(['t_cal', 't1', 't2', 't3', 't4', 't5', 
                               't_slew', 't_cal_2', 't_slew_2'])
        self.times_scan_cumsum = np.cumsum(self.times)
        self.times_sec_dic = dict(zip(self.labels, self.times_scan_cumsum))
    
    def _init_time_parameters(self):
        """Initialize time-related parameters."""
        self.utime_start = self.start_time.jd  # start time in JD
        self.time_session = self.num_scan * self.t_scan / 3600. / 24  # session duration in days
        self.utime_mean = self.utime_start + self.time_session / 2.  # mean UT time
        self.utime_end = self.utime_start + self.time_session  # end UT time