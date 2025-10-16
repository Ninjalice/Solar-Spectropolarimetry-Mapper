import pandas as pd
import astropy.units as u
from lib.antenna import RT32
from lib.weather import Weather
from lib.spiral_sun_observation import SpiralSunObservation
import warnings
from datetime import datetime

# Imprimimos estadísticas resumidas del DataFrame final
pd.set_option('display.float_format', '{:.10f}'.format)
warnings.simplefilter(action='ignore', category=pd.errors.SettingWithCopyWarning)

# RT32 location (Ventspils, Latvia)
rt32_antenna = RT32()
rt32_antenna.az_offset0 = 0
rt32_antenna.el_offset0 = 0

rt32_antenna.set_location(latitude=57.5535171694, longitude=21.8545525000, elevation=20)

# Define constants
year = 2025
month = 5
day = 8
hour_start = 12
minute_start = 1

temperature = u.Quantity(15.0, unit=u.deg_C)
pressure = u.Quantity(1013.25, unit=u.hPa)
relative_humidity = u.Quantity(60.0, unit=u.percent)
obswl =u.Quantity(50000, unit=u.nm) 

weather = Weather(temperature, pressure, relative_humidity, obswl)

observation_datetime = datetime(year, month, day, hour_start, minute_start)
observation = SpiralSunObservation(weather, rt32_antenna, observation_datetime)

fit_file_path = "data/examples/FITS/lnsp4_5ch_250508_120000_125610.fit"
images_path = "."

observation.createImages_GAUSS(fit_file_path , images_path)