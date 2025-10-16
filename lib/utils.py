import numpy as np
import pandas as pd

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import  AltAz, get_sun ,SkyCoord
from astropy.table import Table
from sunpy.coordinates import frames, sun



def RT32_SUN_PARA(utc , location):
    # Define constants
    lat_rt32 = 57.5535171694  # RT32 geographic latitude in degrees
    lon_rt32 = 21.8545525000  # RT32 geographic longitude in degrees

    # Ensure utc is an array
    if np.isscalar(utc):
        utime = np.array([utc])
    else:
        utime = np.array(utc)

    q = np.zeros(len(utc))

    for i in range(len(utime)):
        # Define observer's location
        observer_location = location

        # Create Time object
        observation_time = Time(utc[i], format='jd', location=observer_location)

        # Calculate solar alt-az coordinates
        sun_altaz = get_sun(observation_time).transform_to(AltAz(obstime=observation_time, location=observer_location))
        
        # Calculate local sidereal time
        sidereal_time = observation_time.sidereal_time('mean')

        # Calculate hour angle of the Sun
        ha = sidereal_time - sun_altaz.az
   
        # Calculate tangent of parallax angle
        q_tan = np.sin(ha) / (np.tan(lat_rt32) * np.cos(sun_altaz.alt) - 
                              np.sin(sun_altaz.alt) * np.cos(ha))

        
        # Parallax angle in degrees
        q[i] = -np.rad2deg(np.arctan(q_tan)).value # Positive values are to the left, negative to the right
        
    
    # Return result
    if len(q) == 1:
        return q[0]
    else:
        return q


def time_to_seconds(time):
    hour = time.hour
    minute = time.minute
    second = time.second
    microsecond = time.microsecond
    decimal_hour = hour + minute / 60 + second / 3600 + microsecond / 3600000000
    return decimal_hour * 3600

def seconds_to_time(year, month, day, seconds):
    hours = int(seconds // 3600)
    remaining_seconds = seconds % 3600
    minutes = int(remaining_seconds // 60)
    remaining_seconds %= 60
    remaining_seconds = int(remaining_seconds)

    time = Time(f"{year}-{month:02d}-{day:02d} {hours}:{minutes}:{remaining_seconds}")
    return time.isot

def bintable_to_pandas(file_path, hdu_number):
    try:
        # Read binary table from FITS file
        table = Table.read(file_path, hdu=hdu_number)

        # Convert table to Pandas DataFrame
        df = table.to_pandas()
       
        return df
    except Exception as e:
        print("Error:", e)
        return None
    
def bintable_to_pandas_OLD(file_path, hdu_number):
    try:
        # Read binary table from FITS file
        table = Table.read(file_path, hdu=hdu_number)        
        # Create a list to store the decompressed rows
        rows = []
        # Iterate over each column of the table
        for colname in table.colnames:
            # Check if the column contains an array
            if isinstance(table[colname][0], np.ndarray):
                # Iterate over each element of the array
                for i in range(len(table[colname][0])):
                    # Create a dictionary for the new row
                    new_row_dict = {}
                    # Iterate over each column again to fill the new row
                    for colname_inner in table.colnames:
                        if isinstance(table[colname_inner][0], np.ndarray):
                            new_row_dict[colname_inner] = table[colname_inner][0][i]
                        else:
                            new_row_dict[colname_inner] = table[colname_inner][0]
                    rows.append(new_row_dict)
            else:
                # If the column doesn't contain an array, add the original row
                row_dict = {colname: table[colname][0] for colname in table.colnames}
                rows.append(row_dict)
            # Convert the list of dictionaries to a new table
            table_descompressed = Table(rows)        
            # Convert the filtered table to a Pandas DataFrame
            df = table_descompressed.to_pandas()     

            return df 
        
    except Exception as e:
        print("Error:", e)
        return None
     
def bintable_to_pandas_OLD_BROKEN(file_path, hdu_number):
    try:
        # Read binary table from FITS file
        table = Table.read(file_path, hdu=hdu_number)     
        # Create a list to store the decompressed rows
        rows = []

        # Iterate over each row of the table
        for row_index, row in enumerate(table):
            # Create a dictionary for the original row
            row_dict = {}
            # Iterate over each column
            for colname in table.colnames:
                # Check if the column contains an array and if it's the first row
                if isinstance(row[colname], np.ndarray) and row_index == 0:
                    # If it's the first row and it's an array, add each element as an additional row
                    for i, value in enumerate(row[colname]):
                        new_row_dict = row_dict.copy()  # Copy the original row dictionary
                        new_row_dict[colname] = value  # Update the column value with the array element
                        rows.append(new_row_dict)  # Add the row to the list of rows
                else:
                    # If it's not an array or not the first row, add the value to the original row
                    row_dict[colname] = row[colname]
            # Add the original row to the list of rows
            rows.append(row_dict)

        # Convert the list of dictionaries to a new table
        table_descompressed = Table(rows)

        # Convert the filtered table to a Pandas DataFrame
        df = table_descompressed.to_pandas()       

        return df
    except Exception as e:
        print("Error:", e)
        return None

def getFinalProcessedData(observation, sunPositionDf, data_dfs):

    # Convert Julian_Time to astropy Time objects
    time = Time(sunPositionDf['UTC'], format='jd')

    # Calculate decimal hours with millisecond precision for each value
    decimal_seconds = [time_to_seconds(datetime_obj) for datetime_obj in time.datetime]

    # Add decimal hours to DataFrame
    sunPositionDf['UTC'] = np.round(np.array(decimal_seconds).astype(float), 3)

    columns = sunPositionDf.columns
    interpolated_values = np.empty((sunPositionDf.shape[0] * 1000, sunPositionDf.shape[1]))
    for i, col in enumerate(columns):
        for j in range(len(sunPositionDf)-1):
            interpolated_values[j*1000:(j+1)*1000, i] = np.linspace(sunPositionDf.iloc[j][col], sunPositionDf.iloc[j+1][col], 1001)[0:-1]

    interpolated_df = pd.DataFrame(interpolated_values, columns=columns)
    pd.set_option('display.float_format', '{:.10f}'.format)

    # Dictionary to store processed DataFrames by band
    band_processed_dfs = {}

    for band, data_df in data_dfs.items():
        data_df[f'UTC_{band}'] = data_df[f'UTC_{band}'].astype('float64').round(3)
        
        interpolated_df['UTC'] = interpolated_df['UTC'].round(3)               

        
        merged_df = pd.merge(interpolated_df, data_df, left_on='UTC', right_on=f'UTC_{band}')

        merged_df['SunX'] = merged_df['SunX'] + observation.antenna.az_offset0
        merged_df['SunY'] = merged_df['SunY'] + observation.antenna.el_offset0

      
        t_scan_cumsum = np.linspace(0, observation.num_scan - 1 , observation.num_scan) * observation.t_scan + observation.t_start_seconds

        def filter_dataframe(df, start, end):
            return df[(df['UTC'] >= start) & (df['UTC'] <= end)]

        filtered_dfs = {}
        cal_df_centre = []
        cal_df_sky = []

        columns_to_remove = ["UTC", "SunX", "SunY", f'UTC_{band}']       

        for i, start_seconds in enumerate(t_scan_cumsum):
            start_t_cal_1 = start_seconds
            end_t_cal_1 = start_seconds + observation.times_sec_dic['t_cal']

            start_t_slew_1 = start_seconds + observation.times_sec_dic['t5']
            end_t_slew_1 = start_seconds + observation.times_sec_dic['t_slew']

            start_t_cal_2 = start_seconds + observation.times_sec_dic['t_slew']
            end_t_cal_2 = start_seconds + observation.times_sec_dic['t_cal_2']

            start_t_slew_2 = start_seconds + observation.times_sec_dic['t_cal_2']
            end_t_slew_2 = start_seconds + observation.times_sec_dic['t_slew_2']

            filters = [
                (f'filter_it_{i+1}_t_cal_1', start_t_cal_1, end_t_cal_1),
                (f'filter_it_{i+1}_t_slew_2', start_t_slew_1, end_t_slew_1),
                (f'filter_it_{i+1}_t_cal_2', start_t_cal_2, end_t_cal_2),
                (f'filter_it_{i+1}_t_slew_3', start_t_slew_2, end_t_slew_2)
            ]

            for filter_name, start, end in filters:
                if filter_name.endswith("t_cal_1"):
                    cal_df_centre.append(filter_dataframe(merged_df, start, end).drop(columns=columns_to_remove).mean().values)
                elif filter_name.endswith("t_cal_2"):
                    cal_df_sky.append(filter_dataframe(merged_df, start, end).drop(columns=columns_to_remove).mean().values)
                filtered_dfs[filter_name] = filter_dataframe(merged_df, start, end)

        filtered_indices = set().union(*[filtered_df.index for filtered_df in filtered_dfs.values()])
        rest_of_df = merged_df[~merged_df.index.isin(filtered_indices)]

        cal_df_centre = np.array(cal_df_centre)
        cal_df_sky = np.array(cal_df_sky)

        sun_centre_means = np.mean(cal_df_centre, axis=0) 
        sky_means = np.mean(cal_df_sky, axis=0)      

        rest_of_df[f'RCP_{band}'] = (rest_of_df[f'RCP_{band}'] - sky_means[0]) / (sun_centre_means[0] - sky_means[0])
        rest_of_df[f'LCP_{band}'] = (rest_of_df[f'LCP_{band}'] - sky_means[0]) / (sun_centre_means[0] - sky_means[0])

        rest_of_df[f'STOKE_I_{band}'] = (rest_of_df[f'RCP_{band}'] + rest_of_df[f'LCP_{band}']) / 2
        rest_of_df[f'STOKE_V_{band}'] = (rest_of_df[f'RCP_{band}'] - rest_of_df[f'LCP_{band}']) / 2

        rest_of_df["isoT_time"] = rest_of_df.apply(lambda row: seconds_to_time(observation.year, observation.month, observation.day, row["UTC"]), axis=1)
        band_processed_dfs[band] = rest_of_df

    return band_processed_dfs


def rename_columns(data_df):
    # Rename columns by removing the first two numbers and keeping LCP and RCP prefixes
    renamed_columns = {
        col: f"{col.split()[0]} {col.split()[2]}"
        if col.startswith(('LCP', 'RCP'))
        else col
        for col in data_df.columns
    }
    data_df = data_df.rename(columns=renamed_columns)

     # Band mapping table from numbers to bands
    band_mapping = {
        '01': '4.07GHZ',
        '04': '6.42GHZ',
        '07': '8.40GHZ',
        '09': '9.80GHZ',
        '11': '11.90GHZ'
    }

    # Rename UTC columns
    renamed_columns = {
        col: f'UTC {col.split()[1]} {band_mapping[col.split()[2]]}' 
        if col.startswith('UTC') and col.split()[2] in list(band_mapping.keys())
        else col
        for col in data_df.columns
    }
    data_df = data_df.rename(columns=renamed_columns)

    return data_df

def processData(data_df):   


    data_df = rename_columns(data_df)
    
    # List of bands to process (without the first two numbers)
    bands = ['4.07GHZ', '6.42GHZ', '8.40GHZ', '9.80GHZ', '11.90GHZ']

    # Dictionary to store DataFrames for each band
    band_dfs = {}

    for band in bands:
        # UTC extraction and rounding
        UTC_RCP = np.round(data_df[f'UTC RCP {band}'].dropna() * 3600, 3)
        UTC_LCP = np.round(data_df[f'UTC LCP {band}'].dropna() * 3600, 3)
        RCP = data_df[f'RCP {band}'].dropna()
        LCP = data_df[f'LCP {band}'].dropna()

        # STOKE_I = (RCP.values + LCP.values) / 2 
        # STOKE_V = ( RCP.values - LCP.values) / 2 

        # print(STOKE_I_11_4_11_90GHZ.size)
        # print(STOKE_V_11_4_11_90GHZ.size)
        # print(UTC_RCP_11.size)      

        # Create DataFrame for the current band
        band_df = pd.DataFrame({
            f'UTC_{band}': UTC_RCP.values,          
            f'RCP_{band}': RCP,
            f'LCP_{band}': LCP
        })

        # Store in dictionary
        band_dfs[band] = band_df

    return band_dfs



def process_all_heliocentric_coordinates(band_processed_dfs, observation):
    def process_heliocentric_coordinates(band_df, observation):       

        coordsXHelio = []
        coordsYHelio = []

        # Convert the array of times into a list of Time objects
        times = [Time(t) for t in band_df['isoT_time']]

        # Calculate the sun's positions for each time in the list
        az_sun = observation.sun_location.transform_to(AltAz(obstime=times, location=observation.antenna.location, pressure=observation.weather.pressure , temperature=observation.weather.temperature, relative_humidity=observation.weather.relative_humidity ,obswl=observation.weather.obswl)).az.deg  + observation.antenna.az_offset2
        el_sun = observation.sun_location.transform_to(AltAz(obstime=times, location=observation.antenna.location, pressure=observation.weather.pressure , temperature=observation.weather.temperature, relative_humidity=observation.weather.relative_humidity ,obswl=observation.weather.obswl)).alt.deg + observation.antenna.el_offset2

        az_anten = az_sun + band_df['SunX'] / np.cos(np.deg2rad(el_sun)) / 60.
        el_anten = el_sun + band_df['SunY'] / 60.    

        band_df['az_anten'] = az_anten
        band_df['el_anten'] = el_anten

        # Iterate over the different time moments
        for index, row in band_df.iterrows():
            # Convert the AltAz coordinates from the DataFrame from degrees to radians
            el_deg = row['el_anten'] * u.deg
            az_deg = row['az_anten'] * u.deg

            # Set the observation time
            obstime = row['isoT_time']

            # Convert to heliocentric coordinates
            frame_altaz = AltAz(obstime=Time(obstime), location=observation.antenna.location, pressure=observation.weather.pressure , temperature=observation.weather.temperature, relative_humidity=observation.weather.relative_humidity ,obswl=observation.weather.obswl)
            sun_helio = SkyCoord(alt=el_deg, az=az_deg, observer='earth', distance=sun.earth_distance(obstime), frame=frame_altaz).transform_to(frames.Helioprojective)

            
            # Append the transformed coordinates to the list
            coordsXHelio.append(sun_helio.Tx.value)
            coordsYHelio.append(sun_helio.Ty.value)

        band_df['tx_helio_anten'] = coordsXHelio
        band_df['ty_helio_anten'] = coordsYHelio

        return band_df
    
    # Apply the processing function to each DataFrame in the dictionary
    for band, df in band_processed_dfs.items():
        band_processed_dfs[band] = process_heliocentric_coordinates(df, observation)
    
    return band_processed_dfs
