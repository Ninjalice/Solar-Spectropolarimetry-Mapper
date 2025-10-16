"""
File generation module for spiral sun observations.
Handles the creation of PTF (pointing) files for antenna control.
"""

from astropy.time import Time


class PTFFileGenerator:
    """Generates PTF files for antenna pointing during solar observations."""
    
    def __init__(self, observation_config):
        """
        Initialize the file generator.
        
        Args:
            observation_config: Configuration object containing observation parameters
        """
        self.config = observation_config
        
    def generate_ptf_file(self, path, az_anten, el_anten, utc):
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
        file_path = path + self.config.file_name1
        
        with open(file_path, 'w') as file1:
            # Write header information
            self._write_header(file1)
            
            # Write observation details
            self._write_observation_details(file1)
            
            # Write PTF format sections
            self._write_ptf_sections(file1)
            
            # Write coordinate data
            self._write_coordinate_data(file1, az_anten, el_anten, utc)

        print('-------------------------------------------------------------')
        print('Saved: ', self.config.file_name1, "  ",  len(utc), "  points")
        return True
    
    def _write_header(self, file):
        """Write the PTF file header."""
        file.write("# Table of horizontal coordinates for Sun scanning\n")
        file.write("# Celestial Object ID: SUN\n")
        file.write("# Spiral scanning, number of scans: " + str(self.config.num_scan) + "\n")
        file.write("# Scan time schedule:\n")
        file.write("# Sun center calibration      60 sec\n")
        file.write("# 1 spiral turn              40 sec  step 0->6 arcmin\n")
        file.write("# 2 spiral                   100 sec  step 6->12 arcmin\n")
        file.write("# 3 spiral                   120 sec step 12->18 arcmin\n")
        file.write("# 4 spiral                   120 sec step 18-28 arcmin\n")
        file.write("# 5 spiral                   120 sec step 28-40 arcmin\n")
        file.write("# Slew to sky 4 Rsun        20 sec\n")
        file.write("# Sky calibration           60 sec\n")
        file.write("# Slew to Sun center        20 sec\n")
    
    def _write_observation_details(self, file):
        """Write observation-specific details."""
        file.write(f"# Az@Start Time      : {self.config.az_start:.5f} deg.\n")
        file.write(f"# El@Start Time      : {self.config.el_start:.5f} deg.\n")
        file.write(f"# Start Time (UTC+0) : {self.config.start_time.isot}\n")        
        file.write(f"# End Time, (UTC), hours: {Time(self.config.utime_end, format='jd').isot}\n")
        file.write(f"# Mean observation time (UTC): {Time(self.config.utime_mean, format='jd').isot}\n")
        file.write(f"# Az offset          : {self.config.antenna.az_offset0} deg.\n")
        file.write(f"# El offset          : {self.config.antenna.el_offset0} deg.\n")
        file.write("#!!! Offset of the antenna pointing system must be set 0 !!!\n")
    
    def _write_ptf_sections(self, file):
        """Write PTF format control sections."""
        file.write("\n[Interpolation]\n")
        file.write("Newton\n")
        file.write("\n[Load Mode]\n")
        file.write("New\n")
        file.write("\n[Start Time]\n")
        file.write(f"{self.config.start_time.iso}\n")
        file.write("\n[Table Data]\n")
        file.write("#  Time   Az source [degr.]     El source[degr.]\n")
    
    def _write_coordinate_data(self, file, az_anten, el_anten, utc):
        """Write the time-coordinate data table."""
        for i in range(len(utc)):
            file.write(f"{Time(utc[i], format='jd').isot} \t\t {az_anten[i]:.5f} \t\t {el_anten[i]:.5f}\n")