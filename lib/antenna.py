from astropy.coordinates import EarthLocation

class Antenna:
    def __init__(self, az_offset0=-0.12, el_offset0=0.13, az_offset2=0.05, el_offset2=-0.05, description="Default Antenna" , bands=None, diameter=None):
        if bands is None or diameter is None:
            raise ValueError("Bands and diameter must be provided for the antenna.")
        
        self.az_offset0 = az_offset0
        self.el_offset0 = el_offset0
        self.az_offset2 = az_offset2
        self.el_offset2 = el_offset2
        self.latitude = None
        self.longitude = None
        self.elevation = None
        self.location = None
        self.description = description
        self.bands = []
        self.diameter = None  # in meters
        
    def set_location(self, latitude, longitude, elevation):
        self.latitude = latitude
        self.longitude = longitude
        self.elevation = elevation
        self.location = EarthLocation(lat=self.latitude, lon=self.longitude, height=self.elevation)


class RT32(Antenna):
    def __init__(self, description="RT-32 Antenna"):
        super().__init__(description=description , bands=['4.07GHZ', '6.42GHZ', '8.40GHZ', '9.80GHZ', '11.90GHZ'], diameter=32.0)
    