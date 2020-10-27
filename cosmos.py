import math
import numpy


# This class is responsible for handling the data
class cosmologicaldata( object ):
    def __init__( self, Omega_m, Omega_l, Omega_k, Hubble_h, Redshift ):
        self.data = { }
        self.data[ 'C' ] = float( 299792.458 )  # Light speed in km/s
        self.data[ 'S' ] = float( 0.01 )  # Step of the integration intervals
        self.data[ 'UAtoKM' ] = float( 149597870700 )  # Astronomical unit in km
        self.data[ 'OM' ] = float( Omega_m )  # Omega Matter
        self.data[ 'OL' ] = float( Omega_l )  # Omega Lambda
        self.data[ 'OK' ] = float( Omega_k )  # Omega k
        self.data[ 'h' ] = float( Hubble_h )  # Hubble Constant parameter
        self.data[ 'z' ] = float( Redshift )  # Redshift
        
    def calculate( self ):
        self.data[ 'H0' ] = Hubble_constant( self.data )  # Hubble Constant
        self.data[ 'DH' ] = Hubble_distance( self.data )  # Hubble distance
        self.data[ 'DC' ] = Comoving_distance( self.data )  # Line of sight comoving distance
        self.data[ 'DM' ] = Tcomoving_distance( self.data )  # Transverse comoving distance
        self.data[ 'DA' ] = Angular_diameter_distance( self.data )  # Angular diameter distance
        self.data[ 'AS' ] = Angular_scale( self.data )  # Angular scale
        self.data[ 'DL' ] = Luminosity_distance( self.data )  # Luminosity distance
        self.data[ 'TH' ] = Hubble_time( self.data )  # Hubble time
        self.data[ 'TL' ] = Lookback_time( self.data )  # Lookback time
        self.data[ 'Age' ] = Age( self.data )  # Universe age
        self.data[ 'TE' ] = Emission_age( self.data )  # Galaxy Age
        
    def calculate_globals( self ):
        ''' separate the common function that only need to be calculated once.
            Ensure that these are called after calculate since we are depended on the 
            parameters '''
        self.data[ 'H0' ] = Hubble_constant( self.data )  # Hubble Constant
        self.data[ 'TH' ] = Hubble_time( self.data )  # Hubble time
        return Age( self.data )
        
    def __getitem__( self, key ): return self.data[ key ]


# Down there are the calculationroutines based on the CosmoData Class

# Calculates the Hubble constant given the h parameter
def Hubble_constant( Data ):
    Ho = 100 * Data[ 'h' ]
    return Ho  # Hubble constant in Km/s/Mpc


# Calculates the hubble distance for a given Hubble constant and light speed
def Hubble_distance( Data ):
    Dh = (Data[ 'C' ] / Data[ 'H0' ])
    return Dh  # Hubble distance in meters


# Calculates an E(z) function value for the given parameters
def E( z, Data ):
    _E = math.sqrt( Data[ 'OM' ] * pow( (1 + z), 3 ) + Data[ 'OK' ] * pow( (1 + z), 2 ) + Data[ 'OL' ] )
    return _E


# Calculates the line-of-sight comoving distance
def Comoving_distance( Data ):
    s = 0
    z = 0
    # Algorithm to calculate the integral by diving it into small trapezes with height given by the step
    while (z + Data[ 'S' ] <= Data[ 'z' ]):
        s += (((1 / E( z, Data )) + (1 / E( (z + Data[ 'S' ]), Data ))) * (Data[ 'S' ])) / 2
        z += Data[ 'S' ]
    Dc = Data[ 'DH' ] * s
    return Dc  # Comoving distance in ???


# Calculates the transverse comoving distance
def Tcomoving_distance( Data ):
    Dm = 0
    if float( Data[ 'OK' ] ) > 0:
        Dm = Data[ 'DH' ] * (1 / math.sqrt( Data[ 'OK' ] )) * numpy.sinh(
            math.sqrt( Data[ 'OK' ] ) * (Data[ 'DC' ] / Data[ 'DH' ]) )
    elif float( Data[ 'OK' ] ) == 0:
        Dm = Data[ 'DC' ]
    elif float( Data[ 'OK' ] ) < 0:
        Dm = Data[ 'DH' ] * (1 / math.sqrt( Data[ 'OK' ] ) * (-1)) * numpy.sin(
            math.sqrt( Data[ 'OK' ] * (-1) ) * (Data[ 'DC' ] / Data[ 'DH' ]) )
    return Dm  # Transverse comoving distance in ???


# Calculates the Angular diameter distance
def Angular_diameter_distance( Data ):
    Da = Data[ 'DM' ] / (1 + Data[ 'z' ])
    return Da  # Angular diameter distance in Mpc/rad


# Calculates the Angular Scale. Basically, converts the angular diameter distance into kpc/arcsec
def Angular_scale( Data ):
    As = Data['DA'] * 0.0048481368  # Angular scale in Kpc/arcsec
    return As


# Calculates the luminosity distance
def Luminosity_distance( Data ):
    Dl = (1 + Data[ 'z' ]) * Data[ 'DM' ]
    return Dl  # Luminosity distance in Mpc


# Calculates the lookback time
def Lookback_time( Data ):
    s = 0
    z = 0
    # Algorithm to calculate the integral by diving it into small trapezes with height given by the step
    while (z + Data[ 'S' ] <= Data[ 'z' ]):
        s += ((1 / (1 + z) * (1 / E( z, Data )) + 1 / (1 + z + Data[ 'S' ]) * (1 / E( (z + Data[ 'S' ]), Data ))) * (
        Data[ 'S' ])) / 2
        z += Data[ 'S' ]
    Tl = Data[ 'TH' ] * s
    return Tl  # Lookback time in Gyr


# Calculates the approximate age using the same method used to calculate the lookback time but now with a really high redshift
def Age( Data ):
    s = 0
    z = 0
    # Algorithm to calculate the integral by diving it into small trapezes with height given by the step
    while (z + Data[ 'S' ] <= 1200):
        s += ((1 / (1 + z) * (1 / E( z, Data )) + 1 / (1 + z + Data[ 'S' ]) * (1 / E( (z + Data[ 'S' ]), Data ))) * (
        Data[ 'S' ])) / 2
        z += Data[ 'S' ]
    Age = Data[ 'TH' ] * s
    return Age  # Age in Gyr


# Calculates the Hubble time
def Hubble_time( Data ):
    Th = 977.8131056/Data['H0']  # Lookback time in gyr
    return Th  # Hubble time in Gyr


# Calculates the age of the universe when the photons from a star at a given redshift were emmited
def Emission_age( Data ):
    Te = Data[ 'Age' ] - Data[ 'TL' ]
    return Te