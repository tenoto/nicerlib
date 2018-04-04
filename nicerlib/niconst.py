NUM_OF_MPU = 7
NUM_OF_FPM = 8
NICER_PI_MIN = 0
NICER_PI_MAX = 1500

# Conversion from PI to keV (PI is in units of 10 eV)
PI_TO_KEV = 10.0/1000.0
KEV_TO_PI = 1000.0/10.0

"""
EVENT_FLAGS == xxxxx1: "undershoot" reset
EVENT_FLAGS == xxxx1x: "overshoot" reset
EVENT_FLAGS == xxx1xx: soGware sample
EVENT_FLAGS == xx1xxx: fast signal chain triggered 
EVENT_FLAGS == x1xxxx: slow signal chain triggered 
EVENT_FLAGS == 1xxxxx: first event in MPU packet
"""