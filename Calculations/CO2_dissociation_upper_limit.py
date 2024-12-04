import numpy as np

max_peak_area = 0.15 #this is after dosing for a long time, corresponding to 0.33 coverage? The actual maximum is more like 0.22. this is not consistent with earlier measurements though
max_peak_area = 0.3 #guess for upper limit (it can not be lower than this as this was the actual peak area I measured)
max_coverage = 0.33 #This is for dosing at 100K?
CO_s0 = 0.5 #estimate. keep in mind that the actual sticking probability varies with time, may want to write a function for that
cu_surface_density = 1.5E15 #for palladium, could not find cu yet
pumping_speed = 315 #L/s, from thesis Saskia



#7 CO2 + 12 He
pressure = 1E-7 * 1E-3 #bar. Estimate because QMS data / pressure data are inconsistent for mixed beams

orifice_diam = 0.62 #cm, from LF manual
orifice_area = np.pi*orifice_diam**2

peak_area = 0.3 
peak_area_difference = 0.02 #maximum
dosing_time = 110 * 60 #seconds

# 1 CO2 + 14 He
pressure =  2E-8 * 1E-3 #bar. Estimate because QMS data / pressure data are inconsistent for mixed beams

# 0.5 CO2 + 7 He
pressure =  1E-8 * 1E-3 #bar. Estimate because QMS data / pressure data are inconsistent for mixed beams

# 0.5 CO2 + 7 He
pressure =  1E-8 * 1E-3 #bar. Estimate because QMS data / pressure data are inconsistent for mixed beams

#17 CO2
pressure = 4E-8 * 1E-3

orifice_diam = 0.62 #cm, from LF manual
orifice_area = np.pi*orifice_diam**2

peak_area = 0.3 
peak_area_difference = 0.02 #maximum
dosing_time = 110 * 60 #seconds

"""
Notes

reactivity = surface coverage / total dose = surface coverage / (flux * dosing time)
flux = surface coverage / (reactivity * total time)
"""

def main():
    print('main')

    co2_flux = flux(pressure, pumping_speed, orifice_area)
    print('flux', "{0:.2E}".format(co2_flux))
    print('flux', co2_flux/cu_surface_density, 'ML/s')
    co_coverage = surface_coverage(peak_area, max_peak_area, max_coverage) * cu_surface_density
    print('co_coverage', co_coverage/cu_surface_density)
    co_contamination = co_contamination_level(co_coverage, CO_s0, dosing_time, co2_flux) #took s0 for CO because there was no indication of the adsorption slowing down

    co_coverage_difference = surface_coverage(peak_area_difference, max_peak_area, max_coverage) * cu_surface_density
    co2_reactivity_difference = reactivity(co_coverage_difference, dosing_time, co2_flux)

    print('co contamination', co_contamination)
    print('co2 reactivity difference', co2_reactivity_difference)


def reactivity(coverage, dosing_time, flux):
    print('calculating max_reactivity..')
    return coverage / (flux * dosing_time)


def co_contamination_level(surface_coverage, reactivity, dosing_time, flux): #average reactivity
    print('calculating co_contamination_level..')
    co_flux = surface_coverage / (reactivity*dosing_time)
    return co_flux / flux


def surface_coverage(peak_area, max_peak_area, max_coverage):
    """
    Returns fractional surface coverage
    """
    return peak_area/max_peak_area*max_coverage


def flux(pressure, pumping_speed, area, T=300): 
    """
    Returns number density
    """
    print('calculating flux..')
    R = 0.083 # L bar / mol / K
    avo = 6.022E23
    R /= avo
    return pressure*pumping_speed/(R*T*area)


main()