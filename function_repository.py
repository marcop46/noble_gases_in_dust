import numpy as np



# Verchovsky et al. 2003 PASA
# cross section (cm^2) of the grain for ion capture (equation 1) - using ratio_d_s = d/s
# d = grain diameter (I will use micrometers, but does not really matter)
# s = implantation range

def cross_section (d, ratio_d_s):
    sigma = np.pi / 4.0e0 * d**2 * (1.0 - 1.0/ratio_d_s**2)
    return(sigma)

# concentration (cm^-3), so here I need to be careful the units I am using for d! (Equation 2)
# phi = fluence of ions (cm^-2)
# so, d and phi has to use the same metric scale
def concentration(phi,d,ratio_d_s):
    c = phi * 3.0e0/2.0e0 * (d**2 - d**2/ratio_d_s**2)/d**3
    return(c)

# fraction of ions which stop inside the grain (equation 3)
def fraction_ion_stopped(d,ratio_d_s):
    dum  = cross_section(d,ratio_d_s)
    frac = dum/(np.pi/4.0e0 * d**2)
    return(frac)


# Verchovsky et al. 2003 PASA
# grain size distribution function (equation 10)
def grain_size_distribution(sigma_d, d, average_d):
    f_d = 1.0e0/(sigma_d*np.sqrt(2.*np.pi)) * np.exp(-(d - average_d)**2/(2.0e0*sigma_d**2))
    return(f_d)

# distribution function for implantation energies (or better of implantation ranges) - equation 11
def implantation_range_distribution(sigma_s,s,average_s):
    e_s = 1.0e0/(sigma_s*np.sqrt(2.*np.pi)) * np.exp(-(s - average_s)**2/(2.0e0*sigma_s**2))
    return(e_s)

# concentration integrated (equation 9)
def concentration_integrated(phi, sigma_d,d, average_d,sigma_s,s,average_s, step):
    integral = 0.
    for i in range(len(d)):
        for j in range(len(s)):
            if s[j] <= d[i]: 
                c = phi * 3.0e0/2.0e0 * (d[i]**2 - s[j]**2)/d[i]**3
                e = implantation_range_distribution(sigma_s,s[j],average_s)
                g = grain_size_distribution(sigma_d, d[i], average_d)
                integral = integral + e*g*c * step**2
    return(integral)


def get_g_and_n_components(kj_transpose, ratio_s, ratio_n, ind_, ind_d, factor_concentration):
    
    
    # Fit the regression line for each isotope
    slope_derived = []
    coeff_offset  = [] 
    errors        = []
    for i in range(len(kj_transpose)):
        (m, b), (SSE,), *_ = np.polyfit(kj_transpose[ind_], kj_transpose[i], deg=1, full=True)
        slope_derived.append(m); coeff_offset.append(b); errors.append(SSE)
        
    slope_derived_arr = np.array(slope_derived); coeff_offset_arr = np.array(coeff_offset)
    
    # printing for debugging
    #for i in range(len(kj_xe_transp_xe130)):
    #    print(isotopes_xe[i],', slope=',"%.4e" % slope_derived_arr[i],\
    #          ', coeff=',"%.4e" % coeff_offset_arr[i],\
    #          ',err=',"%.4e" % errors[i])
    
    
    kj_g  = [slope_derived_arr[i] * ratio_s + coeff_offset_arr[i] \
                       for i in range(len(slope_derived_arr))]
    kj_n  = [slope_derived_arr[i] * ratio_n + coeff_offset_arr[i] \
                       for i in range(len(slope_derived_arr))]  
    
    # printing for debugging
    #for i in range(len(kj_g)):
    #    print(isotopes_xe[i], 'G=', "%.4e" % kj_xe_g[i],', N=',"%.4e" % kj_xe_n[i])
    
        
    # for each KJ group, isotopic concentration G and N 
    kj_g_conc = [];kj_n_conc = []
    for i in range(len(kj_g)):
        dum = [kj_g[i] * factor_concentration[j] for j in range(len(factor_concentration))]
        kj_g_conc.append(dum)
        dum = [kj_n[i] * factor_concentration[j] for j in range(len(factor_concentration))]
        kj_n_conc.append(dum)
    
    
    return(kj_g, kj_g_conc, kj_n, kj_n_conc, slope_derived_arr, coeff_offset_arr)


def read_lewis94_data(file_name,element):
    
    # read data from Lewis+1994
    f_ = open(file_name,'r')
    #
    if element == 'xe':
        f_.readline()
        header = f_.readline()
        data = f_.readlines()
        f_.close()
        
        # extra info and most likely we will not need. Reading here for now
        mean_xe   = [float(data[0].split()[i+1].strip()) for i in range(len(header.split('&'))-1)]
        slope_xe  = [float(data[1].split()[i+1].strip()) for i in range(len(header.split('&'))-1)]
        xe_G      = [float(data[2].split()[i+1].strip()) for i in range(len(header.split('&'))-1)]
        xe_N      = [float(data[3].split()[i+1].strip()) for i in range(len(header.split('&'))-1)]
        # reading grain data for different KJ group
        kja = [float(data[4].split()[i+1].strip()) for i in range(len(header.split('&'))-1)]
        kjb = [float(data[5].split()[i+1].strip()) for i in range(len(header.split('&'))-1)]
        kjc = [float(data[6].split()[i+1].strip()) for i in range(len(header.split('&'))-1)]
        kjd = [float(data[7].split()[i+1].strip()) for i in range(len(header.split('&'))-1)]
        kje = [float(data[8].split()[i+1].strip()) for i in range(len(header.split('&'))-1)]
        kjf = [float(data[9].split()[i+1].strip()) for i in range(len(header.split('&'))-1)]
        kjg = [float(data[10].split()[i+1].strip()) for i in range(len(header.split('&'))-1)]
        kjh = [float(data[11].split()[i+1].strip()) for i in range(len(header.split('&'))-1)]
    elif element == 'kr':
        f_.readline()
        header = f_.readline(); f_.readline()
        data = f_.readlines()
        f_.close()
        
        # reading grain data for different KJ group
        kja = [float(data[0].split()[i+1].strip()) for i in range(len(header.split('&'))-1)]
        kjb = [float(data[1].split()[i+1].strip()) for i in range(len(header.split('&'))-1)]
        kjc = [float(data[2].split()[i+1].strip()) for i in range(len(header.split('&'))-1)]
        kjd = [float(data[3].split()[i+1].strip()) for i in range(len(header.split('&'))-1)]
        kje = [float(data[4].split()[i+1].strip()) for i in range(len(header.split('&'))-1)]
        kjf = [float(data[5].split()[i+1].strip()) for i in range(len(header.split('&'))-1)]
        kjg = [float(data[6].split()[i+1].strip()) for i in range(len(header.split('&'))-1)]
        kjh = [float(data[7].split()[i+1].strip()) for i in range(len(header.split('&'))-1)]


    # isotope names and abundance info per KJ group    
    isotopes  = [header.split('&')[i].strip() for i in range(len(header.split('&')))][1:]

    # we may want to allow to load separate samples for different elements? 
    # Lewis+1994 they do, but not fully explained why. Actually, slopes changes if 
    # not all groups of grains are used.... to be explored
    if element == 'xe':
        kj = [kja, kjb, kjc, kjd, kje, kjf, kjg, kjh]
    elif element == 'kr':
        kj = [kja, kjb, kjc, kjd, kje, kjf, kjg, kjh]
    
    if element == 'xe':
        # getting the xe130/xe132 per group and xe132 concentration
        kj_xe30dxe32   = np.array([float(data[i].split()[-2]) for i in range(4,len(kj)+4)])
        kj_xe32        = np.array([float(data[i].split()[-1]) for i in range(4,len(kj)+4)])
        # factor to use and get G and N isotopic concentrations from ratios
        factor_conc = kj_xe30dxe32 * kj_xe32
    elif element == 'kr':
        kj_kr82 = np.array([float(data[i].split()[-1]) for i in range(len(kj))])
        # factor to use and get G and N isotopic concentrations from ratios
        factor_conc = kj_kr82
    # as given normalized in the file
    kj_array        = np.array(kj)
    kj_transp       = np.transpose(kj_array)
    if element == 'xe':
        # I need to use normalized to xe130, so renomarlizing everything to xe130/xe132 
        kj_transp = [kj_transp[i]/kj_transp[4] for i in range(len(kj_transp))]
        
        
    return(isotopes, kj, kj_array, kj_transp, factor_conc)

