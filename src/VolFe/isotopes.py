# isotopes.py

import pandas as pd
import numpy as np

import VolFe.melt_gas as mg
import VolFe.calculations as c

######################
### delta notation ###
######################

def delta_standards(standard,isotope,element):
    if element == "S":
        if standard == "VCDT":
            if isotope == 34:
                reference = 1/22.6436 # 34S/32S Ding et al. (2001)
    elif element == "C":
        if standard == "VPBD":
            if isotope == 13:
                reference = 0.01123720 # 13C/12C International Atomic Energy Agency (1995)
    elif element == "H":
        if standard == "VSMOW":
            if isotope == 2:
                reference = 155.76/1.e6 # 2H/1H Hagemann et al. (1970)
    return reference
 
def ratio2delta(standard,isotope,ratio,element):
    reference = delta_standard(standard,isotope,element)
    d = ((ratio - reference)/reference)*1000.
    return d

def delta2ratio(standard,isotope,d):
    reference = delta_standard(standard,isotope,element)
    ratio = ((d/1000.)*reference) + reference
    return ratio

# newton raphson solver
def newton_raphson(x0,constants,e1,step,eqs,deriv,maxiter=1000):

    def dx(x,eqs):
        f_ = eqs(x,constants)
        result =(abs(0-f_))
        return result
    
    def nr(x0,eqs,deriv):
        f_ = eqs(x0,constants)
        df_ = deriv
        x0 = x0 - step*(f_/df_)
        return x0
    
    # create results table
    results = pd.DataFrame([["guessx","diff","step"]])  
    results.to_csv('results_nr_isotopes.csv', index=False, header=False)
    diff = dx(x0,eqs)
    results1 = pd.DataFrame([[x0,diff,step]]) 
    results = results.append(results1, ignore_index=True)
    results.to_csv('results_nr_isotopes.csv', index=False, header=False) 
    
    delta1 = dx(x0,eqs)
    results = pd.DataFrame([["guessx","diff","step"]])  
    results.to_csv('results_nr_isotopes.csv', index=False, header=False)
    results1 = pd.DataFrame([[x0,delta1,step]]) 
    results = results.append(results1, ignore_index=True)
    results.to_csv('results_nr_isotopes.csv', index=False, header=False) 
    
    while delta1 > e1:
        guessx = x0
        f_ = eqs(x0,constants)
        df_ = deriv(x0,constants)
        x0 = x0 - step*(f_/df_)
        #while x0 < 0.:
        #    step = step/10.
        #    x0 = x0 - step*(f_/df_)
        delta1 = dx(x0,eqs)      
        results1 = pd.DataFrame([[x0,delta1,step]])
        results = results.append(results1, ignore_index=True)
        results.to_csv('results_nr_isotopes.csv', index=False, header=False) 
    return x0  

# two isotopes, two species
def i2s2(element,PT,R_i,melt_wf):
    if element == "S":
        knowns = i2s2_S_melt(PT,R_i,melt_wf)
    a, x, L = knowns # a = fractionation factor A-B, x = mole fraction of S in B, L = mole fraction of light isotope total 
    A = a - 1.
    B = (1. - a)*(L + x) - 1.
    C = a*x*L
    l_b = (-B - ((B**2) - (4.*A*C))**0.5)/(2.*A)
    h_b = x - l_b
    R_b = h_b/l_b
    R_a = a*R_b
    return R_a, R_b

def i2s2_S_melt(PT,R_i,melt_wf):
    a = (1./mg.alpha_SO2_SO4(PT))*mg.alpha_H2S_S(PT)*mg.alpha_gas("S","SO2","H2S",PT) # SO4-S2-
    x = 1.-melt_wf["S6ST"] # mole fraction of S as S2- in the melt
    R_i_ = 1./R_i["S"] # 32S/34S
    L = R_i_/(R_i_ + 1.) # mole fraction of 32S
    return a, x, L

def i2s6_S_alphas(PT): # all alphas against S2- in the melt
    a_b = mg.alpha_H2S_S(PT) # H2S-S
    a_c = (1./mg.alpha_SO2_SO4(PT))*mg.alpha_H2S_S(PT)*mg.alpha_gas("S","SO2","H2S",PT) # SO4-S 
    a_d = mg.alpha_H2S_S(PT)*mg.alpha_gas("S","S2","H2S",PT) # S2-S
    a_e = mg.alpha_H2S_S(PT)*mg.alpha_gas("S","SO2","H2S",PT) # SO2-S
    a_f = mg.alpha_H2S_S(PT)*mg.alpha_gas("S","OCS","H2S",PT) # OCS-S
    return a_b, a_c, a_d, a_e, a_f

def i2s6(element,PT,R,melt_wf,gas_mf,nr_step,nr_tol,guessx): # species distribution is mole fraction of S in each species
    
    if element == "S":
        a_b, a_c, a_d, a_e, a_f = i2s6_S_alphas(PT)
        species_distribution = c.mf_S_species(melt_wf,gas_mf)
        T_a = species_distribution["S2-"]
        T_b = species_distribution["H2S"]
        T_c = species_distribution["SO42-"]
        T_d = species_distribution["S2"]
        T_e = species_distribution["SO2"]
        T_f = species_distribution["OCS"]
        R_i = R["S"]
        
    constants = a_b, a_c, a_d, a_e, a_f, T_a, T_b, T_c, T_d, T_e, T_f, R_i
   
    def isotope_distribution(l_a, constants):
        a_b, a_c, a_d, a_e, a_f, T_a, T_b, T_c, T_d, T_e, T_f, R_i = constants
        R_a = (T_a - l_a)/l_a # 34S/32S
        R_b = a_b*R_a
        R_c = a_c*R_a
        R_d = a_d*R_a
        R_e = a_e*R_a
        R_f = a_f*R_a
        return R_a, R_b, R_c, R_d, R_e, R_f
    
    def av_m_g(element,l_a,constants):
        
        a_b, a_c, a_d, a_e, a_f, T_a, T_b, T_c, T_d, T_e, T_f, R_i = constants
        R_a, R_b, R_c, R_d, R_e, R_f = isotope_distribution(l_a, constants)
        
        # heavy/total isotope ratio
        R_a_ = R_a/(R_a + 1.)
        R_b_ = R_b/(R_b + 1.)
        R_c_ = R_c/(R_c + 1.)
        R_d_ = R_d/(R_d + 1.)
        R_e_ = R_e/(R_e + 1.)
        R_f_ = R_f/(R_f + 1.)
        
        if element == "S":
            h_m = R_a_*T_a + R_c_*T_c # 34S melt (S2- + SO42-)
            l_m = (1. - R_a_)*T_a + (1. - R_c_)*T_c
            h_g = R_b_*T_b + R_d_*T_d + R_e_*T_e + R_f_*T_f # 34S gas (H2S + S2 + SO2 + OCS)
            l_g = (1. - R_b_)*T_b + (1. - R_d_)*T_d + (1. - R_e_)*T_e + (1. - R_f_)*T_f # 34S gas (H2S + S2 + SO2 + OCS)
   
        R_m = h_m/l_m
        R_g = h_g/l_g
        return R_m, R_g
            
    def f(l_a, constants):
        a_b, a_c, a_d, a_e, a_f, T_a, T_b, T_c, T_d, T_e, T_f, R_i = constants
        R_a = (T_a/l_a) - 1.
        l_b = T_b/(1.+a_b*R_a)
        l_c = T_c/(1.+a_c*R_a)
        l_d = T_d/(1.+a_d*R_a)
        l_e = T_e/(1.+a_e*R_a)
        l_f = T_f/(1.+a_f*R_a)
        total = l_a + l_b + l_c + l_d + l_e + l_f
        R_i_ = 1./R_i
        l_t = (R_i_/(R_i_ + 1.))
        return l_t - total
    
    def df(l_a, constants):
        a_b, a_c, a_d, a_e, a_f, T_a, T_b, T_c, T_d, T_e, T_f, R_i = constants
        R_i_ = 1./R_i
        l_t = (R_i_/(R_i_ + 1.))
        result = -T_a*T_b*a_b/(l_a**2*(a_b*(T_a/l_a - 1.0) + 1.0)**2) - T_a*T_c*a_c/(l_a**2*(a_c*(T_a/l_a - 1.0) + 1.0)**2) - T_a*T_d*a_d/(l_a**2*(a_d*(T_a/l_a - 1.0) + 1.0)**2) - T_a*T_f*a_f/(l_a**2*(a_f*(T_a/l_a - 1.0) + 1.0)**2) - 1.
        return result
    
    l_a = newton_raphson(guessx,constants,nr_tol,nr_step,f,df)
    result1 = isotope_distribution(l_a, constants)
    result2 = av_m_g(element,l_a,constants)
    return result1, result2