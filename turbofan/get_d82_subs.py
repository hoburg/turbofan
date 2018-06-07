"""
subs for TASOPT D8.2 engine
"""

def get_D82_subs():
    M4a = .1025
    fan = 1.60474
    lpc  = 4.98
    hpc = 35/8
    
    substitutions = {
        'OPR_{max}': 32,
        'T_{t_{4.1_{max}}}': 1400,

        '\\pi_{tn}': .995,
        '\\pi_{b}': .94,
        '\\pi_{d}': .995,
        '\\pi_{fn}': .985,
        'T_{ref}': 288.15,
        'P_{ref}': 101.325,
        '\\eta_{HPshaft}': .97,
        '\\eta_{LPshaft}': .97,
        '\\eta_{B}': .9827,

        '\\pi_{f_D}': fan,
        '\\pi_{hc_D}': hpc,
        '\\pi_{lc_D}': lpc,

        '\\alpha_{OD}': 6.97,
        '\\alpha_{max}': 6.97,

        'hold_{4a}': 1+.5*(1.313-1)*M4a**2,
        'r_{uc}': .01,
        '\\alpha_c': .19036,
        'T_{t_f}': 435,

        'M_{takeoff}': .9556,

        'G_{f}': 1,

        'h_{f}': 43.003,

        'C_{p_{t1}}': 1280,
        'C_{p_{t2}}': 1184,
        'C_{p_{c}}': 1216,

        'HTR_{f_{SUB}}': 1-.3**2,
        'HTR_{lpc_{SUB}}': 1 - 0.6**2,
     }

    return substitutions
           
