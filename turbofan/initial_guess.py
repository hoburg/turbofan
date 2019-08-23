from gpkit import units

def initialize_guess():
    return {'W_{engine}': 1e4*units('N'),
            'P_{t_0}': 1e1*units('kPa'),
            'T_{t_0}': 1e3*units('K'),
            'h_{t_0}': 1e6*units('J/kg'),
            'P_{t_{1.8}}': 1e1*units('kPa'),
            'T_{t_{1.8}}': 1e3*units('K'),
            'h_{t_{1.8}}': 1e6*units('J/kg'),
            'P_{T_{2}}': 1e1*units('kPa'),
            'T_{T_{2}}': 1e3*units('K'),
            'h_{T_{2}}': 1e6*units('J/kg'),
            'P_{t_{2.1}}': 1e3*units('K'),
            'T_{t_{2.1}}': 1e3*units('K'),
            'h_{t_{2.1}}': 1e6*units('J/kg'),
            'P_{t_{2.5}}': 1e3*units('kPa'),
            'T_{t_{2.5}}': 1e3*units('K'),
            'h_{t_{2.5}}': 1e6*units('J/kg'),
            'P_{t_3}': 1e4*units('kPa'),
            'T_{t_3}': 1e4*units('K'),
            'h_{t_3}': 1e7*units('J/kg'),
            'P_{t_7}': 1e2*units('kPa'),
            'T_{t_7}': 1e3*units('K'),
            'h_{t_7}': 1e6*units('J/kg'),
            'P_{t_4}': 1e4*units('kPa'),
            'h_{t_4}': 1e7*units('J/kg'),
            'T_{t_4}': 1e4*units('K'),
            'P_{t_{4.1}}': 1e4*units('kPa'),
            'T_{t_{4.1}}': 1e4*units('K'),
            'h_{t_{4.1}}': 1e7*units('J/kg'),
            'T_{4.1}': 1e4*units('K'),
            'f': 1e-2,
            'P_{4a}': 1e4*units('kPa'),
            'h_{t_{4.5}}': 1e6*units('J/kg'),
            'P_{t_{4.5}}': 1e3*units('kPa'),
            'T_{t_{4.5}}': 1e4*units('K'),
            'P_{t_{4.9}}': 1e2*units('kPa'),
            'T_{t_{4.9}}': 1e3*units('K'),
            'h_{t_{4.9}}': 1e6*units('J/kg'),
            '\\pi_{HPT}': 1e-1,
            '\\pi_{LPT}': 1e-1,
            'P_{t_5}': 1e2*units('kPa'),
            'T_{t_5}': 1e3*units('K'),
            'h_{t_5}': 1e6*units('J/kg'),
            'P_{8}': 1e2*units('kPa'),
            'P_{t_8}': 1e2*units('kPa'),
            'h_{t_8}': 1e6*units('J/kg'),
            'h_{8}': 1e6*units('J/kg'),
            'T_{t_8}': 1e3*units('K'),
            'T_{8}': 1e3*units('K'),
            'P_{6}': 1e2*units('kPa'),
            'P_{t_6}': 1e2*units('kPa'),
            'T_{t_6': 1e3*units('K'),
            'h_{t_6}': 1e6*units('J/kg'),
            'h_{6}': 1e6*units('J/kg'),
            'F_{8}': 1e2 * units('kN'),
            'F_{6}': 1e2 * units('kN'),
            'F': 1e2 * units('kN'),
            'F_{sp}': 1e-1,
            'TSFC': 1e-1,
            'I_{sp}': 1e4*units('s'),
            'u_{6}': 1e3*units('m/s'),
            'u_{8}': 1e3*units('m/s'),
            'm_{core}': 1e2*units('kg/s'),
            'm_{fan}': 1e3*units('kg/s'),
            '\\alpha': 1e1,
            '\\alpha_{+1}': 1e1,
            'm_{total}': 1e3*units('kg/s'),
            'T_{2}': 1e3*units('K'),
            'P_{2}': 1e2*units('kPa'),
            'u_{2}': 1e3*units('m/s'),
            'h_{2}': 1e6*units('J/kg'),
            'T_{2.5}': 1e3*units('K'),
            'P_{2.5}': 1e2*units('kPa'),
            'u_{2.5}': 1e3*units('m/s'),
            'h_{2.5}': 1e6*units('J/kg'),
            'P_{7}': 1e2*units('kPa'),
            'T_{7}': 1e3*units('K'),
            'u_{7}': 1e3*units('m/s'),
            'P_{5}': 1e2*units('kPa'),
            'T_{5}': 1e3*units('K'),
            'u_{5}': 1e3*units('m/s'),
            'P_{atm}': 1e2*units('kPa'),
            'T_{atm}': 1e3*units('K'),
            'V': 1e3*units('knot'),
            'a': 1e3*units('m/s'),
    }
