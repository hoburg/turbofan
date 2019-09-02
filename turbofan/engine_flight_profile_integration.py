"""Simple commercial aircraft flight profile and aircraft model"""
from __future__ import absolute_import
from builtins import range
from numpy import pi
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, Vectorize
from gpkit.constraints.sigeq import SignomialEquality as SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import Tight as TCS
import matplotlib.pyplot as plt
from gpkit.small_scripts import mag
from .simple_ac_imports import Aircraft, CruiseSegment, ClimbSegment, FlightState
from .get_parametric_studies_subs import get_parametric_studies_subs

"""
Models requird to minimze the aircraft total fuel weight. Rate of climb equation taken from John
Anderson's Aircraft Performance and Design (eqn 5.85).
Inputs
-----
- Number of passtengers
- Passegner weight [N]
- Fusealge area per passenger (recommended to use 1 m^2 based on research) [m^2]
- Engine weight [N]
- Number of engines
- Required mission range [nm]
- Oswald efficiency factor
- Max allowed wing span [m]
- Cruise altitude [ft]
"""

class StateLinking(Model):
    """
    link all the state model variables
    """
    def setup(self, climbstate, cruisestate, enginestate, Nclimb, Ncruise):
        statevarkeys = ['p_{sl}', 'T_{sl}', 'L_{atm}', 'M_{atm}', 'P_{atm}', 'R_{atm}',
                        '\\rho', 'T_{atm}', '\\mu', 'T_s', 'C_1', 'h', 'hft', 'V', 'a', 'R', '\\gamma', 'M']
        constraints = []
        for i in range(len(statevarkeys)):
            varkey = statevarkeys[i]
            for i in range(Nclimb):
                constraints.extend([
                    climbstate[varkey][i] == enginestate[varkey][i]
                    ])
            for i in range(Ncruise):
                constraints.extend([
                    cruisestate[varkey][i] == enginestate[varkey][i+Nclimb]
                    ])           
        
        return constraints

class Mission(Model):
    """
    mission class, links together all subclasses
    """
    def setup(self, Nclimb, Ncruise, substitutions = None, **kwargs):
        eng = 0
        
        # vectorize
        with Vectorize(Nclimb + Ncruise):
            enginestate = FlightState()

        ac = Aircraft(Nclimb, Ncruise, enginestate, eng)
        
        #Vectorize
        with Vectorize(Nclimb):
            climb = ClimbSegment(ac)

        with Vectorize(Ncruise):
            cruise = CruiseSegment(ac)

        statelinking = StateLinking(climb.state, cruise.state, enginestate, Nclimb, Ncruise)

        #declare new variables
        W_ftotal = Variable('W_{f_{total}}', 'N', 'Total Fuel Weight')
        W_fclimb = Variable('W_{f_{climb}}', 'N', 'Fuel Weight Burned in Climb')
        W_fcruise = Variable('W_{f_{cruise}}', 'N', 'Fuel Weight Burned in Cruise')
        W_total = Variable('W_{total}', 'N', 'Total Aircraft Weight')
        CruiseAlt = Variable('CruiseAlt', 'ft', 'Cruise Altitude [feet]')
        ReqRng = Variable('ReqRng', 'nautical_miles', 'Required Cruise Range')
        W_dry = Variable('W_{dry}', 'N', 'Aircraft Dry Weight')

        RCmin = Variable('RC_{min}', 'ft/min', 'Minimum allowed climb rate')

        dhfthold = Variable('dhfthold', 'ft', 'Hold Variable')

        h = climb['h']
        hftClimb = climb['hft']
        dhft = climb['dhft']
        hftCruise = cruise['hft']

        #make overall constraints
        constraints = []

        constraints.extend([
            #weight constraints
            TCS([ac['W_{e}'] + ac['W_{payload}'] + ac['numeng'] * ac['W_{engine}'] + ac['W_{wing}'] <= W_dry]),
            TCS([W_dry + W_ftotal <= W_total]),

            climb['W_{start}'][0] == W_total,
            climb['W_{end}'][-1] == cruise['W_{start}'][0],

            # similar constraint 1
            TCS([climb['W_{start}'] >= climb['W_{end}'] + climb['W_{burn}']]),
            # similar constraint 2
            TCS([cruise['W_{start}'] >= cruise['W_{end}'] + cruise['W_{burn}']]),

            climb['W_{start}'][1:] == climb['W_{end}'][:-1],
            cruise['W_{start}'][1:] == cruise['W_{end}'][:-1],

            TCS([W_dry <= cruise['W_{end}'][-1]]),

            TCS([W_ftotal >=  W_fclimb + W_fcruise]),
            TCS([W_fclimb >= sum(climb['W_{burn}'])]),
            TCS([W_fcruise >= sum(cruise['W_{burn}'])]),

            #altitude constraints
            hftCruise == CruiseAlt,
            TCS([hftClimb[1:Nclimb] >= hftClimb[:Nclimb-1] + dhft[:Nclimb-1]]),
            TCS([hftClimb[0] >= dhft[0]]),
            hftClimb[-1] <= hftCruise,

            #compute the dh
            dhfthold == hftCruise[0]/Nclimb,

            dhft == dhfthold,

            #set the range for each cruise segment, doesn't take credit for climb
            #down range disatnce covered
            cruise.cruiseP['Rng'] == ReqRng/(Ncruise),
            
            #compute fuel burn from TSFC
            cruise['W_{burn}'] == ac['numeng']*ac.engine['TSFC'][Nclimb:] * cruise['thr'] * ac.engine['F'][Nclimb:],              
            climb['W_{burn}'] == ac['numeng']*ac.engine['TSFC'][:Nclimb] * climb['thr'] * ac.engine['F'][:Nclimb],

            CruiseAlt >= 30000*units('ft'),
            
            #min climb rate constraint
            climb['RC'] >= RCmin,
            ])

        M2 = .8
        M25 = .6
        M4a = .1025
        M0 = .8

        engineclimb = [
            ac.engine.engineP['M_2'][:Nclimb] == climb['M'],
            ac.engine.engineP['M_{2.5}'][:Nclimb] == M25,
            ac.engine.engineP['hold_{2}'] == 1+.5*(1.398-1)*M2**2,
            ac.engine.engineP['hold_{2.5}'] == 1+.5*(1.354-1)*M25**2,
            ac.engine.engineP['c1'] == 1+.5*(.401)*M0**2,

            #constraint on drag and thrust
            ac['numeng']*ac.engine['F'][:Nclimb] >= climb['D'] + climb['W_{avg}'] * climb['\\theta'],

            #climb rate constraints
            TCS([climb['excessP'] + climb.state['V'] * climb['D'] <=  climb.state['V'] * ac['numeng'] * ac.engine['F'][:Nclimb]]),
            ]

        M25 = .6

        enginecruise = [
            ac.engine.engineP['M_2'][Nclimb:] == cruise['M'],
            ac.engine.engineP['M_{2.5}'][Nclimb:] == M25,
            
            #steady level flight constraint on D 
            cruise['D'] == ac['numeng'] * ac.engine['F'][Nclimb:],

            #breguet range eqn
            TCS([cruise['z_{bre}'] >= (ac.engine['TSFC'][Nclimb:] * cruise['thr']*
            cruise['D']) / cruise['W_{avg}']]),
            ]
        
        return constraints + ac + climb + cruise + enginecruise + engineclimb + enginestate + statelinking

def test():
    substitutions = get_parametric_studies_subs()    
    substitutions.update({'ReqRng': 2000}) #('sweep', np.linspace(500,2000,4)),
           
    mission = Mission(2, 2)
    m = Model(mission['W_{f_{total}}'], mission, substitutions)
    sol = m.localsolve(verbosity = 0)

if __name__ == '__main__':
    plotRC = True
    plotR = False
    plotAlt = False

    fs = 18
    
    substitutions = get_parametric_studies_subs()    
    substitutions.update({'ReqRng': 2000}) #('sweep', np.linspace(500,2000,4)),
           
    mission = Mission(2, 2)
    m = Model(mission['W_{f_{total}}'], mission, substitutions)
    sol = m.localsolve(solver='mosek', verbosity = 1)

    if plotR == True:
        substitutions = get_parametric_studies_subs()    
        substitutions.update({'ReqRng': ('sweep', np.linspace(1000,3000,23)),})
          
        mission = Mission(2, 2)
        m = Model(mission['W_{f_{total}}'], mission, substitutions, x0=x0)
        solRsweep = m.localsolve(solver='mosek', verbosity = 1, skipsweepfailures=True)

        plt.plot(solRsweep('ReqRng'), solRsweep('W_{f_{total}}'), '-r', linewidth=2.0)
        plt.xlabel('Mission Range [nm]')
        plt.ylabel('Total Fuel Burn [N]')
        plt.title('Fuel Burn vs Range')
        plt.savefig('engine_Rsweeps/fuel_burn_range.pdf')
        plt.show()

        plt.plot(solRsweep('ReqRng'), solRsweep('b'), '-r', linewidth=2.0)
        plt.xlabel('Mission Range [nm]')
        plt.ylabel('Wing Span [m]')
        plt.title('Wing Span vs Range')
        plt.savefig('engine_Rsweeps/b_range.pdf')
        plt.show()

        plt.plot(solRsweep('ReqRng'), solRsweep('AR'), '-r', linewidth=2.0)
        plt.xlabel('Mission Range [nm]')
        plt.ylabel('Wing Aspect Ratio')
        plt.title('Wing Aspect Ratio vs Range')
        plt.savefig('engine_Rsweeps/AR_range.pdf')
        plt.show()

        plt.plot(solRsweep('ReqRng'), solRsweep('CruiseAlt'), '-r', linewidth=2.0)
        plt.xlabel('Mission Range [nm]')
        plt.ylabel('Cruise Altitude [ft]')
        plt.title('Cruise Altitude vs Range')
        plt.savefig('engine_Rsweeps/cruise_altitude_range.pdf')
        plt.show()

        plt.plot(solRsweep('ReqRng'), solRsweep('A_{2}'), '-r', linewidth=2.0)
        plt.xlabel('Mission Range [nm]')
        plt.ylabel('Fan Area [m^$2$]')
        plt.title('Fan Area vs Range')
        plt.savefig('engine_Rsweeps/fan_area_R.pdf')
        plt.show()
        
        irc = []
        f = []
        f6 = []
        f8 = []
        totsfc = []
        cruisetsfc = []
        maxm = []
        maxF = []
        cruiseF = []
        
        i=0
        while i < len(solRsweep('RC')):
            irc.append(mag(solRsweep('RC')[i][0]))
            f.append(mag(solRsweep('F')[i][0]))
            f6.append(mag(solRsweep('F_{6}')[i][0]))
            f8.append(mag(solRsweep('F_{8}')[i][0]))
            totsfc.append(mag(solRsweep('TSFC')[i][0]))
            cruisetsfc.append(mag(solRsweep('TSFC')[i][2]))
            maxm.append(max(mag(solRsweep('m_{total}')[i])))
            maxF.append(max(mag(solRsweep('F')[i])))
            cruiseF.append(mag(solRsweep('F')[i][2]))
            i+=1

        plt.plot(solRsweep('ReqRng'), cruiseF, '-r', linewidth=2.0)
        plt.xlabel('Mission Range [nm]')
        plt.ylabel('Initial Cruise Thrust [N]')
        plt.title('Initial Cruise Thrust vs Range')
        plt.savefig('engine_Rsweeps/max_m_range.pdf')
        plt.show()

        plt.plot(solRsweep('ReqRng'), maxm, '-r', linewidth=2.0)
        plt.xlabel('Mission Range [nm]')
        plt.ylabel('Max Engine Mass Flow [kg/s]')
        plt.title('Max Engine Mass Flow vs Range')
        
        plt.savefig('engine_Rsweeps/max_m_range.pdf', bbox_inches="tight")
        plt.show()

        plt.plot(solRsweep('ReqRng'), maxF, '-r', linewidth=2.0)
        plt.xlabel('Mission Range [nm]', fontsize = fs)
        plt.ylabel('Max Engine Thrust [N]', fontsize = fs)
        plt.title('Max Engine Thrust vs Range', fontsize = fs)
        plt.tick_params(axis='both', which='major', labelsize=fs)
        plt.ylim((10000,25000))
        plt.savefig('engine_Rsweeps/max_F_range.pdf', bbox_inches="tight")
        plt.show()

        plt.plot(solRsweep('ReqRng'), totsfc, '-r', linewidth=2.0)
        plt.plot(solRsweep('ReqRng'), cruisetsfc, '-b', linewidth=2.0)
        plt.legend(['Initial Climb', 'Initial Cruise'], loc=4, fontsize = fs)
        plt.ylim((0,.5))
        plt.xlabel('Mission Range [nm]', fontsize = fs)
        plt.ylabel('TSFC [1/hr]', fontsize = fs)
        plt.title('TSFC vs Range', fontsize = fs)
        plt.tick_params(axis='both', which='major', labelsize=fs)
        plt.savefig('engine_Rsweeps/TSFC_range.pdf', bbox_inches="tight")
        plt.show()

        plt.plot(solRsweep('ReqRng'), irc, '-r', linewidth=2.0)
        plt.xlabel('Mission Range [nm]')
        plt.ylabel('Initial Rate of Climb [ft/min]')
        plt.title('Initial Rate of Climb vs Range')
        plt.savefig('engine_Rsweeps/initial_RC_range.pdf')
        plt.show()

        plt.plot(solRsweep('ReqRng'), f, '-r', linewidth=2.0)
        plt.xlabel('Mission Range [nm]', fontsize = fs)
        plt.ylabel('Initial Thrsut [N]', fontsize = fs)
        plt.title('Initial Thrust vs Range', fontsize = fs)
        plt.title('Max Engine Mass Flow vs Range', fontsize = fs)
        plt.savefig('engine_Rsweeps/intitial_thrust.pdf', bbox_inches="tight")
        plt.show()

        plt.plot(solRsweep('ReqRng'), f6, '-r', linewidth=2.0)
        plt.xlabel('Mission Range [nm]')
        plt.ylabel('Initial Core Thrsut [N]')
        plt.title('Initial Core Thrust vs Range')
        plt.savefig('engine_Rsweeps/initial_F6_range.pdf')
        plt.show()

        plt.plot(solRsweep('ReqRng'), f8, '-r', linewidth=2.0)
        plt.xlabel('Mission Range [nm]')
        plt.ylabel('Initial Fan Thrsut [N]')
        plt.title('Initial Fan Thrust vs Range')
        plt.savefig('engine_Rsweeps/initial_F8_range.pdf')
        plt.show()

        plt.plot(solRsweep('ReqRng'), f8, '-r', linewidth=2.0)
        plt.plot(solRsweep('ReqRng'), f6, '-b', linewidth=2.0)
        plt.legend(['Initial Fan Thrust', 'Initial Core Thrust'], loc=2, fontsize = fs)
        plt.xlabel('Mission Range [nm]', fontsize = fs)
        plt.ylabel('Initial Thrust [N]', fontsize = fs)
        plt.title('Initial Thrust vs Range', fontsize = fs)
        plt.ylim((0,20000))
        plt.tick_params(axis='both', which='major', labelsize=fs)
        plt.savefig('engine_Rsweeps/initial_F8_F6_range.pdf', bbox_inches="tight")
        plt.show()

        plt.plot(solRsweep('ReqRng'), solRsweep('W_{engine}'), '-r', linewidth=2.0)
        plt.xlabel('Mission Range [nm]', fontsize = fs)
        plt.ylabel('Engine Weight [N]', fontsize = fs)
        plt.title('Engine Weight vs Range', fontsize = fs)
        plt.ylim((0,15000))
        plt.tick_params(axis='both', which='major', labelsize=fs)
        plt.savefig('engine_Rsweeps/engine_weight_range.pdf', bbox_inches="tight")
        plt.show()

        plt.plot(solRsweep('ReqRng'), solRsweep('A_{2}'), '-r', linewidth=2.0)
        plt.xlabel('Mission Range [nm]')
        plt.ylabel('Fan Area [m^$2$]')
        plt.title('Fan Area vs Range')
        plt.ylim((.5,.7))
        plt.savefig('engine_Rsweeps/fan_area_range.pdf')
        plt.show()

        plt.plot(solRsweep('ReqRng'), solRsweep('A_5'), '-r', linewidth=2.0)
        plt.xlabel('Mission Range [nm]')
        plt.ylabel('$A_5$ [m^$2$]')
        plt.title('$A_5$ vs Range')
        plt.savefig('engine_Rsweeps/a5_range.pdf')
        plt.show()

        plt.plot(solRsweep('ReqRng'), solRsweep('A_{2.5}'), '-r', linewidth=2.0)
        plt.xlabel('Mission Range [nm]')
        plt.ylabel('$A_{2.5}$ [m^$2$]')
        plt.title('$A_{2.5}$ vs Range')
        plt.savefig('engine_Rsweeps/a25_range.pdf')
        plt.show()

        plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['M_{takeoff}'], '-r', linewidth=2.0)
        plt.ylabel('Sensitivity to $M_{takeoff}$')
        plt.xlabel('Mission Range [nm]')
        plt.title('Sensitivity to $M_{takeoff}$ vs Range')
        plt.savefig('engine_Rsweeps/mtakeoff_sens_range.pdf')
        plt.show()

        plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['\\\pi_{f_D}'], '-r', linewidth=2.0)
        plt.ylabel('Sensitivity to $\\pi_{f_D}$')
        plt.xlabel('Mission Range [nm]')
        plt.title('Sensitivity to $\\pi_{f_D}$ vs Range')
        plt.savefig('engine_Rsweeps/pifd_sens_range.pdf')
        plt.show()

        plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['\\\pi_{lc_D}'], '-r', linewidth=2.0)
        plt.ylabel('Sensitivity to $\\pi_{lc_D}$')
        plt.xlabel('Mission Range [nm]')
        plt.title('Sensitivity to $\\pi_{lc_D}$ vs Range')
        plt.savefig('engine_Rsweeps/pilcD_sens_range.pdf')
        plt.show()

        plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['\\\pi_{hc_D}'], '-r', linewidth=2.0)
        plt.ylabel('Sensitivity to $\\pi_{hc_D}$')
        plt.xlabel('Mission Range [nm]')
        plt.title('Sensitivity to $\\pi_{hc_D}$ vs Range')
        plt.savefig('engine_Rsweeps/pihcD_sens_range.pdf')
        plt.show()

        plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['T_{t_f}'], '-r', linewidth=2.0)
        plt.ylabel('Sensitivity to $T_{t_f}$')
        plt.xlabel('Mission Range [nm]')
        plt.title('Sensitivity to $T_{t_f}$ vs Range')
        plt.savefig('engine_Rsweeps/ttf_sens_range.pdf')
        plt.show()

        plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['\\alpha_c'], '-r', linewidth=2.0)
        plt.ylabel('Sensitivity to $\\alpha_c$')
        plt.xlabel('Mission Range [nm]')
        plt.title('Sensitivity to $\\alpha_c$ vs Range')
        plt.savefig('engine_Rsweeps/alphac_sens_range.pdf')
        plt.show()

        plt.plot(solRsweep('ReqRng'), solRsweep('\\\pi_{f_D}'), '-r', linewidth=2.0)
        plt.ylabel('$\\pi_{f_D}$')
        plt.xlabel('Mission Range [nm]')
        plt.title('Fan Design Pressure Ratio vs Initial Rate of Climb')
        plt.savefig('engine_Rsweeps/pifD_R.pdf')
        plt.show()

    if plotAlt == True:
        substitutions = get_parametric_studies_subs()    
        substitutions.update({'ReqRng': 2000, 'CruiseAlt': ('sweep', np.linspace(30000,40000,20))}) #('sweep', np.linspace(500,2000,4)),
    
        mmission = Mission(2, 2)
        m = Model(mission['W_{f_{total}}'], mission, substitutions)
        solAltsweep = m.localsolve(solver='mosek', verbosity = 4, skipsweepfailures=True)

        irc = []
        f = []
        f6 = []
        f8 = []
        i=0
        while i < len(solAltsweep('RC')):
            irc.append(mag(solAltsweep('RC')[i][0]))
            f.append(mag(solAltsweep('F')[i][0]))
            f6.append(mag(solAltsweep('F_{6}')[i][0]))
            f8.append(mag(solAltsweep('F_{8}')[i][0]))
            i+=1

        plt.plot(solAltsweep('CruiseAlt'), irc, '-r')
        plt.xlabel('Mission Range [nm]')
        plt.ylabel('Initial Rate of Climb [ft/min]')
        plt.title('Initial Rate of Climb vs Cruise Altitude')
        plt.savefig('engine_Altsweeps/initial_RC_alt.pdf')
        plt.show()

        plt.plot(solAltsweep('CruiseAlt'), f, '-r')
        plt.xlabel('Mission Range [nm]')
        plt.ylabel('Initial Thrsut [N]')
        plt.title('Initial Thrust vs Cruise Altitude')
        plt.savefig('engine_Altsweeps/intitial_thrust_alt.pdf')
        plt.show()

        plt.plot(solAltsweep('CruiseAlt'), f6, '-r')
        plt.xlabel('Mission Range [nm]')
        plt.ylabel('Initial Core Thrsut [N]')
        plt.title('Initial Core Thrust vs Cruise Altitude')
        plt.savefig('engine_Altsweeps/initial_F6_alt.pdf')
        plt.show()

        plt.plot(solAltsweep('CruiseAlt'), f8, '-r')
        plt.xlabel('Mission Range [nm]')
        plt.ylabel('Initial Fan Thrsut [N]')
        plt.title('Initial Fan Thrust vs Cruise Altitude')
        plt.savefig('engine_Altsweeps/initial_F8_alt.pdf')
        plt.show()

        plt.plot(solAltsweep('CruiseAlt'), solAltsweep('W_{f_{total}}'), '-r')
        plt.xlabel('Cruise Alt [ft]')
        plt.ylabel('Total Fuel Burn [N]')
        plt.title('Fuel Burn vs Cruise Altitude')
        plt.savefig('engine_Altsweeps/fuel_alt.pdf')
        plt.show()

        plt.plot(solAltsweep('CruiseAlt'), solAltsweep('W_{engine}'), '-r')
        plt.xlabel('Cruise Alt [ft]')
        plt.ylabel('Engine Weight [N]')
        plt.title('Engine WEight vs Cruise Altitude')
        plt.savefig('engine_Altsweeps/weight_engine_alt.pdf')
        plt.show()

        plt.plot(solAltsweep('CruiseAlt'), solAltsweep('A_{2}'), '-r')
        plt.xlabel('Cruise Alt [ft]')
        plt.ylabel('Fan Area [m^$2$]')
        plt.title('Fan Area vs Cruise Altitude')
        plt.savefig('engine_Altsweeps/fan_area_alt.pdf')
        plt.show()

        plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['M_{takeoff}'], '-r')
        plt.ylabel('Sensitivity to $M_{takeoff}$')
        plt.xlabel('Cruise Alt [ft]')
        plt.title('Fan Area vs Cruise Altitdue')
        plt.savefig('engine_Altsweeps/m_takeoff_sens_alt.pdf')
        plt.show()

        plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['\\pi_{f_D}'], '-r')
        plt.ylabel('Sensitivity to $\\pi_{f_D}$')
        plt.ylabel('Fan Area [m^$2$]')
        plt.title('Fan Area vs Cruise Altitude')
        plt.savefig('engine_Altsweeps/pifD_sens_alt.pdf')
        plt.show()

        plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['\\pi_{lc_D}'], '-r')
        plt.ylabel('Sensitivity to $\\pi_{lc_D}$')
        plt.xlabel('Cruise Alt [ft]')
        plt.title('Fan Area vs Cruise Altitdue')
        plt.savefig('engine_Altsweeps/pilcD_sens_alt.pdf')
        plt.show()

        plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['\\pi_{hc_D}'], '-r')
        plt.ylabel('Sensitivity to $\\pi_{hc_D}$')
        plt.xlabel('Cruise Alt [ft]')
        plt.title('Fan Area vs Cruise Altitude')
        plt.savefig('engine_Altsweeps/pihcD_sens_alt.pdf')
        plt.show()

        plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['T_{t_f}'], '-r')
        plt.ylabel('Sensitivity to $T_{t_f}$')
        plt.xlabel('Cruise Alt [ft]')
        plt.title('Fan Area vs Cruise Altitude')
        plt.savefig('engine_Altsweeps/Ttf_sens_alt.pdf')
        plt.show()

        plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['\\alpha_c'], '-r')
        plt.ylabel('Sensitivity to $\\alpha_c$')
        plt.xlabel('Cruise Alt [ft]')
        plt.title('Fan Area vs Cruise Altitude')
        plt.savefig('engine_Altsweeps/alpha_c_sens_alt.pdf')
        plt.show()

    if plotRC == True:
        substitutions = get_parametric_studies_subs()    
        substitutions.update({'ReqRng': 2000, 'RC_{min}': ('sweep', np.linspace(500,3500,10))}) #('sweep', np.linspace(500,2000,4)),

        mission = Mission(2, 2)
        m = Model(mission['W_{f_{total}}'], mission, substitutions)
        solRCsweep = m.localsolve(solver='mosek', verbosity = 2, skipsweepfailures=True)

        i = 0

        f = []
        f6 = []
        f8 = []
        crtsfc = []
        itsfc = []
        maxm = []
        maxF = []
        cruiseF = []

        while i < len(solRCsweep('RC')):
            f.append(mag(solRCsweep('F')[i][0]))
            f6.append(mag(solRCsweep('F_{6}')[i][0]))
            f8.append(mag(solRCsweep('F_{8}')[i][0]))
            crtsfc.append(mag(solRCsweep('TSFC')[i][2]))
            itsfc.append(mag(solRCsweep('TSFC')[i][0]))
            maxm.append(max(mag(solRCsweep('m_{total}')[i])))
            maxF.append(max(mag(solRCsweep('F')[i])))
            cruiseF.append(mag(solRCsweep('F')[i][2]))
            i+=1

        plt.plot(solRCsweep('RC_{min}'), cruiseF, '-r', linewidth=2.0)
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
        plt.ylabel('Initial Cruise Thrust [N]')
        plt.title('Initial Cruise Thrust vs Range')
        plt.xlim((500, 3500))
        plt.savefig('engine_RCsweeps/max_m_range.pdf')
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), maxm, '-r', linewidth=2.0)
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
        plt.ylabel('Max Engine Mass Flow [kg/s]')
        plt.title('Max Engine Mass Flow vs Range')
        plt.ylim((0,.180))
        plt.xlim((500, 3500))
        plt.savefig('engine_RCsweeps/max_m_range.pdf')
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), maxF, '-r', linewidth=2.0)
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]', fontsize = fs)
        plt.ylabel('Max Engine Thrust [N]',fontsize = fs)
        plt.title('Max Engine Thrust vs Minimum Initial Climb Rate',fontsize = fs)
        plt.tick_params(axis='both', which='major', labelsize=fs)
        plt.xlim((500, 3500))
        plt.ylim((0,55000))
        plt.savefig('engine_RCsweeps/max_F_RC.pdf', bbox_inches="tight")
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), solRCsweep('CruiseAlt'), '-r', linewidth=2.0)
        plt.ylabel('Cruise Altitude [ft]')
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
        plt.title('Cruise Altitude vs Initial Rate of Climb')
        plt.xlim((500, 3500))
        plt.savefig('engine_RCsweeps/cralt_RC.pdf')
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), itsfc, '-r', linewidth=2.0)
        plt.ylabel('Initial Climb TSFC [1/hr]')
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
        plt.title('Initial Climb TSFC vs Initial Rate of Climb')
        plt.xlim((500, 3500))
        plt.savefig('engine_RCsweeps/itsfc_RC.pdf')
        plt.ylim((0,.7))
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), itsfc, '-r', linewidth=2.0)
        plt.plot(solRCsweep('RC_{min}'), crtsfc, '-b', linewidth=2.0)
        plt.legend(['Initial Climb', 'Initial Cruise'], loc=2)
        plt.ylabel('TSFC [1/hr]')
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
        plt.title('Initial Climb and Cruise TSFC vs Minimum Initial Rate of Climb')
        plt.xlim((500, 3500))
        plt.savefig('engine_RCsweeps/tsfc_RC.pdf')
        plt.ylim((0,.7))
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), crtsfc, '-r', linewidth=2.0)
        plt.ylabel('Initial Cruise TSFC [1/hr]')
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
        plt.title('Initial Cruise TSFC vs Minimum Initial Rate of Climb')
        plt.xlim((500, 3500))
        plt.ylim((0,.7))
        plt.savefig('engine_RCsweeps/crtsfc_RC.pdf')
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), f, '-r', linewidth=2.0)
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
        plt.ylabel('Initial Thrsut [N]')
        plt.title('Initial Thrust vs Initial Rate of Climb')
        plt.xlim((500, 3500))
        plt.savefig('engine_RCsweeps/intitial_thrust_RC.pdf')
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), f6, '-r', linewidth=2.0)
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
        plt.ylabel('Initial Core Thrsut [N]')
        plt.title('Initial Core Thrust vs Initial Rate of Climb')
        plt.xlim((500, 3500))
        plt.savefig('engine_RCsweeps/initial_F6_RC.pdf')
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), f8, '-r', linewidth=2.0)
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
        plt.ylabel('Initial Fan Thrsut [N]')
        plt.title('Initial Fan Thrust vs Initial Rate of Climb')
        plt.xlim((500, 3500))
        plt.savefig('engine_RCsweeps/initial_F8_RC.pdf')
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), f8, '-r', linewidth=2.0)
        plt.plot(solRCsweep('RC_{min}'), f6, '-b', linewidth=2.0)
        plt.legend(['Initial Fan Thrust', 'Initial Core Thrust'], loc=2,fontsize = fs)
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]',fontsize = fs)
        plt.ylabel('Initial Thrsut [N]',fontsize = fs)
        plt.title('Initial Fan & Core Thrust vs Min. Initial Climb Rate',fontsize = fs)
        plt.tick_params(axis='both', which='major', labelsize=fs)
        plt.xlim((500, 3500))
        plt.ylim((0,30000))
        plt.savefig('engine_RCsweeps/initial_F8_F6_RC.pdf', bbox_inches="tight")
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), solRCsweep('W_{f_{total}}'), '-r', linewidth=2.0)
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
        plt.ylabel('Total Fuel Burn [N]')
        plt.title('Fuel Burn vs Initial Rate of Climb')
        plt.xlim((500, 3500))
        plt.ylim((30000,90000))
        plt.savefig('engine_RCsweeps/fuel_RC.pdf')
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), solRCsweep('W_{engine}'), '-r', linewidth=2.0)
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
        plt.ylabel('Engine Weight [N]')
        plt.title('Engine Weight vs Minimum Initial Rate of Climb')
        plt.xlim((500, 3500))
        plt.ylim((0, 19000))
        plt.savefig('engine_RCsweeps/weight_engine_RC.pdf')
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), solRCsweep('A_{2}'), '-r', linewidth=2.0)
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
        plt.ylabel('Fan Area [m$^2$]')
        plt.title('Fan Area vs  Minimum Initial Rate of Climb')
        plt.ylim((0, 0.8))
        plt.xlim((500, 3500))
        plt.savefig('engine_RCsweeps/fan_area_RC.pdf')
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), solRCsweep('\\pi_{f_D}'), '-r', linewidth=2.0)
        plt.ylabel('$\\pi_{f_D}$')
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
        plt.title('Fan Design Pressure Ratio vs Min. Initial Rate of Climb')
        plt.xlim((500, 3500))
        plt.savefig('engine_RCsweeps/pifD_RC.pdf')
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), solRCsweep['sensitivities']['constants']['M_{takeoff}'], '-r', linewidth=2.0)
        plt.ylabel('Sensitivity to $M_{takeoff}$')
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
        plt.title('Core Mass Flow Bleed vs Initial Rate of Climb')
        plt.xlim((500, 3500))
        plt.savefig('engine_RCsweeps/m_takeoff_sens_RC.pdf')
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), solRCsweep['sensitivities']['constants']['\\pi_{f_D}'], '-r', linewidth=2.0)
        plt.ylabel('Sensitivity to $\\pi_{f_D}$')
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
        plt.ylim((0, .6))
        plt.title('Fan Design Pressure Ratio Sensitivity vs Min. Initial Climb Rate')
        plt.xlim((500, 3500))
        plt.savefig('engine_RCsweeps/pifD_sens_RC.pdf')
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), solRCsweep['sensitivities']['constants']['\\pi_{lc_D}'], '-r', linewidth=2.0)
        plt.ylabel('Sensitivity to $\\pi_{lc_D}$')
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
        plt.title('LPC Design Pressure Ratio Sensitivity vs Initial Rate of Climb')
        plt.xlim((500, 3500))
        plt.savefig('engine_RCsweeps/pilcD_sens_RC.pdf')
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), solRCsweep['sensitivities']['constants']['\\pi_{hc_D}'], '-r', linewidth=2.0)
        plt.ylabel('Sensitivity to $\\pi_{hc_D}$')
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
        plt.title('HPC Design Pressure Ratio Sensitivity vs Initial Rate of Climb')
        plt.xlim((500, 3500))
        plt.savefig('engine_RCsweeps/pihcD_sens_RC.pdf')
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), solRCsweep['sensitivities']['constants']['T_{t_f}'], '-r', linewidth=2.0)
        plt.ylabel('Sensitivity to $T_{t_f}$')
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
        plt.title('Input Fuel Temp Sensitivity vs Initial Rate of Climb')
        plt.xlim((500, 3500))
        plt.savefig('engine_RCsweeps/Ttf_sens_alt.pdf')
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), solRCsweep['sensitivities']['constants']['\\alpha_c'], '-r')
        plt.ylabel('Sensitivity to $\\alpha_c$')
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
        plt.title('Cooling Flow BPR Sensitivity vs Initial Rate of Climb')
        plt.xlim((500, 3500))
        plt.savefig('engine_RCsweeps/alpha_c_sens_alt.pdf')
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), solRCsweep('b'), '-r')
        plt.ylabel('Wing Span [m]')
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
        plt.title('Wing Span vs Initial Rate of Climb')
        plt.xlim((500, 3500))
        plt.savefig('engine_RCsweeps/b_RC.pdf')
        plt.show()

        plt.plot(solRCsweep('RC_{min}'), solRCsweep('AR'), '-r')
        plt.ylabel('Wing Aspect Ratio')
        plt.xlabel('Minimum Initial Rate of Climb [ft/min]')
        plt.title('Wing Aspect Ratio vs Initial Rate of Climb')
        plt.xlim((500, 3500))
        plt.savefig('engine_RCsweeps/AR_RC.pdf')
        plt.show()
