"""Simple commercial aircraft flight profile and aircraft model used in a multi-mission
optimization example"""
from __future__ import absolute_import
from builtins import range
from numpy import pi
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, Vectorize
from gpkit.constraints.sigeq import SignomialEquality as SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import Tight as TCS
import matplotlib.pyplot as plt
from simple_ac_imports import Aircraft, CruiseSegment, ClimbSegment, FlightState
from get_parametric_studies_subs import get_parametric_studies_subs
from engine_validation import Engine
Engine.Ttmax = False

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

class FleetMission(Model):
    """
    mission class, links together all subclasses
    """
    def setup(self, Nclimb, Ncruise, Nfleet):
        eng = 0

        #two level vectorization to make a fleet
        with Vectorize(Nfleet):
            # vectorize
            with Vectorize(Nclimb + Ncruise):
                enginestate = FlightState()

        ac = Aircraft(Nclimb, Ncruise, enginestate, eng, Nfleet)

        #two level vectorization to make a fleet
        with Vectorize(Nfleet):
            #Vectorize
            with Vectorize(Nclimb):
                climb = ClimbSegment(ac)

            with Vectorize(Ncruise):
                cruise = CruiseSegment(ac)

        statelinking = StateLinking(climb.state, cruise.state, enginestate, Nclimb, Ncruise)

        with Vectorize(Nfleet):
            #declare new variables
            W_ftotal = Variable('W_{f_{total}}', 'N', 'Total Fuel Weight')
            W_fclimb = Variable('W_{f_{climb}}', 'N', 'Fuel Weight Burned in Climb')
            W_fcruise = Variable('W_{f_{cruise}}', 'N', 'Fuel Weight Burned in Cruise')
            W_total = Variable('W_{total}', 'N', 'Total Aircraft Weight')
            CruiseAlt = Variable('CruiseAlt', 'ft', 'Cruise Altitude [feet]')
            ReqRng = Variable('ReqRng', 'nautical_miles', 'Required Cruise Range')
            W_dry = Variable('W_{dry}', 'N', 'Aircraft Dry Weight')

            dhfthold = Variable('dhfthold', 'ft', 'Hold Variable')

        W_ffleet = Variable('W_{f_{fleet}}', 'N', 'Total Fleet Fuel Burn')

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
            climb['RC'] >= 500*units('ft/min'),
            ])

        fleetfuel = [
            #compute the fleet fuel burn
            W_ffleet >= .375*W_ftotal[0] + .375*W_ftotal[1] + .125*W_ftotal[2] + .125*W_ftotal[3],
            ]

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

        ranges = [
            ReqRng[0] == 500*units('nautical_miles'),
            ReqRng[1] == 1000*units('nautical_miles'),
            ReqRng[2] == 1500*units('nautical_miles'),
            ReqRng[3] == 2000*units('nautical_miles'),
            ]

        return constraints + ac + climb + cruise + enginecruise + engineclimb + enginestate + statelinking + ranges + fleetfuel

def test():
    substitutions = get_parametric_studies_subs()

    #dict of initial guesses
    x0 = {
        'W_{engine}': 1e4*units('N'),
        'P_{t_0}': 1e1*units('kPa'),
        'T_{t_0}': 1e3*units('K'),
        'h_{t_0}': 1e6*units('J/kg'),
        'P_{t_{1.8}}': 1e1*units('kPa'),
        'T_{t_{1.8}}': 1e3*units('K'),
        'h_{t_{1.8}}': 1e6*units('J/kg'),
        'P_{t_2}': 1e1*units('kPa'),
        'T_{t_2}': 1e3*units('K'),
        'h_{t_2}': 1e6*units('J/kg'),
        'P_{t_2.1}': 1e3*units('K'),
        'T_{t_2.1}': 1e3*units('K'),
        'h_{t_2.1}': 1e6*units('J/kg'),
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
        '\pi_{HPT}': 1e-1,
        '\pi_{LPT}': 1e-1,
        'P_{t_5}': 1e2*units('kPa'),
        'T_{t_5}': 1e3*units('K'),
        'h_{t_5}': 1e6*units('J/kg'),
        'P_8': 1e2*units('kPa'),
        'P_{t_8}': 1e2*units('kPa'),
        'h_{t_8}': 1e6*units('J/kg'),
        'h_8': 1e6*units('J/kg'),
        'T_{t_8}': 1e3*units('K'),
        'T_{8}': 1e3*units('K'),
        'P_6': 1e2*units('kPa'),
        'P_{t_6}': 1e2*units('kPa'),
        'T_{t_6': 1e3*units('K'),
        'h_{t_6}': 1e6*units('J/kg'),
        'h_6': 1e6*units('J/kg'),
        'F_8': 1e2 * units('kN'),
        'F_6': 1e2 * units('kN'),
        'F': 1e2 * units('kN'),
        'F_{sp}': 1e-1,
        'TSFC': 1e-1,
        'I_{sp}': 1e4*units('s'),
        'u_6': 1e3*units('m/s'),
        'u_8': 1e3*units('m/s'),
        'm_{core}': 1e2*units('kg/s'),
        'm_{fan}': 1e3*units('kg/s'),
        '\\alpha': 1e1,
        'alphap1': 1e1,
        'm_{total}': 1e3*units('kg/s'),
        'T_2': 1e3*units('K'),
        'P_2': 1e2*units('kPa'),
        'u_2': 1e3*units('m/s'),
        'h_{2}': 1e6*units('J/kg'),
        'T_{2.5}': 1e3*units('K'),
        'P_{2.5}': 1e2*units('kPa'),
        'u_{2.5}': 1e3*units('m/s'),
        'h_{2.5}': 1e6*units('J/kg'),
        'P_{7}': 1e2*units('kPa'),
        'T_{7}': 1e3*units('K'),
        'u_7': 1e3*units('m/s'),
        'P_{5}': 1e2*units('kPa'),
        'T_{5}': 1e3*units('K'),
        'u_5': 1e3*units('m/s'),
        'P_{atm}': 1e2*units('kPa'),
        'T_{atm}': 1e3*units('K'),
        'V': 3e3*units('knot'),
        'a': 1e3*units('m/s'),
    }

    mission = FleetMission(2, 2, 4)
    m = Model(mission['W_{f_{fleet}}'], mission, substitutions)
    sol = m.localsolve(solver='mosek', verbosity = 4)

if __name__ == '__main__':
    substitutions = get_parametric_studies_subs()

    #dict of initial guesses
    x0 = {
        'W_{engine}': 1e4*units('N'),
        'P_{t_0}': 1e1*units('kPa'),
        'T_{t_0}': 1e3*units('K'),
        'h_{t_0}': 1e6*units('J/kg'),
        'P_{t_{1.8}}': 1e1*units('kPa'),
        'T_{t_{1.8}}': 1e3*units('K'),
        'h_{t_{1.8}}': 1e6*units('J/kg'),
        'P_{t_2}': 1e1*units('kPa'),
        'T_{t_2}': 1e3*units('K'),
        'h_{t_2}': 1e6*units('J/kg'),
        'P_{t_2.1}': 1e3*units('K'),
        'T_{t_2.1}': 1e3*units('K'),
        'h_{t_2.1}': 1e6*units('J/kg'),
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
        '\pi_{HPT}': 1e-1,
        '\pi_{LPT}': 1e-1,
        'P_{t_5}': 1e2*units('kPa'),
        'T_{t_5}': 1e3*units('K'),
        'h_{t_5}': 1e6*units('J/kg'),
        'P_8': 1e2*units('kPa'),
        'P_{t_8}': 1e2*units('kPa'),
        'h_{t_8}': 1e6*units('J/kg'),
        'h_8': 1e6*units('J/kg'),
        'T_{t_8}': 1e3*units('K'),
        'T_{8}': 1e3*units('K'),
        'P_6': 1e2*units('kPa'),
        'P_{t_6}': 1e2*units('kPa'),
        'T_{t_6': 1e3*units('K'),
        'h_{t_6}': 1e6*units('J/kg'),
        'h_6': 1e6*units('J/kg'),
        'F_8': 1e2 * units('kN'),
        'F_6': 1e2 * units('kN'),
        'F': 1e2 * units('kN'),
        'F_{sp}': 1e-1,
        'TSFC': 1e-1,
        'I_{sp}': 1e4*units('s'),
        'u_6': 1e3*units('m/s'),
        'u_8': 1e3*units('m/s'),
        'm_{core}': 1e2*units('kg/s'),
        'm_{fan}': 1e3*units('kg/s'),
        '\\alpha': 1e1,
        'alphap1': 1e1,
        'm_{total}': 1e3*units('kg/s'),
        'T_2': 1e3*units('K'),
        'P_2': 1e2*units('kPa'),
        'u_2': 1e3*units('m/s'),
        'h_{2}': 1e6*units('J/kg'),
        'T_{2.5}': 1e3*units('K'),
        'P_{2.5}': 1e2*units('kPa'),
        'u_{2.5}': 1e3*units('m/s'),
        'h_{2.5}': 1e6*units('J/kg'),
        'P_{7}': 1e2*units('kPa'),
        'T_{7}': 1e3*units('K'),
        'u_7': 1e3*units('m/s'),
        'P_{5}': 1e2*units('kPa'),
        'T_{5}': 1e3*units('K'),
        'u_5': 1e3*units('m/s'),
        'P_{atm}': 1e2*units('kPa'),
        'T_{atm}': 1e3*units('K'),
        'V': 3e3*units('knot'),
        'a': 1e3*units('m/s'),
    }

    mission = FleetMission(2, 2, 4)
    m = Model(mission['W_{f_{fleet}}'], mission, substitutions)
    sol = m.localsolve(solver='mosek', verbosity = 1)
