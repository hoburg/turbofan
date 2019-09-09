"SP Implementation of the TASOPT engine model"
from __future__ import print_function
from __future__ import absolute_import
from gpkit import (Model, Variable, SignomialsEnabled, units,
                   Vectorize, SignomialEquality)
from gpkit.constraints.tight import Tight
from gpkit.small_scripts import mag
import numpy as np

# Import engine components and maps
from turbofan.maps import FanMap, HPCMap, LPCMap
from turbofan.turbine import Turbine
from turbofan.combustor import Combustor
from turbofan.compressor import Compressor

# Import substitution files
from turbofan.subs import get_737800_subs, get_D82_subs, get_cfm56_subs, get_ge90_subs
from turbofan.test_missions import (TestMissionCFM, TestMissionTASOPT,
                           TestMissionGE90, TestMissionD82, diffs)
from turbofan.initial_guess import initialize_guess

#Cp and gamma values estimated from https://www.ohio.edu/mechanical/thermo/property_tables/air/air_Cp_{c}v.html

class Engine(Model):
    """
    Tasopt engine model
    ________
    INPUTS
    res 7 = 0 = Thrust constrained engine, 1 = burner exit temp/turbine entry temp constrained engine
    cooling = True = cooling model, False = no cooling model
    N = number of discrete flight segments
    state = state model discretized into N segments
    eng = 0 = CFM56 vals, 1 = TASOPT 737-800 vals, 2 = GE90 vals, 3 = TASOPT D8.2 vals, 4 = TASOPT 777-300ER vals
    Nfleet - number of discrete missions in a fleet mission optimization problem, default is 0
    """
    Ttmax = True
    OPRmax = True
    vals = []

    def setup(self, cooling, N, state, eng, Nfleet=0, BLI=False, goption=1):
        """
        setup method for the engine model
        """
        self.constants = EngineConstants(BLI)
        self.setvals(eng, BLI, goption)
        self.compressor = Compressor(self.constants.vals['fexp1'],
                                     self.constants.vals['lpcexp1'],
                                     self.constants.vals['hpcexp1'])
        self.combustor = Combustor(self.constants.vals['ccexp1'],
                                   self.constants.vals['ccexp2'])
        self.turbine = Turbine(self.constants.vals['hptexp1'],
                               self.constants.vals['lptexp1'])
        self.fanmap = FanMap()
        # self.fanmap = FanMap(faneta = vals['faneta'], fgamma = vals['fgamma'])
        self.lpcmap = LPCMap()
        self.hpcmap = HPCMap()
        self.thrust = Thrust()
        self.sizing = Sizing()
        self.state = state

        if Nfleet != 0:
            with Vectorize(Nfleet):
                with Vectorize(N):
                    #-------------------Specified Thrust or Tt4-----------------------
                    self.engineP = self.dynamic(self.state, BLI)
                    OPR = Variable('OPR', '-', 'Overall Pressure Ratio')
        else:
            with Vectorize(N):
                #-------------------Specified Thrust or Tt4-----------------------
                self.engineP = self.dynamic(self.state, BLI)
                OPR = Variable('OPR', '-', 'Overall Pressure Ratio')

        models = [self.compressor , self. combustor, self. turbine, self. thrust, self.fanmap, self.lpcmap, self.hpcmap, self.sizing, self.state, self.engineP]

        #engine weight
        W_engine = Variable('W_{engine}', 'N', 'Weight of a Single Turbofan Engine')

        #engine fan diameter
        df = Variable('d_{f}', 'm', 'Fan Diameter')
        dlpc = Variable('d_{LPC}', 'm', 'LPC Diameter')
        HTRfSub = Variable('HTR_{f_{SUB}}', '-', '1 - HTRf^2')
        HTRlpcSub = Variable('HTR_{lpc_{SUB}}', '-', '1 - HTRlpc^2')

        #OPRmax
        OPRmax = Variable('OPR_{max}', '-', 'Maximum OPR')
        Tt41max = Variable('T_{t_{4.1_{max}}}', 'K', 'Max turbine inlet temperature')

        # Writing constraints
        with SignomialsEnabled():

            weight = [Tight([W_engine/units('kg') >= (self.engineP['m_{total}']/(self.engineP['\\alpha_{+1}']))*
                                          ((1/(100*units('lb/s')))*9.81*units('m/s^2'))*
                                          (1684.5+17.7*(self.engineP['\\pi_{f}']*self.engineP['\\pi_{lc}']*self.engineP['\\pi_{hc}'])/30+
                                           1662.2*(self.engineP['\\alpha']/5)**1.2)])]

            diameter = [df == (4 * self.sizing['A_{2}']/(np.pi * HTRfSub))**.5,
                        dlpc == (4 * self.sizing['A_{2.5}']/(np.pi * HTRlpcSub))**.5,
                        ]

            fmix = [
                # compute f with mixing
                Tight([self.combustor['\\eta_{B}'] * self.engineP['f'] * self.combustor['h_{f}'] >= (1-self.combustor['\\alpha_c'])*self.engineP['h_{t_4}'] - \
                       (1-self.combustor['\\alpha_c'])*self.engineP['h_{t_3}']+self.combustor['C_{p_{fuel}']*self.engineP['f']*(self.engineP['T_{t_4}']-self.combustor['T_{t_f}'])]),
                # compute Tt41...mixing causes a temperature drop
                # had to include Tt4 here to prevent it from being pushed down to zero
                # relaxed SigEq
                Tight([self.engineP['h_{t_{4.1}}']*self.engineP['fp1'] <= ((1-self.combustor['\\alpha_c']+self.engineP['f'])*self.engineP['h_{t_4}'] + self.combustor['\\alpha_c']*self.engineP['h_{t_3}'])]),

                self.engineP['P_{t_4}'] == self.combustor['\\pi_{b}'] * self.engineP['P_{t_3}'],   #B.145
                ]

            fnomix = [
                #only if mixing = false
                #compute f without mixing, overestimation if there is cooling
                Tight([self.combustor['\\eta_{B}'] * self.engineP['f'] * self.combustor['h_{f}'] + self.engineP['h_{t_3}'] >= self.engineP['h_{t_4}']]),

                self.engineP['P_{t_4}'] == self.combustor['\\pi_{b}'] * self.engineP['P_{t_3}'],   #B.145
                ]

            shaftpower = [
                # HPT shaft power balance
                Tight([self.constants['M_{takeoff}']*self.turbine['\eta_{HPshaft}']*(self.engineP['fp1'])*(self.engineP['h_{t_{4.1}}']-self.engineP['h_{t_{4.5}}']) >= self.engineP['h_{t_3}'] - self.engineP['h_{t_{2.5}}']]),    #B.161

                #LPT shaft power balance
                Tight([self.constants['M_{takeoff}']*self.turbine['\eta_{LPshaft}']*(self.engineP['fp1'])*
                (self.engineP['h_{t_{4.5}}'] - self.engineP['h_{t_{4.9}}']) >= self.engineP['h_{t_{2.5}}']-self.engineP['h_{t_{1.8}}']+self.engineP['\\alpha']*(self.engineP['h_{t_{2.1}}'] - self.engineP['h_{t_{2}}'])]),    #B.165
                ]

            hptexit = [
                #HPT Exit states (station 4.5)
                self.engineP['P_{t_{4.5}}'] == self.engineP['\\pi_{HPT}'] * self.engineP['P_{t_{4.1}}'],
                self.engineP['\\pi_{HPT}'] == (self.engineP['T_{t_{4.5}}']/
                                               self.engineP['T_{t_{4.1}}'])**(self.constants.vals['hptexp1']),      #turbine efficiency is 0.9
                ]

            fanmap = [
                self.engineP['\\pi_{f}']*(1.7/self.fanmap['\\pi_{f_D}']) == (1.05*self.engineP['N_{f}']**.0871)**10,
                (self.engineP['\\pi_{f}']*(1.7/self.fanmap['\\pi_{f_D}'])) <= 1.1*(1.06 * (self.engineP['m_{tild_f}'])**0.137)**10,
                (self.engineP['\\pi_{f}']*(1.7/self.fanmap['\\pi_{f_D}'])) >= .9*(1.06 * (self.engineP['m_{tild_f}'])**0.137)**10,

                #define mbar
                self.engineP['m_{f}'] == self.engineP['m_{fan}']*((self.engineP['T_{t_{2}}']/self.constants['T_{ref}'])**.5)/(self.engineP['P_{t_{2}}']/self.constants['P_{ref}']),    #B.280

                self.engineP['\\pi_{f}'] >= 1,

                ]

            lpcmap = [
                self.engineP['\\pi_{lc}']*(26/self.lpcmap['\pi_{lc_D}']) == (1.38 * (self.engineP['N_{1}'])**0.566)**10,
                self.engineP['\\pi_{lc}']*(26/self.lpcmap['\pi_{lc_D}']) <= 1.1*(1.38 * (self.engineP['m_{tild_lc}'])**0.122)**10,
                self.engineP['\\pi_{lc}']*(26/self.lpcmap['\pi_{lc_D}']) >= .9*(1.38 * (self.engineP['m_{tild_lc}'])**0.122)**10,

                self.engineP['\\pi_{lc}'] >= 1,

                #define mbar..technially not needed b/c constrained in res 2 and/or 3
                self.engineP['m_{lc}'] == self.engineP['m_{core}']*((self.engineP['T_{t_{2}}']/self.constants['T_{ref}'])**.5)/(self.engineP['P_{t_{2}}']/self.constants['P_{ref}']),    #B.280
                ]

            hpcmap = [
                self.engineP['\\pi_{hc}']*(26/self.hpcmap['\\pi_{hc_D}']) == (1.38 * (self.engineP['N_{2}'])**0.566)**10,
                self.engineP['\\pi_{hc}']*(26/self.hpcmap['\\pi_{hc_D}']) >= .9*(1.38 * (self.engineP['m_{tild_{hc}}'])**0.122)**10,
                self.engineP['\\pi_{hc}']*(26/self.hpcmap['\\pi_{hc_D}']) <= 1.1*(1.38 * (self.engineP['m_{tild_{hc}}'])**0.122)**10,

                self.engineP['m_{hc}'] == self.engineP['m_{core}']*((self.engineP['T_{t_{2.5}}']/self.constants['T_{ref}'])**.5)/(self.engineP['P_{t_{2.5}}']/self.constants['P_{ref}']),    #B.280

                self.engineP['\\pi_{hc}'] >= 1,
                ]

            opr = [
                OPR == self.engineP['\\pi_{hc}']*self.engineP['\\pi_{lc}']*self.engineP['\\pi_{f}'],
                ]

            thrust = [
                self.engineP['P_{t_8}'] == self.engineP['P_{t_7}'], #B.179
                self.engineP['T_{t_8}'] == self.engineP['T_{t_7}'], #B.180

                self.engineP['P_{t_6}'] == self.engineP['P_{t_5}'], #B.183
                self.engineP['T_{t_6}'] == self.engineP['T_{t_5}'], #B.184

                Tight([self.engineP['F_{6}']/(self.constants['M_{takeoff}']*self.engineP['m_{core}']) + (self.engineP['fp1'])*self.state['V'] <= (self.engineP['fp1'])*self.engineP['u_{6}']]),

                #ISP
                self.engineP['I_{sp}'] == self.engineP['F_{sp}']*self.state['a']*(self.engineP['\\alpha_{+1}'])/(self.engineP['f']*self.constants['g']),  #B.192
                ]

            res1 = [
                #residual 1 Fan/LPC speed
                self.engineP['N_{f}']*self.sizing['G_{f}'] == self.engineP['N_{1}'],
                self.engineP['N_{1}'] <= 1.1,
                self.engineP['N_{2}'] <= 1.1,
                ]

                #note residuals 2 and 3 differ from TASOPT, by replacing mhc with mlc
                #in residual 4 I was able to remove the LPC/HPC mass flow equality
                #in residual 6 which allows for convergence

            res2 = [
                #residual 2 HPT mass flow
                self.sizing['m_{htD}'] == (self.engineP['fp1'])*self.engineP['m_{hc}']*self.constants['M_{takeoff}']*
                (self.engineP['P_{t_{2.5}}']/self.engineP['P_{t_{4.1}}'])*
                (self.engineP['T_{t_{4.1}}']/self.engineP['T_{t_{2.5}}'])**.5,
                ]

            res3 = [
                #residual 3 LPT mass flow
                (self.engineP['fp1'])*self.engineP['m_{lc}']*self.constants['M_{takeoff}']*
                (self.engineP['P_{t_{1.8}}']/self.engineP['P_{t_{4.5}}'])*
                (self.engineP['T_{t_{4.5}}']/self.engineP['T_{t_{1.8}}'])**.5
                == self.sizing['m_{ltD}'],
                ]

            res4 = [
                #residual 4
                (self.engineP['P_{7}']/self.engineP['P_{t_7}']) == (self.engineP['T_{7}']/self.engineP['T_{t_7}'])**(3.5),
                Tight([(self.engineP['T_{7}']/self.engineP['T_{t_7}'])**-1 >= 1 + .2 * self.engineP['M_7']**2]),
                ]

            res5 = [
                #residual 5 core nozzle mass flow
                (self.engineP['P_{5}']/self.engineP['P_{t_5}']) == (self.engineP['T_{5}']/self.engineP['T_{t_5}'])**(3.583979),
                Tight([(self.engineP['T_{5}']/self.engineP['T_{t_5}'])**-1 >= 1 + .2 * self.engineP['M_5']**2]),
                ]


            massflux = [
                #compute core mass flux
                self.constants['M_{takeoff}'] * self.engineP['m_{core}'] == self.engineP['\\rho_5'] * self.sizing['A_{5}'] * self.engineP['u_{5}']/(self.engineP['fp1']),

                #compute fan mas flow
                self.engineP['m_{fan}'] == self.engineP['\\rho_7']*self.sizing['A_{7}']*self.engineP['u_{7}'],

                Tight([self.engineP['m_{total}'] >= self.engineP['m_{fan}'] + self.engineP['m_{core}']]), # [SP] # [SigEq]
                ]

            #component area sizing
            fanarea = [
                #fan area
                self.engineP['P_{2}'] == self.engineP['P_{t_{2}}']*(self.engineP['hold_{2}'])**(-3.512),
                self.engineP['T_{2}'] == self.engineP['T_{t_{2}}'] * self.engineP['hold_{2}']**-1,
                self.sizing['A_{2}'] == self.engineP['m_{fan}']/(self.engineP['\\rho_{2}']*self.engineP['u_{2}']),     #B.198
                ]

            HPCarea = [
                #HPC area
                self.engineP['P_{2.5}'] == self.engineP['P_{t_{2.5}}']*(self.engineP['hold_{2.5}'])**(-3.824857),
                self.engineP['T_{2.5}'] == self.engineP['T_{t_{2.5}}'] * self.engineP['hold_{2.5}']**-1,
                self.sizing['A_{2.5}'] == self.engineP['m_{core}']/(self.engineP['\\rho_{2.5}']*self.engineP['u_{2.5}']),     #B.203
                ]

            if eng == 0:
                """CFM56 vals"""
                onDest = [
                    #estimate relevant on design values
                    self.sizing['m_{htD}'] <= 1.3*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}'] *((1400.0/288)**.5)/(1527/101.325),
                    self.sizing['m_{htD}'] >= .7*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}'] *((1400.0/288)**.5)/(1527/101.325),
                    self.sizing['m_{ltD}'] <= 1.3*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}'] *((1038.8/288)**.5)/(589.2/101.325),
                    self.sizing['m_{ltD}'] >= .7*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}'] *((1038.8/288)**.5)/(589.2/101.325),
                    self.lpcmap['m_{lc_D}'] >= .7*self.sizing['m_{coreD}']*((292.57/288)**.5)/(84.25/101.325),
                    self.lpcmap['m_{lc_D}'] <= 1.3*self.sizing['m_{coreD}'] *((292.57/288)**.5)/(84.25/101.325),
                    self.hpcmap['m_{hc_D}'] >= .7*self.sizing['m_{coreD}'] *((362.47/288)**.5)/(163.02/101.325),
                    self.hpcmap['m_{hc_D}'] <= 1.3*self.sizing['m_{coreD}'] *((362.47/288)**.5)/(163.02/101.325),
                    self.fanmap['\\bar{m}_{fan_{D}}'] >= .7 * self.sizing['\\alpha_{OD}'] * self.sizing['m_{coreD}'] *((250.0/288)**.5)/(50/101.325),
                    self.fanmap['\\bar{m}_{fan_{D}}'] <= 1.3 * self.sizing['\\alpha_{OD}'] * self.sizing['m_{coreD}']* ((250.0/288)**.5)/(50/101.325),
                    ]
            if eng == 1:
                """TASOPT 737-800 vals"""
                onDest = [
                    self.sizing['m_{htD}'] <= 1.3*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}'] *((1400.0/288)**.5)/(1498/101.325),
                    self.sizing['m_{htD}'] >= .7*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}']*((1400.0/288)**.5)/(1498/101.325),
                    self.sizing['m_{ltD}'] <= 1.3*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}'] *((1144.8/288)**.5)/(788.5/101.325),
                    self.sizing['m_{ltD}'] >= .7*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}']*((1144.8/288)**.5)/(788.5/101.325),
                    self.lpcmap['m_{lc_D}'] >= .7*self.sizing['m_{coreD}']*((294.5/288)**.5)/(84.25/101.325),
                    self.lpcmap['m_{lc_D}'] <= 1.3*self.sizing['m_{coreD}'] *((294.5/288)**.5)/(84.25/101.325),
                    self.hpcmap['m_{hc_D}'] >= .7*self.sizing['m_{coreD}'] *((482.7/288)**.5)/(399.682/101.325),
                    self.hpcmap['m_{hc_D}'] <= 1.3*self.sizing['m_{coreD}'] *((482.7/288)**.5)/(399.682/101.325),
                    self.fanmap['\\bar{m}_{fan_{D}}'] >= .7 * self.sizing['\\alpha_{OD}'] * self.sizing['m_{coreD}'] *((250.0/288)**.5)/(50/101.325),
                    self.fanmap['\\bar{m}_{fan_{D}}'] <= 1.3 * self.sizing['\\alpha_{OD}'] * self.sizing['m_{coreD}']* ((250.0/288)**.5)/(50/101.325),
                    ]
            if eng == 2:
                """GE90 vals"""
                onDest = [
                    #estimate relevant on design values
                    self.sizing['m_{htD}'] <= 1.3*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}'] *((1400.0/288)**.5)/(1527/101.325),
                    self.sizing['m_{htD}'] >= .7*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}'] *((1400.0/288)**.5)/(1527/101.325),
                    self.sizing['m_{ltD}'] <= 1.3*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}'] *((1038.8/288)**.5)/(589.2/101.325),
                    self.sizing['m_{ltD}'] >= .7*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}'] *((1038.8/288)**.5)/(589.2/101.325),
                    self.lpcmap['m_{lc_D}'] >= .7*self.sizing['m_{coreD}']*((292.57/288)**.5)/(84.25/101.325),
                    self.lpcmap['m_{lc_D}'] <= 1.3*self.sizing['m_{coreD}'] *((292.57/288)**.5)/(84.25/101.325),
                    self.hpcmap['m_{hc_D}'] >= .7*self.sizing['m_{coreD}'] *((362.47/288)**.5)/(163.02/101.325),
                    self.hpcmap['m_{hc_D}'] <= 1.3*self.sizing['m_{coreD}'] *((362.47/288)**.5)/(163.02/101.325),
                    self.fanmap['\\bar{m}_{fan_{D}}'] >= .3 * self.sizing['\\alpha_{OD}'] * self.sizing['m_{coreD}'] *((250.0/288)**.5)/(50/101.325),
                    self.fanmap['\\bar{m}_{fan_{D}}'] <= 1.7 * self.sizing['\\alpha_{OD}'] * self.sizing['m_{coreD}']* ((250.0/288)**.5)/(50/101.325),
                ]

            if eng == 3:
                """TASOPT D8.2 vals"""
                if BLI:
                   onDest = [
                        self.sizing['m_{htD}'] <= 1.3*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}'] *((1400.0/288)**.5)/(1433.49/101.325),
                        self.sizing['m_{htD}'] >= .7*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}']*((1400.0/288)**.5)/(1433.49/101.325),
                        self.sizing['m_{ltD}'] <= 1.3*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}'] *((1121.85/288)**.5)/(706.84/101.325),
                        self.sizing['m_{ltD}'] >= .7*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}']*((1121.85/288)**.5)/(706.84/101.325),
                        self.lpcmap['m_{lc_D}'] >= .7*self.sizing['m_{coreD}']*((289.77/288)**.5)/(65.79434/101.325),
                        self.lpcmap['m_{lc_D}'] <= 1.3*self.sizing['m_{coreD}'] *((289.77/288)**.5)/(65.79434/101.325),
                        self.hpcmap['m_{hc_D}'] >= .7*self.sizing['m_{coreD}'] *((481.386/288)**.5)/(327.66/101.325),
                        self.hpcmap['m_{hc_D}'] <= 1.3*self.sizing['m_{coreD}'] *((481.386/288)**.5)/(327.66/101.325),
                        self.fanmap['\\bar{m}_{fan_{D}}'] >= .7 * self.sizing['\\alpha_{OD}'] * self.sizing['m_{coreD}'] *((250.0/288)**.5)/(41.0/101.325),
                        self.fanmap['\\bar{m}_{fan_{D}}'] <= 1.3 * self.sizing['\\alpha_{OD}'] * self.sizing['m_{coreD}']* ((250.0/288)**.5)/(41.0/101.325),
                        ]
                else:
                    onDest = [
                        self.sizing['m_{htD}'] <= 1.3*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}'] *((1400.0/288)**.5)/(1598.32/101.325),
                        self.sizing['m_{htD}'] >= .7*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}']*((1400.0/288)**.5)/(1598.32/101.325),
                        self.sizing['m_{ltD}'] <= 1.3*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}'] *((1142.6/288)**.5)/(835.585/101.325),
                        self.sizing['m_{ltD}'] >= .7*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}']*((1142.6/288)**.5)/(835.585/101.325),
                        self.lpcmap['m_{lc_D}'] >= .7*self.sizing['m_{coreD}']*((289.77/288)**.5)/(80.237/101.325),
                        self.lpcmap['m_{lc_D}'] <= 1.3*self.sizing['m_{coreD}'] *((289.77/288)**.5)/(80.237/101.325),
                        self.hpcmap['m_{hc_D}'] >= .7*self.sizing['m_{coreD}'] *((481.386/288)**.5)/(399.58/101.325),
                        self.hpcmap['m_{hc_D}'] <= 1.3*self.sizing['m_{coreD}'] *((481.386/288)**.5)/(399.58/101.325),
                        self.fanmap['\\bar{m}_{fan_{D}}'] >= .7 * self.sizing['\\alpha_{OD}'] * self.sizing['m_{coreD}'] *((250.0/288)**.5)/(50/101.325),
                        self.fanmap['\\bar{m}_{fan_{D}}'] <= 1.3 * self.sizing['\\alpha_{OD}'] * self.sizing['m_{coreD}']* ((250.0/288)**.5)/(50/101.325),
                        ]

            if eng == 4:
                """TASOPT 777-300ER"""
                onDest = [
                    #estimate relevant on design values
                    self.sizing['m_{htD}'] <= 1.3*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}'] *((1600.0/288)**.5)/(2094.48/101.325),
                    self.sizing['m_{htD}'] >= .7*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}'] *((1600.0/288)**.5)/(2094.48/101.325),
                    self.sizing['m_{ltD}'] <= 1.3*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}'] *((1277.99/288)**.5)/(1022.19/101.325),
                    self.sizing['m_{ltD}'] >= .7*self.engineP['fp1']*self.constants['M_{takeoff}']*self.sizing['m_{coreD}'] *((1277.99/288)**.5)/(1022.19/101.325),
                    self.lpcmap['m_{lc_D}'] >= .7*self.sizing['m_{coreD}']*((289.26/288)**.5)/(79.79/101.325),
                    self.lpcmap['m_{lc_D}'] <= 1.3*self.sizing['m_{coreD}'] *((289.26/288)**.5)/(79.79/101.325),
                    self.hpcmap['m_{hc_D}'] >= .7*self.sizing['m_{coreD}'] *((481.15/288)**.5)/(398.95/101.325),
                    self.hpcmap['m_{hc_D}'] <= 1.3*self.sizing['m_{coreD}'] *((481.15/288)**.5)/(398.95/101.325),
                    self.fanmap['\\bar{m}_{fan_{D}}'] >= .3 * self.sizing['\\alpha_{OD}'] * self.sizing['m_{coreD}'] *((250.0/288)**.5)/(50/101.325),
                    self.fanmap['\\bar{m}_{fan_{D}}'] <= 1.7 * self.sizing['\\alpha_{OD}'] * self.sizing['m_{coreD}']* ((250.0/288)**.5)/(50/101.325),
                    ]

        if self.OPRmax:
            opr.extend([OPR <= OPRmax])

        constraints = [weight, diameter, shaftpower, hptexit,
                       fanmap, lpcmap, hpcmap, opr, thrust,
                       res1, res2, res3, res4, res5,
                       massflux, fanarea, HPCarea, onDest]

        if cooling:
            constraints += fmix
        else:
            constraints += fnomix

        return models, constraints

    def dynamic(self, state, BLI):
        """
        creates an instance of the engine performance model
        """
        return EnginePerformance(self, state, BLI)

    def setvals(self, eng, BLI, goption):
        """

        :param eng: 0,..., 4 determines the component efficiencies depending
                    on engine choice
        :param BLI: boolean whether the engine ingests boundary layer
        :param goption: sets gamma value
        :return: dictionary of exponents and components efficiencies
        """
        vd = {}
        vd['eng'] = eng
        vd['BLI'] = BLI
        if BLI:
            vd['fan_eta_reduct'] = 0.96
        else:
            vd['fan_eta_reduct'] = 1.0

        # Setting engine exponents in sequence of...
        # Component gas constants
        gammakeys = ['fgamma', 'lpcgamma', 'hpcgamma', 'ccgamma', 'lptgamma', 'hptgamma',
                     'sta6gamma', 'sta8gamma']
        # Component efficiencies
        etakeys = ['faneta', 'LPCeta', 'HPCeta', 'HPTeta', 'LPTeta']
        # Exponents
        expkeys = ['fexp1', 'lpcexp1', 'hpcexp1',
                   'ccexp1', 'ccexp2', 'lptexp1',
                   'hptexp1', 'fanexexp', 'turbexexp']
        goption0 = [1.401, 1.398, 1.354, 1.313, 1.306, 1.299, 1.4, 1.387]
        goption1 = [1.401, 1.398, 1.354, 1.313, 1.354, 1.318, 1.4, 1.387]
        # Different engine numbers set the component efficiencies for the different engines
        etas = dict(zip([0,1,2,3,4],[
                        [.9005, .9306, .9030, .8731, .8851], #CFM56
                        [.8948, .8800, .8700, .8990, .8890], #TASOPT CFM56 for 737-800
                        [.9153, .9037, .9247, .9121, .9228], #GE90
                        [.9300, .9200, .8900, .9100, .9200], #D8, 37 size
                        [.9100, .9000, .8900, .9000, .9000]]))#TASOPT 777-300ER
        # Note: ccgamma = gamma value of air @ 1400 K
        if goption == 0:
            vd.update(dict.zip(gammakeys, goption0))
        else:
            vd.update(dict(zip(gammakeys, goption1)))

        # Updating component efficiencies depending on engine
        vd.update(dict(zip(etakeys, etas[eng])))

        exps = [# fan exponent
                (vd['fgamma'] - 1)/(vd['faneta'] * vd['fgamma']),
                # compressor exponents
                (vd['lpcgamma'] - 1)/(vd['LPCeta'] * vd['lpcgamma']),
                (vd['hpcgamma'] - 1)/(vd['HPCeta'] * vd['hpcgamma']),
                # combustor cooling exponents
                vd['ccgamma']/(1 - vd['ccgamma']),
                -vd['ccgamma']/(1 - vd['ccgamma']),
                # turbine exponents
                vd['lptgamma'] * vd['LPTeta'] / (vd['lptgamma'] - 1),
                vd['hptgamma'] * vd['HPTeta'] / (vd['hptgamma'] - 1),
                # fan exit exponent (station 8)
                (vd['sta8gamma'] - 1)/ vd['sta8gamma'],
                # turbine exit exponent (station 6)
                (vd['sta6gamma'] - 1) / vd['sta6gamma']]
        vd.update(dict(zip(expkeys, exps)))
        self.constants.vals = vd

class EnginePerformance(Model):
    """
    Engine performance model
    """
    def setup(self, engine, state, BLI, **kwargs):

        #create the subcomponent performance models
        self.compP = engine.compressor.dynamic(engine.constants, state, BLI)
        self.combP = engine.combustor.dynamic(engine.constants, state)
        self.turbineP = engine.turbine.dynamic(engine.constants)
        self.thrustP = engine.thrust.dynamic(engine.constants, state, BLI)
        self.fanmapP = engine.fanmap.dynamic(engine.constants)
        self.lpcmapP = engine.lpcmap.dynamic(engine.constants)
        self.hpcmapP = engine.hpcmap.dynamic(engine.constants)
        self.sizingP = engine.sizing.dynamic(engine.constants, engine.compressor, engine.fanmap, engine.lpcmap, engine.hpcmap, state)

        models = [self.compP, self.combP, self.turbineP, self.thrustP, self.fanmapP, self.lpcmapP, self.hpcmapP, self.sizingP]

        return models

class EngineConstants(Model):
    """
    Class of constants used in the engine model
    """
    def setup(self, BLI):
        #-----------------------air properties------------------
        #ambient
        R = Variable('R', 287, 'J/kg/K', 'R', constant=True)

        #gravity
        g = Variable('g', 9.81, 'm/(s^2)', 'Gravitational Acceleration', constant=True)

        #-------------------------reference temp and pressure--------------------
        Tref = Variable('T_{ref}', 'K', 'Reference Stagnation Temperature')
        Pref = Variable('P_{ref}', 'kPa', 'Reference Stagnation Pressure')

        #---------------------------efficiencies & takeoffs-----------------------
        Mtakeoff = Variable('M_{takeoff}', '-', '1 Minus Percent mass flow loss for de-ice, pressurization, etc.')

        #----------------BLI pressure loss factor-------------------
        if BLI:
            fBLIP = Variable('f_{BLI_{P}}', '-', 'BLI Stagnation Pressure Loss Ratio')
            fBLIV = Variable('f_{BLI_{V}}', '-', 'BLI Velocity Loss Ratio')



class Thrust(Model):
    """"
    thrust sizing model
    """
    def setup(self):
        #define new variables
        #fan and exhaust
        Cptex = Variable('C_{p_{tex}', 1029, 'J/kg/K', "Cp Value for Combustion Products at Core Exhaust") #500K, gamma = 1.387
        Cpfanex = Variable('C_{p_{fex}', 1005, 'J/kg/K', "Cp Value for Air at 300K") #gamma =1.4        #heat of combustion of jet fuel

        #max by pass ratio
        alpha_max = Variable('\\alpha_{max}', '-', 'Maximum bypass ratio')

    def dynamic(self, engine, state, BLI):
        """
        creates an instance of the thrust performance model
        """
        return ThrustPerformance(self, engine, state, BLI)

class ThrustPerformance(Model):
    """
    thrust performance model
    """
    def setup(self, thrust, engine, state, BLI):
        self.thrust = thrust
        self.engine = engine

        #define new variables
        #------------------fan exhaust (station 8) statge variables------------
        P8 = Variable('P_{8}', 'kPa', 'Fan Exhaust Static Pressure')
        Pt8 = Variable('P_{t_8}', 'kPa', 'Fan Exhaust Stagnation Pressure')
        ht8 = Variable('h_{t_8}', 'J/kg', 'Fan Exhaust Stagnation Enthalpy')
        h8 = Variable('h_{8}', 'J/kg', 'Fan Exhasut Static Enthalpy')
        Tt8 = Variable('T_{t_8}', 'K', 'Fan Exhaust Stagnation Temperature (8)')
        T8 = Variable('T_{8}', 'K', 'Fan Exhaust Static Temperature (8)')

        #-----------------core exhaust (station 6) state variables-------------
        P6 = Variable('P_{6}', 'kPa', 'Core Exhaust Static Pressure')
        Pt6 = Variable('P_{t_6}', 'kPa', 'Core Exhaust Stagnation Pressure')
        Tt6 = Variable('T_{t_6}', 'K', 'Core Exhaust Stagnation Temperature (6)')
        T6 = Variable('T_{6}', 'K', 'Core Exhaust Static Temperature (6)')
        ht6 = Variable('h_{t_6}', 'J/kg', 'Core Exhaust Stagnation Enthalpy')
        h6 = Variable('h_6', 'J/kg', 'Core Exhaust Static Enthalpy')

        #thrust variables
        F8 = Variable('F_{8}', 'N', 'Fan Thrust')
        F6 = Variable('F_{6}', 'N', 'Core Thrust')
        F = Variable('F', 'N', 'Total Thrust')
        Fsp = Variable('F_{sp}', '-', 'Specific Net Thrust')
        Isp = Variable('I_{sp}', 's', 'Specific Impulse')
        TSFC = Variable('TSFC', '1/hr', 'Thrust Specific Fuel Consumption')

        #exhaust speeds
        u6 = Variable('u_{6}', 'm/s', 'Core Exhaust Velocity')
        u8 = Variable('u_{8}', 'm/s', 'Fan Exhaust Velocity')

        #mass flows
        mCore = Variable('m_{core}', 'kg/s', 'Core Mass Flow')
        mFan = Variable('m_{fan}', 'kg/s', 'Fan Mass Flow')

        #------------------By-Pass Ratio (BPR)----------------------------
        alpha = Variable('\\alpha', '-', 'By Pass Ratio')
        alphap1 = Variable('\\alpha_{+1}', '-', '1 plus BPR')

        #constraints
        constraints = []

        with SignomialsEnabled():
            #exhaust and thrust constraints
            constraints.extend([
                P8 == state["P_{atm}"],
                h8 == self.thrust['C_{p_{fex}'] * T8,
                Tight([u8**2 + 2*h8 <= 2*ht8]),
                (P8/Pt8)**(self.engine.vals['fanexexp']) == T8/Tt8,
                ht8 == self.thrust['C_{p_{fex}'] * Tt8,

                #core exhaust
                P6 == state["P_{atm}"],   #B.4.11 intro
                (P6/Pt6)**(self.engine.vals['turbexexp']) == T6/Tt6,
                Tight([u6**2 + 2*h6 <= 2*ht6]),
                h6 == self.thrust['C_{p_{tex}'] * T6,
                ht6 == self.thrust['C_{p_{tex}'] * Tt6,

                #constrain the new BPR
                alpha == mFan / mCore,
                SignomialEquality(alphap1, alpha + 1),
                alpha <= self.thrust['\\alpha_{max}'],

                #SIGNOMIAL
                Tight([F <= F6 + F8]),

                Fsp == F/((alphap1)*mCore*state['a']),   #B.191

                #TSFC
                TSFC == 1/Isp,
                ])

        if BLI:
            constraints.extend([
                #overall thrust values
                Tight([F8/(alpha * mCore) + state['V']*engine['f_{BLI_{V}}'] <= u8]),  #B.188
                ])

        else:
            constraints.extend([
                #overall thrust values
                Tight([F8/(alpha * mCore) + state['V'] <= u8]),  #B.188
                ])

        return constraints

class Sizing(Model):
    """"
    engine sizing model
    """
    def setup(self):
        #define new variables
        #gear ratio, set to 1 if no gearing present
        Gf = Variable('G_{f}', '-', 'Gear Ratio Between Fan and LPC')

        mhtD = Variable('m_{htD}', 'kg/s', 'Design HPT Corrected Mass Flow (see B.225)')
        mltD = Variable('m_{ltD}', 'kg/s', 'Design LPT Corrected Mass Flow (see B.226)')

        #on design by pass ratio
        alpha_max = Variable('\\alpha_{OD}', '-', 'By Pass Ratio')

        #-------------------------Areas------------------------
        A2 = Variable('A_{2}', 'm^2', 'Fan Area')
        A25 = Variable('A_{2.5}', 'm^2', 'HPC Area')
        A5 = Variable('A_{5}', 'm^2', 'Core Exhaust Nozzle Area')
        A7 = Variable('A_{7}', 'm^2', 'Fan Exhaust Nozzle Area')

        mCoreD = Variable('m_{coreD}', 'kg/s', 'Estimated on Design Mass Flow')

        #constraints
        constraints = []

        constraints.extend([
            #-------------------------Areas------------------------
            Tight([A5 + A7 <= A2]),
            ])

        return constraints

    def dynamic(self, engine, compressor, fanmap, lpcmap, hpcmap, state):
        """
        creates an instance of the engine sizing performance model
        """
        return SizingPerformance(self, engine, compressor, fanmap, lpcmap, hpcmap, state)

class SizingPerformance(Model):
    """
    engine sizing performance model
    """
    def setup(self, sizing, engine, compressor, fanmap, lpcmap, hpcmap, state):
        self.sizing = sizing
        self.engine = engine
        self.compressor = compressor
        self.fanmap = fanmap
        self.lpcmap = lpcmap
        self.hpcmap = hpcmap

        #new variables
        #exhaust mach numbers
        a5 = Variable('a_{5}', 'm/s', 'Speed of Sound at Station 5')
        a7 = Variable('a_{7}', 'm/s', 'Speed of Sound at Station 7')

        #mass flows
        mtot = Variable('m_{total}', 'kg/s', 'Total Engine Mass Flux')

        #-------------------fan face variables---------------------
        rho2 = Variable('\\rho_{2}', 'kg/m^3', 'Air Static Density at Fan Face')
        T2 = Variable('T_{2}', 'K', 'Air Static Temperature at Fan Face')
        P2 = Variable('P_{2}', 'kPa', 'Air Static Pressure at Fan Face')
        u2 = Variable('u_{2}', 'm/s', 'Air Speed at Fan Face')
        h2 = Variable('h_{2}', 'J/kg', 'Static Enthalpy at the Fan Inlet (2)')
        M2 = Variable('M_2', '-', 'Fan Face/LPC Face Axial Mach Number')

        #------------------HPC face variables---------------------
        rho25 = Variable('\\rho_{2.5}', 'kg/m^3', 'Static Air Density at HPC Face')
        T25 = Variable('T_{2.5}', 'K', 'Static Air Temperature at HPC Face')
        P25 = Variable('P_{2.5}', 'kPa', 'Static Air Pressure at HPC Face')
        u25 = Variable('u_{2.5}', 'm/s', 'Air Speed at HPC Face')
        h25 = Variable('h_{2.5}', 'J/kg', 'Static Enthalpy at the LPC Exit (2.5)')
        M25 = Variable('M_{2.5}', '-', 'HPC Face Axial Mach Number')

        #fan exhuast states
        P7 = Variable('P_{7}', 'kPa', 'Fan Exhaust Static Pressure (7)')
        T7 = Variable('T_{7}', 'K', 'Static Temperature at the Fan Nozzle Exit (7)')
        u7 = Variable('u_{7}', 'm/s', 'Station 7 Exhaust Velocity')
        M7 = Variable('M_7', '-', 'Station 7 Mach Number')
        rho7 = Variable('\\rho_7', 'kg/m^3', 'Air Static Density at Fam Exhaust Exit (7)')

        #core exhaust states
        P5 = Variable('P_{5}', 'kPa', 'Core Exhaust Static Pressure (5)')
        T5 = Variable('T_{5}', 'K', 'Static Temperature at the Turbine Nozzle Exit (5)')
        M5 = Variable('M_5', '-', 'Station 5 Mach Number')
        u5 = Variable('u_{5}', 'm/s', 'Station 5 Exhaust Velocity')
        rho5 = Variable('\\rho_5', 'kg/m^3', 'Air Static Density at Core Exhaust Exit (5)')

        #constraints
        constraints = []

        #sizing constraints that hold the engine together
        constraints.extend([
            #residual 4
            P7 >= state["P_{atm}"],
            M7 <= 1,
            a7 == (1.4*self.engine['R']*T7)**.5,
            a7*M7 == u7,
            rho7 == P7/(self.engine['R']*T7),

            #residual 5 core nozzle mass flow
            P5 >= state["P_{atm}"],
            M5 <= 1,
            a5 == (1.387*self.engine['R']*T5)**.5,
            a5*M5 == u5,
            rho5 == P5/(self.engine['R']*T5),

            #component area sizing
            #fan area
            h2 == self.compressor['C_{p_{1}'] * T2,
            rho2 == P2/(self.engine['R'] * T2),  #B.196
            u2 == M2*(self.compressor['C_{p_{1}']*self.engine['R']*T2/(781.*units('J/kg/K')))**.5,  #B.197

            #HPC area
            h25 == self.compressor['C_{p_{2}'] * T25,
            rho25 == P25/(self.engine['R']*T25),
            u25 == M25*(self.compressor['C_{p_{2}']*self.engine['R']*T25/(781.*units('J/kg/K')))**.5,   #B.202
        ])

        return constraints

class TestState(Model):
    """
    state class only to be used for testing purposes
    """
    def setup(self):
        #define variables
        p_atm = Variable("P_{atm}", "kPa", "air pressure")
        T_atm = Variable("T_{atm}", "K", "air temperature")

        V = Variable('V', 'kts', 'Aircraft Flight Speed')
        a = Variable('a', 'm/s', 'Speed of Sound')
        R = Variable('R', 287, 'J/kg/K', 'Air Specific Heat')
        gamma = Variable('\\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        M = Variable('M', '-', 'Mach Number')

        #make constraints
        constraints = []

        constraints.extend([
            V == M * a,
            a  == (gamma * R * T_atm)**.5,
            ])

        return constraints

def test():
    """
    Test the cfm56 engine model
    """
    with Vectorize(2):
        state = TestState()

    engine = Engine(True, 2, state, 0)

    mission = TestMissionCFM(engine)

    substitutions = get_cfm56_subs()
    m = Model((10*engine.engineP.thrustP['TSFC'][0]+engine.engineP.thrustP['TSFC'][1]), [engine, mission], substitutions)
    sol = m.localsolve(verbosity = 2)

if __name__ == "__main__":
    """
    eng = 0 is CFM56, set N = 2
    eng = 1 is TASOPT 737-800, set N = 3
    eng = 2 is GE90, set N = 2
    eng = 3 is TASOPT D8.2, set N=2
    """
    eng = 0

    if eng == 0 or eng == 2 or eng == 3:
        N = 2
    if eng == 1:
        N = 3

    with Vectorize(N):
        state = TestState()

    engine = Engine(True, N, state, eng)

    if eng == 0:
        mission = TestMissionCFM(engine)
        substitutions = get_cfm56_subs()

    if eng == 1:
        mission = TestMissionTASOPT(engine)
        substitutions = get_737800_subs()

    if eng == 2:
        mission = TestMissionGE90(engine)
        substitutions = get_ge90_subs()

    if eng == 3:
        mission = TestMissionD82(engine)
        substitutions = get_D82_subs()

    #dict of initial guesses
    x0 = initialize_guess()

    #select the proper objective based off of the number of flight segments
    if eng == 0 or eng == 2 or eng == 3:
        m = Model(engine['W_{engine}']**0.5*np.dot([10.,1.],engine.engineP.thrustP['TSFC']), [engine, mission], substitutions, x0=x0)
    if eng == 1:
        m = Model(engine['W_{engine}']**0.5*np.dot([10., 1., 1.], engine.engineP.thrustP['TSFC']), [engine, mission], substitutions, x0=x0)

    #update substitutions and solve
    m.substitutions.update(substitutions)
    sol = m.localsolve(verbosity=2, mutategp=False) #x0=x0)
    sol.savetxt()

    #print out various percent differences in TSFC and engine areas
    diffs(sol, eng)
