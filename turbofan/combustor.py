from gpkit import Variable, Model, SignomialsEnabled, SignomialEquality
from gpkit.constraints.tight import Tight

class Combustor(Model):
    """"
    combustor model
    """
    def setup(self, ccexp1, ccexp2):
        self.ccexp1 = ccexp1
        self.ccexp2 = ccexp2
        #define new variables
        Cpc = Variable('C_{p_{c}}', 1216, 'J/kg/K', "Cp Value for Fuel/Air Mix in Combustor") #1400K, gamma equals 1.312
        Cpfuel = Variable('C_{p_{fuel}', 2010, 'J/kg/K', 'Specific Heat Capacity of Kerosene (~Jet Fuel)')
        hf = Variable('h_{f}', 43.003, 'MJ/kg', 'Heat of Combustion of Jet Fuel')     #http://hypeRbook.com/facts/2003/EvelynGofman.shtml...prob need a better source

        #-------------------------diffuser pressure ratios--------------------------
        pib = Variable('\\pi_{b}', '-', 'Burner Pressure Ratio')

        #---------------------------efficiencies & takeoffs-----------------------
        etaB = Variable('\\eta_{B}', '-', 'Burner Efficiency')

        #------------------------Variables for cooling flow model---------------------------
        #cooling flow bypass ratio
        ac = Variable('\\alpha_c', '-', 'Total Cooling Flow Bypass Ratio')
        #variables for cooling flow velocity
        ruc = Variable('r_{uc}', '-', 'User Specified Cooling Flow Velocity Ratio')

        hold4a = Variable('hold_{4a}', '-', '1+(gamma-1)/2 * M_4a**2')

        Ttf = Variable('T_{t_f}', 'K', 'Incoming Fuel Total Temperature')

    def dynamic(self, engine, state):
        """
        creates an instance of the fan map performance model
        """
        return CombustorPerformance(self, engine, state)

class CombustorPerformance(Model):
    """
    combustor performance constraints
    """
    def setup(self, combustor, engine, state, mixing = True):
        self.combustor = combustor
        self.engine = engine

        R = Variable('R', 287, 'J/kg/K', 'Air Specific Heat')

        #define new variables
        #--------------------------combustor exit (station 4) stagnation states------------------
        Pt4 = Variable('P_{t_4}', 'kPa', 'Stagnation Pressure at the Combustor Exit (4)')
        ht4 = Variable('h_{t_4}', 'J/kg', 'Stagnation Enthalpy at the Combustor Exit (4)')
        Tt4 = Variable('T_{t_4}', 'K', 'Combustor Exit (Station 4) Stagnation Temperature')

        #--------------------High Pressure Turbine inlet state variables (station 4.1)-------------------------
        Pt41 = Variable('P_{t_{4.1}}', 'kPa', 'Stagnation Pressure at the Turbine Inlet (4.1)')
        Tt41 = Variable('T_{t_{4.1}}', 'K', 'Stagnation Temperature at the Turbine Inlet (4.1)')
        ht41 = Variable('h_{t_{4.1}}', 'J/kg', 'Stagnation Enthalpy at the Turbine Inlet (4.1)')
        u41 = Variable('u_{4.1}', 'm/s', 'Flow Velocity at Station 4.1')
        T41 = Variable('T_{4.1}', 'K', 'Static Temperature at the Turbine Inlet (4.1)')

        #------------------------Variables for cooling flow model---------------------------
        # define the (f+1) variable, limits the number of signomials
        # for station 4a
        u4a = Variable('u_{4a}', 'm/s', 'Flow Velocity at Station 4a')
        M4a = Variable('M_{4a}', .1025, '-', 'User Specified Station 4a Mach #')
        P4a = Variable('P_{4a}', 'kPa', 'Static Pressure at Station 4a (4a)')
        uc = Variable('u_c', 'm/s', 'Cooling Airflow Speed at Station 4a')

        #---------------------------fuel flow fraction f--------------------------------
        f = Variable('f', '-', 'Fuel Air Mass Flow Fraction')
        fp1 = Variable('fp1', '-', 'f + 1')

        #make the constraints
        constraints = []

        with SignomialsEnabled():
            #combustor constraints
            constraints.extend([
                #flow through combustor
                ht4 == self.combustor['C_{p_{c}}'] * Tt4,

                #compute the station 4.1 enthalpy
                ht41 == self.combustor['C_{p_{c}}'] * Tt41,

                # making f+1 GP compatible --> needed for convergence
                SignomialEquality(fp1,f+1),
                # Tight([fp1 <= f+1]),
                # Relaxation of this SE makes problem hit iteration limit
                ])

            #mixing constraints
            if mixing:
                constraints.extend([
                    fp1*u41 == (u4a*(fp1)*self.combustor['\\alpha_c']*uc)**.5,
                    #this is a stagnation relation, loosened SigEq
                    # SignomialEquality(T41, Tt41-.5*(u41**2)/self.combustor['C_{p_{c}}']),
                    T41 <= Tt41-.5*(u41**2)/self.combustor['C_{p_{c}}'],

                    #here we assume no pressure loss in mixing so P41=P4a
                    Pt41 == P4a*(Tt41/T41)**(self.combustor.ccexp1),

                    #compute station 4a quantities, assumes a gamma value of 1.313 (air @ 1400K)
                    u4a == M4a*((1.313*R*Tt4)**.5)/self.combustor['hold_{4a}'],
                    uc == self.combustor['r_{uc}']*u4a,
                    P4a == Pt4*self.combustor['hold_{4a}']**(self.combustor.ccexp2),
                    ])
            #combustor constraints with no mixing
            else:
                constraints.extend([
                    Pt41 == Pt4,
                    Tt41 == Tt4,
                    ])

        return constraints
