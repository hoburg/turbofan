from gpkit import Variable, Model

class Turbine(Model):
    """"
    Turbine model
    """
    def setup(self, hptexp1, lptexp1):
        self.hptexp1 = hptexp1
        self.lptexp1 = lptexp1
        #define new variables
        #turbines
        Cpt1 = Variable('C_{p_{t1}}', 1280, 'J/kg/K', "Cp Value for Combustion Products in HP Turbine") #1300K gamma = 1.318
        Cpt2 = Variable('C_{p_{t2}}', 1184, 'J/kg/K', "Cp Value for Combustion Products in LP Turbine") #800K gamma = 1.354

        #-------------------------diffuser pressure ratios--------------------------
        pitn = Variable('\\pi_{tn}', '-', 'Turbine Nozzle Pressure Ratio')

        #---------------------------efficiencies & takeoffs-----------------------
        # Note: HP and LP shaft efficiencies smear losses for electrical power
        etaHPshaft = Variable('\eta_{HPshaft}', '-', 'Power Transmission Efficiency of High Pressure Shaft')
        etaLPshaft = Variable('\eta_{LPshaft}', '-', 'Power Transmission Efficiency of Low Pressure Shaft')

    def dynamic(self, engine):
        """
        creates an instance of the fan map performance model
        """
        return TurbinePerformance(self, engine)

class TurbinePerformance(Model):
    """
    combustor performance constraints
    """
    def setup(self, turbine, engine):
        self.turbine = turbine
        self.engine = engine

        #define new variables
        #------------------HPT exit stagnation states (station 4.5)------------------
        ht45 = Variable('h_{t_{4.5}}', 'J/kg', 'Stagnation Enthalpy at the HPT Exit (4.5)')
        Pt45 = Variable('P_{t_{4.5}}', 'kPa', 'Stagnation Pressure at the HPT Exit (4.5)')
        Tt45 = Variable('T_{t_{4.5}}', 'K', 'Stagnation Temperature at the HPT Exit (4.5)')

        #-------------------LPT exit (station 4.5) stagnation states--------------------
        Pt49 = Variable('P_{t_{4.9}}', 'kPa', 'Stagnation Pressure at the LPT Exit (49)')
        Tt49 = Variable('T_{t_{4.9}}', 'K', 'Stagnation Temperature at the LPT Exit (49)')
        ht49 = Variable('h_{t_{4.9}}', 'J/kg', 'Stagnation Enthalpy at the LPT Exit (49)')

        #------------------------turbo machinery pressure ratios--------------
        pihpt = Variable('\\pi_{HPT}', '-', 'HPT Pressure Ratio')
        pilpt = Variable('\\pi_{LPT}', '-', 'LPT Pressure Ratio')

        #------------------turbine nozzle exit stagnation states (station 5)------------
        Pt5 = Variable('P_{t_5}', 'kPa', 'Stagnation Pressure at the Turbine Nozzle Exit (5)')
        Tt5 = Variable('T_{t_5}', 'K', 'Stagnation Temperature at the Turbine Nozzle Exit (5)')
        ht5 = Variable('h_{t_5}', 'J/kg', 'Stagnation Enthalpy at the Turbine Nozzle Exit (5)')

        #make the constraints
        constraints = []

        #turbine constraints
        constraints.extend([
            #HPT Exit states (station 4.5)
            ht45 == self.turbine['C_{p_{t1}}'] * Tt45,

            #LPT Exit States
            Pt49 == pilpt * Pt45,
            pilpt == (Tt49/Tt45)**(self.turbine.lptexp1),    #turbine efficiency is 0.9
            ht49 == self.turbine['C_{p_{t2}}'] * Tt49,

            #turbine nozzle exit states
            Pt5 == self.turbine['\\pi_{tn}'] * Pt49, #B.167
            Tt5 == Tt49,    #B.168
            ht5 == ht49     #B.169
            ])

        return constraints
