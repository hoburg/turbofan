from gpkit import Variable, Model

class Compressor(Model):
    """"
    Compressor model
    """
    def setup(self, fexp1, lpcexp1, hpcexp1):
        self.fexp1 = fexp1
        self.lpcexp1 = lpcexp1
        self.hpcexp1 = hpcexp1
        #define new variables
        # fan, LPC, HPC
        Cp1 = Variable('C_{p_{1}', 1008, 'J/kg/K', "Cp Value for Air at 350K")#gamma = 1.398
        Cp2 = Variable('C_{p_{2}', 1099, 'J/kg/K', "Cp Value for Air at 800K") #gamma = 1.354

        #-------------------------diffuser pressure ratios--------------------------
        pid = Variable('\pi_{d}', '-', 'Diffuser Pressure Ratio')
        pifn = Variable('\\pi_{fn}', '-', 'Fan Duct Pressure Loss Ratio')

        gammaAir = Variable('\\gamma_{air}', 1.4, '-', 'Specific Heat Ratio for Ambient Air')
        Cpair = Variable('C_{p_{air}', 1003, 'J/kg/K', "Cp Value for Air at 250K")

    def dynamic(self, engine, state, BLI):
        """
        creates an instance of the compressor performance model
        """
        return CompressorPerformance(self, engine, state, BLI)

class CompressorPerformance(Model):
    """
    compressor performance constraints
    """
    def setup(self, compressor, engine, state, BLI):
        self.compressor = compressor
        self.engine = engine

        #define new variables
        #--------------------------free stream stagnation states--------------------------
        Pt0 = Variable('P_{t_0}', 'kPa', 'Free Stream Stagnation Pressure')
        Tt0 = Variable('T_{t_0}', 'K', 'Free Stream Stagnation Temperature')
        ht0 = Variable('h_{t_0}', 'J/kg', 'Free Stream Stagnation Enthalpy')

        #--------------------------diffuser exit stagnation states------------------------
        Pt18 = Variable('P_{t_{1.8}}', 'kPa', 'Stagnation Pressure at the Diffuser Exit (1.8)')
        Tt18 = Variable('T_{t_{1.8}}', 'K', 'Stagnation Temperature at the Diffuser Exit (1.8)')
        ht18 = Variable('h_{t_{1.8}}', 'J/kg', 'Stagnation Enthalpy at the Diffuser Exit (1.8)')

        #--------------------------fan inlet (station 2) stagnation states---------------------------
        Pt2 = Variable('P_{t_{2}}', 'kPa', 'Stagnation Pressure at the Fan Inlet (2)')
        Tt2 = Variable('T_{t_{2}}', 'K', 'Stagnation Temperature at the Fan Inlet (2)')
        ht2 = Variable('h_{t_{2}}', 'J/kg', 'Stagnation Enthalpy at the Fan Inlet (2)')

        #--------------------------fan exit (station 2.1) stagnation states---------------------
        Pt21 = Variable('P_{t_{2.1}}', 'kPa', 'Stagnation Pressure at the Fan Exit (2.1)')
        Tt21 = Variable('T_{t_{2.1}}', 'K', 'Stagnation Temperature at the Fan Exit (2.1)')
        ht21 = Variable('h_{t_{2.1}}', 'J/kg', 'Stagnation Enthalpy at the Fan Exit (2.1)')

        #-------------------------LPC exit (station 2.5) stagnation states-------------------
        Pt25 = Variable('P_{t_{2.5}}', 'kPa', 'Stagnation Pressure at the LPC Exit (2.5)')
        Tt25 = Variable('T_{t_{2.5}}', 'K', 'Stagnation Temperature at the LPC Exit (2.5)')
        ht25 = Variable('h_{t_{2.5}}', 'J/kg', 'Stagnation Enthalpy at the LPC Exit (2.5)')

        #--------------------------HPC exit stagnation states (station 3)---------------------
        Pt3 = Variable('P_{t_3}', 'kPa', 'Stagnation Pressure at the HPC Exit (3)')
        Tt3 = Variable('T_{t_3}', 'K', 'Stagnation Temperature at the HPC Exit (3)')
        ht3 = Variable('h_{t_3}', 'J/kg', 'Stagnation Enthalpy at the HPC Exit (3)')

        #---------------------------fan nozzle exit (station 7) stagnation states---------------
        Pt7 = Variable('P_{t_7}', 'kPa', 'Stagnation Pressure at the Fan Nozzle Exit (7)')
        Tt7 = Variable('T_{t_7}', 'K', 'Stagnation Temperature at the Fan Nozzle Exit (7)')
        ht7 = Variable('h_{t_7}', 'J/kg', 'Stagnation Enthalpy at the Fan Nozzle Exit (7)')

        #------------------------turbo machinery pressure ratios--------------
        pif = Variable('\\pi_{f}', '-', 'Fan Pressure Ratio')
        pilc = Variable('\\pi_{lc}', '-', 'LPC Pressure Ratio')
        pihc = Variable('\\pi_{hc}', '-', 'HPC Pressure Ratio')

        hold25 = Variable('hold_{2.5}', '-', '1+(gamma-1)/2 * M_2.5**2')
        hold2 = Variable('hold_{2}', '-', '1+(gamma-1)/2 * M_2**2')
        c1 = Variable('c1', '-', 'Constant in Stagnation Eqn')

        diffuser = [
            #free stream stagnation values
             #https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
            Tt0 == state["T_{atm}"] / (c1) ** (-1),
            ht0 == self.compressor['C_{p_{air}'] * Tt0,

            #diffuser exit stagnation values (station 1.8)
            Pt18 == self.compressor['\pi_{d}'] * Pt0,  #B.113
            Tt18 == Tt0,        #B.114
            ht18 == ht0,        #B.115
            ]

        if BLI:
            diffuser.extend([
                Pt0 == self.engine['f_{BLI_{P}}']*state["P_{atm}"] / (c1 ** -3.5),
                ])

        if not BLI:
            diffuser.extend([
                Pt0 == state["P_{atm}"] / (c1 ** -3.5),
                ])

        fan = [
            #fan inlet constraints (station 2)
            Tt2 == Tt18,    #B.120
            ht2 == ht18,    #B.121
            Pt2 == Pt18,

            #fan exit constraints (station 2.1)
            Pt21 == pif * Pt2,  #16.50
            Tt21 == Tt2 * pif ** (self.compressor.fexp1),   #16.50
            ht21 == self.compressor['C_{p_{air}'] * Tt21,   #16.50

            #fan nozzle exit (station 7)
            Pt7 == self.compressor['\\pi_{fn}'] * Pt21,     #B.125
            Tt7 == Tt21,    #B.126
            ht7 == ht21,    #B.127
            ]

        lpc = [
            #LPC exit (station 2.5)
            Pt25 == pilc * pif * Pt2,
            Tt25 == Tt2 * (pif*pilc) ** (self.compressor.lpcexp1),
            ht25 == Tt25 * self.compressor['C_{p_{1}'],
            ]

        hpc = [
            Pt3 == pihc * Pt25,
            Tt3 == Tt25 * pihc ** (self.compressor.hpcexp1),
            ht3 == self.compressor['C_{p_{2}'] * Tt3
            ]

        return diffuser, fan, lpc, hpc
