from gpkit import Variable, Model

class FanMap(Model):
    """"
    Fan map model
    """
    def setup(self):
        #define new variables
        #------------------Fan map variables----------------
        mFanBarD = Variable('\\bar{m}_{fan_{D}}', 'kg/s', 'Fan On-Design Corrected Mass Flow')
        piFanD = Variable('\\pi_{f_D}', '-', 'On-Design Pressure Ratio')

    def dynamic(self, engine):
        """
        creates an instance of the fan map performance model
        """
        return FanMapPerformance(self, engine)

class FanMapPerformance(Model):
    """
    Fan map perfomrance constraints
    """
    def setup(self, fanmap, engine):
        self.fanmap = fanmap
        self.engine = engine

        #define new variables
        #-----------------------Fan Map Variables--------------------
        #Mass Flow Variables
        mf = Variable('m_{f}', 'kg/s', 'Fan Corrected Mass Flow')
        mtildf = Variable('m_{tild_f}', '-', 'Fan Normalized Mass Flow')

        #----------------------Compressor Speeds--------------------
        #Speed Variables...by setting the design speed to be 1 since only ratios are
        #imporant I was able to drop out all the other speeds
        Nf = Variable('N_{f}', '-', 'Fan Speed')

        #make the constraints
        constraints = []

        #fan map
        constraints.extend([
            #define mtild
            mtildf == mf/self.fanmap['\\bar{m}_{fan_{D}}'],   #B.282
            ])

        return constraints

class LPCMap(Model):
    """"
    LPC map model
    """
    def setup(self):
        #define new variables
        #-----------------LPC map variables-------------------
        mlcD = Variable('m_{lc_D}', 'kg/s', 'On Design LPC Corrected Mass Flow')
        pilcD = Variable('\pi_{lc_D}', '-', 'LPC On-Design Pressure Ratio')

    def dynamic(self, engine):
        """
        creates an instance of the HPC map performance model
        """
        return LPCMapPerformance(self, engine)

class LPCMapPerformance(Model):
    """
    LPC map perfomrance constraints
    """
    def setup(self, lpcmap, engine):
        self.lpcmap = lpcmap
        self.engine = engine

        #define new variables
        #-------------------------LPC Map Variables-------------------------
        #Mass Flow Variables
        mlc = Variable('m_{lc}', 'kg/s', 'LPC Corrected Mass Flow')
        mtildlc = Variable('m_{tild_lc}', '-', 'LPC Normalized Mass Flow')

        #----------------------Compressor Speeds--------------------
        #Speed Variables...by setting the design speed to be 1 since only ratios are
        #imporant I was able to drop out all the other speeds
        N1 = Variable('N_{1}', '-', 'LPC Speed')

        #make the constraints
        constraints = []

        #LPC map
        constraints.extend([
            #define mtild
            mtildlc == mlc/self.lpcmap['m_{lc_D}'],   #B.282
        ])

        return constraints

class HPCMap(Model):
    """"
    HPC map model
    """
    def setup(self):
        #define new variables
        #-----------------HPC map variables-----------
        mhcD = Variable('m_{hc_D}', 'kg/s', 'On Design HPC Corrected Mass Flow')
        pihcD = Variable('\\pi_{hc_D}', '-', 'HPC On-Design Pressure Ratio')
        mhcD = Variable('m_{hc_D}', 'kg/s', 'On Design HPC Corrected Mass Flow')

    def dynamic(self, engine):
        """
        creates an instance of the HPC map performance model
        """
        return HPCMapPerformance(self, engine)

class HPCMapPerformance(Model):
    """
    HPC map perfomrance constraints
    """
    def setup(self, hpcmap, engine):
        self.hpcmap = hpcmap
        self.engine = engine

        #define new variables
        #--------------------------HPC Map Variables------------------
        #Mass Flow Variables
        mhc = Variable('m_{hc}', 'kg/s', 'HPC Corrected Mass Flow')
        mtildhc = Variable('m_{tild_{hc}}', '-', 'HPC Normalized Mass Flow')

        #----------------------Compressor Speeds--------------------
        #Speed Variables...by setting the design speed to be 1 since only ratios are
        #imporant I was able to drop out all the other speeds
        N2 = Variable('N_{2}', '-', 'HPC Speed')

        #make the constraints
        constraints = []

        #HPC map
        constraints.extend([
            #define mtild
            mtildhc == mhc/self.hpcmap['m_{hc_D}'],   #B.282
            ])

        return constraints
