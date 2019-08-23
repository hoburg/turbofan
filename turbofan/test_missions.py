from gpkit import Variable, Model, units

# Test missions specify thrust and ambient conditions
# at different flight states, for different engines.

class TestMissionCFM(Model):
    def setup(self, engine):
        self.engine = engine
        M2 = .6
        M25 = .6
        M4a = .1025
        M0 = .8

        climb = [
            engine['F_{spec}'][0] == 5496.4 * 4.4 * units('N'),
            engine['F_{spec}'][1] == 5961.9 * 4.4 * units('N'),

            engine.state["T_{atm}"] == 218*units('K'),


            engine.state['P_{atm}'][0] == 23.84*units('kPa'),    #36K feet
            engine.state['M'][0] == M0,
            engine.engineP['M_2'][0] == M2,
            engine.engineP['M_{2.5}'][0] == M25,
            engine.engineP['hold_{2}'] == 1+.5*(1.398-1)*M2**2,
            engine.engineP['hold_{2.5}'] == 1+.5*(1.354-1)*M25**2,
            engine.engineP['c1'] == 1+.5*(.401)*M0**2,
            ]

        M2 = .6
        M25 = .6
        M4a = .1025
        M0 = .8

        cruise = [
            engine.state['P_{atm}'][1] == 23.84*units('kPa'),    #36K feet
            engine.state['M'][1] == M0,
            engine.engineP['M_2'][1] == M2,
            engine.engineP['M_{2.5}'][1] == M25,

            engine['W_{engine}'] <= 5216*units('lbf'),
            ]

        return climb, cruise

class TestMissionTASOPT(Model):
    def setup(self, engine):
        self.engine = engine
        M2 = .6
        M25 = .6
        M4a = .1025
        M0 = .8025

        toclimb = [
            engine['F_{spec}'][0] == 94.971 * units('kN'),
            engine['F_{spec}'][1] == 30.109* units('kN'),
            engine['F_{spec}'][2] == 22.182 * units('kN'),

            engine.state['P_{atm}'][1] == 23.84*units('kPa'),    #36K feet
            engine.state["T_{atm}"][1] == 218*units('K'),
            engine.state['M'][1] == M0,
            engine.engineP['M_2'][1] == M2,
            engine.engineP['M_{2.5}'][1] == M25,
            engine.engineP['hold_{2}'][1] == 1+.5*(1.398-1)*M2**2,
            engine.engineP['hold_{2.5}'][1] == 1+.5*(1.354-1)*M25**2,
            engine.engineP['c1'][1] == 1+.5*(.401)*M0**2,
            ]

        M2 = .6
        M25 = .6
        M4a = .1025
        M0 = .8

        cruise = [
            engine.state['P_{atm}'][2] == 23.92*units('kPa'),    #36K feet
            engine.state["T_{atm}"][2] == 219.4*units('K'),
            engine.state['M'][2] == M0,
            engine.engineP['M_2'][2] == M2,
            engine.engineP['M_{2.5}'][2] == M25,
            engine.engineP['hold_{2}'][2] == 1+.5*(1.398-1)*M2**2,
            engine.engineP['hold_{2.5}'][2] == 1+.5*(1.354-1)*M25**2,
            engine.engineP['c1'][2] == 1+.5*(.401)*M0**2,
            ]

        M2 = .6
        M25 = .6
        M4a = .1025
        M0 = .2201

        rotation = [
            engine.state['P_{atm}'][0] == 101.325*units('kPa'),
            engine.state["T_{atm}"][0] == 288*units('K'),
            engine.state['M'][0] == M0,
            engine.engineP['M_2'][0] == M2,
            engine.engineP['M_{2.5}'][0] == M25,
            engine.engineP['hold_{2}'][0] == 1+.5*(1.401-1)*M2**2,
            engine.engineP['hold_{2.5}'][0] == 1+.5*(1.401-1)*M25**2,
            engine.engineP['c1'][0] == 1+.5*(.401)*M0**2,

            engine['W_{engine}'] <= 1.1*7870.7*units('lbf'),
            ]

        return rotation, toclimb, cruise

class TestMissionGE90(Model):
    def setup(self, engine):
        self.engine = engine
        M2 = .65
        M25 = .6
        M4a = .1025
        M0 = .85

        climb = [
            engine['F_{spec}'][1] == 19600.4*units('lbf'),
            engine['F_{spec}'][0] == 16408.4 * units('lbf'),

            engine.state['P_{atm}'][1] == 23.84*units('kPa'),    #36K feet
            engine.state["T_{atm}"][1] == 218*units('K'),
            engine.state['M'][1] == M0,
            engine.engineP['M_2'][1] == M2,
            engine.engineP['M_{2.5}'][1] == M25,
            engine.engineP['hold_{2}'][1] == 1+.5*(1.398-1)*M2**2,
            engine.engineP['hold_{2.5}'][1] == 1+.5*(1.354-1)*M25**2,
            engine.engineP['c1'][1] == 1+.5*(.401)*M0**2,
            ]

        M0 = .65
        M2 = .6
        M25 = .45
        M4a = .1025

        cruise = [
            engine.state['P_{atm}'][0] == 23.84*units('kPa'),    #36K feet
            engine.state["T_{atm}"][0] == 218*units('K'),
            engine.state['M'][0] == M0,
            engine.engineP['M_2'][0] == M2,
            engine.engineP['M_{2.5}'][0] == M25,
            engine.engineP['hold_{2}'][0] == 1+.5*(1.398-1)*M2**2,
            engine.engineP['hold_{2.5}'][0] == 1+.5*(1.354-1)*M25**2,
            engine.engineP['c1'][0] == 1+.5*(.401)*M0**2,

            engine['W_{engine}'] <= 77399*units('N'),
            ]

        return climb, cruise

class TestMissionD82(Model):
    def setup(self, engine):
        self.engine = engine
        M2 = .6
        M25 = .6
        M4a = .1025
        M0 = .5

        engineclimb = [
            engine['F_{spec}'][1] == 13.798*units('kN'),
            engine['F_{spec}'][0] == 11.949 * units('kN'),

            engine.state['P_{atm}'][1] == 23.84*units('kPa'),    #36K feet
            engine.state["T_{atm}"][1] == 218*units('K'),
            engine.state['M'][1] == M0,
            engine.engineP['M_{2.5}'][1] == M25,
            engine.engineP['hold_{2}'][1] == 1+.5*(1.398-1)*M2**2,
            engine.engineP['hold_{2.5}'][1] == 1+.5*(1.354-1)*M25**2,
            engine.engineP['c1'][1] == 1+.5*(.401)*M0**2,
            ]

        M2 = .6
        M25 = .6
        M4a = .1025
        M0 = .72

        enginecruise = [
            engine.state['P_{atm}'][0] == 23.84*units('kPa'),    #36K feet
            engine.state["T_{atm}"][0] == 218*units('K'),
            engine.state['M'][0] == M0,
            engine.engineP['M_2'][0] == M2,
            engine.engineP['M_{2.5}'][0] == M25,
            engine.engineP['hold_{2}'][0] == 1+.5*(1.398-1)*M2**2,
            engine.engineP['hold_{2.5}'][0] == 1+.5*(1.354-1)*M25**2,
            engine.engineP['c1'][0] == 1+.5*(.401)*M0**2,

            engine['W_{engine}'] <= 10502.8*units('lbf'),
            ]


        return engineclimb, enginecruise
