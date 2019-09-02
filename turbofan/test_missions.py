from gpkit import Variable, Model, units
from gpkit.small_scripts import mag

# Test missions specify thrust and ambient conditions
# at different flight states, for different engines.

class TestMissionCFM(Model):
    def setup(self, engine):
        self.engine = engine
        M2 = .6
        M25 = .6
        M0 = .8

        climb = [
            engine['F'][0] == 5496.4 * 4.4 * units('N'),
            engine['F'][1] == 5961.9 * 4.4 * units('N'),

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
        M0 = .8025

        toclimb = [
            engine['F'][0] == 94.971 * units('kN'),
            engine['F'][1] == 30.109* units('kN'),
            engine['F'][2] == 22.182 * units('kN'),

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
        M0 = .85

        climb = [
            engine['F'][1] == 19600.4*units('lbf'),
            engine['F'][0] == 16408.4 * units('lbf'),

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
        M0 = .5

        engineclimb = [
            engine['F'][1] == 13.798*units('kN'),
            engine['F'][0] == 11.949 * units('kN'),

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

def diffs(sol, eng):
    if eng == 0:
        tocerror = 100*(mag(sol('TSFC')[1]) - .6941)/.6941
        cruiseerror = 100*(mag(sol('TSFC')[0]) - .6793)/.6793
        weighterror =  100*(mag(sol('W_{engine}').to('lbf'))-5216)/5216

        print("Cruise error")
        print(cruiseerror)
        print("TOC error")
        print(tocerror)
        print("Weight Error")
        print(weighterror)

    if eng == 1:
        rotationerror = 100*(mag(sol('TSFC')[0]) - .48434)/.48434
        tocerror = 100*(mag(sol('TSFC')[1]) - .65290)/.65290
        cruiseerror = 100*(mag(sol('TSFC')[2]) - .64009)/.64009

        print(rotationerror, tocerror, cruiseerror)

##        print 100*(mag(sol('A_{2}').to('m^2'))-1.6026)/1.6026
##        print 100*(mag(sol('A_{7}').to('m^2'))-.7423)/.7423
##        print 100*(mag(sol('A_{5}').to('m^2'))-.2262)/.2262
        print("----weight---")
        print(100*(mag(sol('W_{engine}').to('lbf'))-7870.7)/7870.7)

        print("TO Tt4.1")
        print(100*(mag(sol('T_{t_{4.1}}')[0]) - 1658.7)/1658.7)
        print("TOC Tt4.1")
        print(100*(mag(sol('T_{t_{4.1}}')[1]) - 1605.4)/1605.4)
        print("Cruise Tt4.1")
        print(100*(mag(sol('T_{t_{4.1}}')[2]) - 1433.8)/1433.8)

        print("Cooling deltas")
        print("TO")
        print(100*(mag(sol('T_{t_4}')[0]-sol('T_{t_{4.1}}')[0]) - 174.3)/174.3)
        print("TOC")
        print(100*(mag(sol('T_{t_4}')[1]-sol('T_{t_{4.1}}')[1]) - 178.4)/178.4)
        print("Cruise")
        print(100*(mag(sol('T_{t_4}')[2]-sol('T_{t_{4.1}}')[2]) - 153.2)/153.2)

    if eng == 2:
        tocerror = 100*(mag(sol('TSFC')[1]) - 0.5846)/0.5846
        cruiseerror = 100*(mag(sol('TSFC')[0]) - 0.5418)/0.5418
        weighterror =  100*(mag(sol('W_{engine}').to('lbf'))-17400)/17400

        print(tocerror, cruiseerror, weighterror)

    if eng == 3:
        print("Sorry, no data available for the D8.2.")


#code for estimating on design parameters
##Pt0 = 50
##Tt0 = 250
##Pt3 = Pt0*lpc*fan*hpc
##Pt21 = fan * Pt0
##Pt25 = Pt0 * fan * lpc
##
##Tt21 = Tt0 * (fan)**(.4/(1.4*.9153))
##Tt25 = Tt21 * (lpc)**(.4/(1.4*.9037))
##Tt3 = Tt25 * (hpc)**(.4/(1.4*.9247))
##
##Tt41 = 1400
##
##Tt45 = Tt41 - (Tt3 - Tt25)
##
##Tt49 = Tt45 - (Tt25 - Tt21)
##
##piHPT = (Tt45/Tt41)**(.9121*1.4/.4)
##
##piLPT = (Tt49/Tt45)**(.9228*1.4/.4)
##
##Pt45 = piHPT * Pt3
