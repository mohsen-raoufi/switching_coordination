import SwitchingCoordination as sc

from importlib import reload
reload(sc)


# initialize a parameter dictionary
params=sc.InitParams(N=6,switchingRate=0.1,
                     refTime=0.0,noiseStd=.0,
                     avgFrequency=0.1, writeFile=False, 
                     simTime = 50,
                     showAnimation=True, saveAnimation=True)


#perform a single simulation
outData, data = sc.SingleSimulation(params)