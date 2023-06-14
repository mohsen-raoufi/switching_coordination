import pandas as pd
import SwitchingCoordination as sc


from importlib import reload
reload(sc)


# initialize a parameter dictionary
params=sc.InitParams(N=3,switchingRate=0.5,
                     refTime=5.0,noiseStd=.5,
                     avgFrequency=0.0, writeFile=False, showAnimation=True)




#perform a single simulation
outData, data = sc.SingleSimulation(params)