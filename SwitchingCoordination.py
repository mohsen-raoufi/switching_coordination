import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx

def InitParams(N=3,couplingStrength=1.0,noiseStd=0.1,switchingRate=1.0,
               refTime=1.0,dt=0.1,simTime=100.0,outTime=1.0,avgFrequency=0.0,stdFrequency=0.0,writeFile=False,showAnimation=False):

    ''' Initialize parameter dictionary'''

    params=dict()
    params['N']=N # number of agents
    params['couplingStrength']=couplingStrength # alignment strength
    params['noiseStd']=noiseStd
    params['noiseAmplitude']=noiseStd*np.sqrt(dt) # noise strength 
    params['switchingRate']=switchingRate # rate of randomly switching neighbors
    params['refTime']=refTime # refractory time -> "blind time" after switching
    params['avgFrequency']=avgFrequency # for Kuramoto model -> average eigenfrequency of agents
    params['stdFrequency']=stdFrequency # for Kuramoto model -> std. deviation of eigenfrequencies
    params['dt']=dt # numerical time step
    params['simTime']=simTime # simulation time in time units
    params['simSteps']=int(simTime/dt) # number simulation steps
    params['outTime']=outTime; # time intervall between outputs saved to outData
    params['outStep']=int(outTime/dt) # time steps between outputs to outData
    params['writeFile']=writeFile # write results to file
    params['showAnimation']=showAnimation # bool: show animation of the graph + phase of the agents

    return params


def InitData(params):
    ''' Initialize dictionaries keeping the simulation data'''
    
    data=dict()
    
    # phase/polar angle or heading of the agent
    data['phi']=2*np.pi*np.random.random(params['N'])
    # set Kuramoto frequencies -> if avg and std set to 0 then all omega=0 -> XY model, simple directional alignment
    if(params['avgFrequency']!=0 or params['stdFrequency']!=0):
        data['omega']=np.random.normal(loc=params['avgFrequency'],scale=params['stdFrequency'],size=params['N']) 
    else:
        data['omega']=np.zeros(params['N'])

    # index of the neighbor which the agent pays attention to,  each agent has only one neighbor 
    data['neighbor']=np.int32(np.zeros(params['N']))
    # initial random assignment of a neighbor
    nArray=np.arange(params['N'])
    for i in range(params['N']):
        data['neighbor'][i]=np.random.choice(np.delete(nArray,i))
    
    # timer keeping track of refractory time. If timer>0 refractory phase -> no coupling.
    data['timer']=-1*np.ones(params['N'])
    # array keeping track of coupling strength: If timer>0 coupling=0, else coupling=params['couplingStrength']
    data['coupling']=np.ones(params['N'])


    # initialize outData dict for saving only every n-th time step 
    outData=dict()
    outData['t']=[0]
    outData['phi']=[2*np.pi*np.random.random(params['N'])]
    outData['neighbor']=[data['neighbor']]
    outData['timer']=[-params['dt']*np.ones(params['N'])]
    outData['order']=[np.abs(np.mean(np.exp(1j * data['phi'])))]

    return data,outData

def UpdatePhiTimer(phi,omega,timer,neighbor,coupling,N,dt,noiseAmplitude):
    ''' update the phi according to the Euler scheme, update timer'''
    noise = np.random.normal(loc=0,scale=noiseAmplitude,size=N)
    dphi  = (omega+coupling*np.sin(phi[neighbor] - phi))*dt + noise

    phi+= dphi
    #print('--------------')
    #print(timer)
    timer-=dt
    #print(timer)
    #print('--------------')
    phi=np.mod(phi,2*np.pi) 

    return phi,timer

def UpdateNetwork(neighbor,timer,coupling,switchingRate,dt,N,refTime,couplingStrength):
    ''' update network, rewire with probability switchingRate*dt per time step per agent'''
    rnd=np.random.random(size=N)
    switchArray=rnd<(switchingRate*dt)

    #print(rnd,switchArray)
    #print(neighbor)
    
    nArray=np.arange(N)
    for idx in np.where(switchArray)[0]:
        # neighbor[idx]=np.random.choice(np.delete(nArray,np.int32([neighbor[idx],idx])))
        neighbor[idx]=np.random.choice(nArray)#,np.int32([neighbor[idx],idx])))

    timer[switchArray]=refTime
    coupling[:]=couplingStrength
    coupling[timer>0]=0.0;

    return neighbor,timer,coupling

def UpdateOutData(params,data,outData,t):
    ''' append current time step results to outData '''
    outData['t'].append(t)
    outData['phi'].append(np.copy(data['phi']))
    outData['neighbor'].append(np.copy(data['neighbor'][:]))
    outData['timer'].append(np.copy(data['timer'][:]))
    outData['order'].append(np.abs(np.mean(np.exp(1j * data['phi']))))

    return 

def GenerateOutputString(params):
    ''' generate output string for file name'''
    nameString='N'+str(params['N'])+'-K'+str(params['couplingStrength'])+'-R'+str(params['switchingRate'])+'-sigma'+str(params['noiseStd'])
    nameString=nameString.replace('.','_')

    return nameString

def SaveResultsToFile(params,outData):
    ''' save order parameter to file '''
    saveDict={k:v for k,v in outData.items() if k in ['t','order'] }
    df=pd.DataFrame.from_dict(saveDict)
    nameString=GenerateOutputString(params)
    df.to_csv('results/order-'+nameString+'.csv',index=False)
    
    return

def make_network(params, data):
    ''' return the network given the list of neighbors '''
    # make an empty directed graph with N nodes
    G = nx.empty_graph(n=params['N'],create_using=nx.DiGraph())

    # go through the neighbor list and add the corresponding directed edge
    for node_i in range(G.number_of_nodes()):
        G.add_edge(*(node_i, data['neighbor'][node_i]))
    
    return G

def draw_animation_frame(params, data, time, fig, ax, node_pos):
    ''' return a figure showing the agents on a graph + their states as color '''
    
    G = make_network(params=params, data=data)

    # Background nodes
    # nx.draw_networkx_edges(G, ax=ax, edge_color="gray")
    # null_nodes = nx.draw_networkx_nodes(G, pos=pos, nodelist=set(G.nodes()) - set(path), node_color="white",  ax=ax)
    # null_nodes.set_edgecolor("black")

    # Query nodes
    # query_nodes = nx.draw_networkx_nodes(G, pos=pos, nodelist=path, node_color=idx_colors[:len(path)], ax=ax)
    # query_nodes.set_edgecolor("white")
    # nx.draw_networkx_labels(G, pos=pos, labels=dict(zip(path,path)),  font_color="white", ax=ax)
    # edgelist = [path[k:k+2] for k in range(len(path) - 1)]
    # nx.draw_networkx_edges(G, pos=pos, edgelist=edgelist, width=idx_weights[:len(path)], ax=ax)


    ax.clear()

    nx.draw(G, with_labels=True, pos=node_pos) #, ax=ax)
    ax.set_title("Frame %d:  order:   "%(time+1))
    plt.pause(0.0001)
    
    # ax.set_title("Frame %d:    "%(time+1))
    # ax.set_xticks([])
    # ax.set_yticks([])

    # return ax


def SingleSimulation(params,data=[]):
    ''' perform a single run'''

    #initialize data if not passed to function
    if(len(data)==0):
        data,outData=InitData(params)

    # make fig and ax for animation visualization
    if(params['showAnimation']):
        fig, ax = plt.subplots();
        G_null = nx.empty_graph(n=params['N'],create_using=nx.DiGraph())
        node_pos = nx.circular_layout(G_null)


    # perform time loop for simple Euler scheme integration
    for t in range(1,params['simSteps']):

        data['phi'], data['timer'] = UpdatePhiTimer(data['phi'],data['omega'],data['timer'],data['neighbor'],data['coupling'],
                                                    params['N'],params['dt'],params['noiseAmplitude'])

        data['neighbor'],data['timer'],data['coupling'] = UpdateNetwork(data['neighbor'],data['timer'],data['coupling'],
                                                                        params['switchingRate'],params['dt'],params['N'],
                                                                        params['refTime'],params['couplingStrength'])
        
        if(params['showAnimation']):
            draw_animation_frame(params=params, data=data, time=t, ax=ax, fig=fig, node_pos=node_pos)
        

        #write outData 
        if(t % params['outStep']==0):
            UpdateOutData(params,data,outData,t*params['dt'])


    # save results to file
    if(params['writeFile']):
        SaveResultsToFile(params,outData)
        

    return outData, data



