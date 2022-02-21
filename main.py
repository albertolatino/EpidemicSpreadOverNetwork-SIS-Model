from igraph import *
from functions import compute_threshold, simulation, saveGraph
import time
import numpy as np
# Infected get recover with probabiliti teta
# An infected node attemps to infect each neighbor with probability β
# the updates at time depend only on the state at t-1

# Beta -> Infection Rate
# Teta -> virus death rate

### PARAMETERS ###
gamma = 0.2
beta = 0.4
tmax = 100 # between 20 and 50
network_size = 1000
p0 = 100 # inital number of infected nodes
epsilon = 0.5
experiments = 10
##################

### Setting ###
make_video = False # disable to save time and avoid video
###############

#Fully Connected Network
g1 = Graph.Full(n=network_size)
#Scale-Free Barabasi
g2 = Graph.Barabasi(n=network_size,m=5, directed=False)
#Watts Strogatz
g3 = Graph.Watts_Strogatz(dim = 1, size=network_size, nei = 1 , p = 0.5 )#E' GIA UNA SMALL WORLD????
#Tree Network
g4 = Graph.Tree(network_size, 2)
#Star Network
g5 = Graph.Star(network_size)
#Regular Lattice
g6 = Graph.Lattice(dim=[31, 31], circular=False)#1000 nodes is ok 31x31???
#Erdos Renyi
g7 = Graph.Erdos_Renyi(network_size, p=0.5, directed=False, loops=False)

nets = [
        {"name": "FullGraph","net":g1},
        {"name": "Barabasi","net":g2},
        {"name": "Watts_Strogatz","net":g3},
        {"name": "Tree","net":g4},
        {"name": "Star","net":g5},
        {"name": "Lattice","net":g6},
        {"name": "Erdos_Renyi","net":g7}
        ]

start = time.time()

print("Start Simulation of the nets")

if(make_video):
    print("Video enabled, this option can slow the entire process")

for net in nets:
    name = net["name"]
    print(f"Network : {name}")
    g = net["net"]

    th = np.real(compute_threshold(g))
    bg = round(beta/gamma,4)

    beta_min = gamma * (th) * (1 - epsilon)
    beta_plus = gamma * (th) * (1 + epsilon)
    beta_equal = gamma *(th)

    equal_avg = []
    above_avg = []
    below_avg = []
    for i in range(1,experiments):

     infected_equal = simulation(g,p0,gamma,beta_equal,tmax,name = name,video=make_video)
     infected_below = simulation(g, p0, gamma, beta_min, tmax, name=name, video=make_video)
     infected_above = simulation(g, p0, gamma, beta_plus, tmax, name=name, video=make_video)
     equal_avg.append(infected_equal)
     above_avg.append(infected_above)
     below_avg.append(infected_below)

    equal = np.array(equal_avg)
    below = np.array(below_avg)
    above = np.array(above_avg)

    saveGraph(tmax, np.mean(equal, axis=0), np.mean(below, axis=0),np.mean(above, axis=0) , gamma, beta, name)





    print(f"Done in {round(time.time()-start,4)}")

for net in nets:
    name = net["name"]
    print(f"Network : {name}")
    g = net["net"]
    gamma = 0.4
    beta = 0.1
    th = np.real(compute_threshold(g))
    bg = round(beta / gamma, 4)

    infected_equal = simulation(g, p0, gamma, beta, tmax, name=name, video=make_video)
    saveGraph(tmax, infected_equal,None, None, gamma, beta, name)

    print(f"    1/λ : {th}, β/γ {bg}")
    if (bg > th):
        print(" Epidemic occur β/γ > 1/λ ")
    else:
        print(" Epidemic not occur β/γ < 1/λ")


