from igraph import *
import random


# Infected get recover with probabiliti teta
# An infected node attemps to infect each neighbor with probability β
# the updates at time depend only on the state at t-1


### PARAMETERS ###
teta = 0.2 # Teta -> virus death rate
beta = 0.3 # Beta -> Infection Rate
p0 = 5 # inital number of infected nodes
tmax = 30 # between 20 and 50
network_size = 1000
p0 = 5 # inital number of infected nodes
##################

g = Graph.Barabasi(n=1000,m=400, directed=False)





def init_state(n):
    initial_state = dict() # key-mapped data structure
    
    initial_state.update(n.nodes()) # each node is a key
    initial_state.update()
    return initial_state
    
    

def simulation(g: Graph, p0: int, teta: float, beta: float, tmax: int):

    #Get indices of nodes
    vertices = g.vs.indices

    #Initialize all the node as NOT infected, creating a new attribute
    for i in vertices:
        g.vs[i]["Infected"] = 0

    #Set p0 initial nodes to infected
    infected =random.sample(vertices,p0)
    for i in infected:
        g.vs[i]["Infected"] = 1

    #Get list of indexes of not infected nodes
    not_infected = list(set(vertices) - set(infected))


    for t in range(0,tmax):
        pass
      #Perform Bernoulli (maybe binomial since we are using a list and not a single node), then based on the outcome update the infected/not infected lists
      #DO it on Parallel








def get_infected(network):
    pass

if __name__ == "__main__":
    net = nx.barabasi_albert_graph(network_size,8)
    init_state(net)