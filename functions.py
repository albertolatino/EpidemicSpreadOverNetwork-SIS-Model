from igraph import *
import random
from scipy.stats import bernoulli
import time
import matplotlib.pyplot as plt
import cairo
import os
os.environ["IMAGEIO_FFMPEG_EXE"] = "/usr/bin/ffmpeg"
import glob
import moviepy.video.io.ImageSequenceClip
from igraph.drawing.text import TextDrawer
import numpy as np
import os

data = []

def initialization(g: Graph,p0: int,gamma: float, beta: float):
    # Get indices of nodes
    vertices = g.vs.indices
    # Initialize all the node as NOT infected, creating a new attribute
    g.vs.set_attribute_values("Infected",0)
    # Set p0 initial nodes to infected
    infected = list(random.sample(vertices, p0))
    g.vs[infected]["Infected"] = 1
    # Get list of indexes of not infected nodes
    not_infected = list(set(vertices) - set(infected))
    
    recover = bernoulli(gamma)
    infect = bernoulli(beta)
    
    return [g,not_infected,infected,recover,infect]

def infect_neighbors(node: int,neighbors: list,infect):
    # for the infected node take the neighbors and then perform bernulli
    inf = infect.rvs(len(neighbors)) # generate bernoulli samples for each of the infected
    return [neighbors[i] for i in range(0,len(inf)) if inf[i] == 1] # return id of infected


def saveVideoFrame(g,name,not_infected,infected,t):

    palette=ClusterColoringPalette(6)
    
    g.vs[not_infected]["color"] = 4
    g.vs[infected]["color"] = 5
    
    g.vs[not_infected]["label"] = "H"
    g.vs[infected]["label"] = "I"
    
    p = Plot(target = f"video/{name}_{t}.png",bbox=(1200, 1300),palette=palette,background="white")
    p.add(g,bbox=(20, 70, 1180, 1270))
    p.redraw()
    ctx = cairo.Context(p.surface)
    ctx.set_font_size(30)
    drawer = TextDrawer(ctx,f"{name} T = {t}",halign=TextDrawer.CENTER)
    drawer.draw_at(0, 40, width=1200)
    p.save()
    
    return

def deleteFrames():
    print("cancel files")
    files = glob.glob("/home/colo/Desktop/UPC-FIB/ComplexAndSocialNet/CSN-Lab7/video/*.png")
    for f in files:
        try:
            os.remove(f)
        except OSError as e:
            print("Error in deleting frames of the videos")
    pass

def makeFinalVideo(name: str = "final",n_image: int = 30):
    image_folder='/home/colo/Desktop/UPC-FIB/ComplexAndSocialNet/CSN-Lab7/video'
    fps=1
    
    image_files = [os.path.join(image_folder,img)
               for img in os.listdir(image_folder)
               if img.endswith(".png")]
    
    # order the sequence of images !!
    image_files_ordered = [None] * n_image
    
    for img in image_files:
        i = 1 + img.rindex("_")
        f = img.rindex(".")
        img_index = int(img[i:f])
        image_files_ordered[img_index] = img
    

    clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files_ordered, fps=fps)
    clip.write_videofile(f'final_video/{name}.mp4')

    ## Cancellare frame una volta fatto il video
    deleteFrames()


def simulation(g: Graph, p0: int, gamma: float, beta: float, tmax: int, name: str,video: bool):
    
    prop_infected_over_t = [] # for creating the graph
    prop_healthy_over_t = []
    
    [g,not_infected,infected,recover_bern,infect_bern] = initialization(g,p0,gamma,beta) # initialize simulation
    
    net_size = g.vcount()
    new_infected = []
    new_recovered = []

    #Simulation
    start_time = time.time()

    if(video):
        saveVideoFrame(g,name,not_infected,infected,0) # first frame of the video

    prop_healthy_over_t.append(len(not_infected) / net_size)
    prop_infected_over_t.append(len(infected) / net_size)

    for t in range(1,tmax):

        #Select Nodes that will be marked as infected at next T (based on Bernoulli trials)
        for id in infected:
            #Get neighbors
            neighbors_i = g.neighbors(id) 
            # remove nodes that are already infected
            neighbors_i = list(set(neighbors_i) - set(infected))
            # remove nodes that just became healthy at the current time
            neighbors_i = list(set(neighbors_i) - set(new_recovered))
            # infect neighbors 
            infected_neigh_i = infect_neighbors(id,neighbors_i,infect_bern)
            # add the new infected nodes 
            new_infected = list(set.union(set(new_infected),set(infected_neigh_i)))
            
            if(recover_bern.rvs() == 1):
                new_recovered.append(id)
                infected.remove(id) #Removed to avoid it could recover again
    
        ## State for t + 1
        
        ## Infected population
        
        infected = list(set(infected) - set(new_recovered)) # remove recovered population
        infected = list(set.union(set(infected),set(new_infected))) # add new infected population
        
        ## not Infected population
        not_infected = list(set(not_infected) - set(new_infected)) # remove the infected at this iteration 
        not_infected = list(set.union(set(not_infected),set(new_recovered))) # add recovered to the not_infected population

        
        new_infected = []
        new_recovered = []
        # Collecting data about the simulation at time t
        prop_infected_over_t.append(len(infected) / net_size)
        prop_healthy_over_t.append(len(not_infected) / net_size)

        if(video):
            saveVideoFrame(g,name,not_infected,infected,t)
    
        
    print(f"Simulation done is {round(time.time() - start_time,4)} sec")


    if(video):
        makeFinalVideo(name)

    ## SALVARE GRAFICO QUI
    #saveGraph(tmax,prop_infected_over_t,prop_healthy_over_t,gamma,beta,name)

    #Fai retur di infected e poi plotta nel main all together
    return prop_infected_over_t

def saveGraph(tmax,infected_equal,infected_below,infected_above,gamma,beta,name):

    x = list(range(0,tmax))
    plt.plot(x, infected_equal,label="Infected",color='b')
    if(infected_below is not None and infected_above is not None):
        plt.plot(x, infected_below, label="Infected", color = 'g')
        plt.plot(x, infected_above, label="Infected", color = 'r')
    plt.xlabel("Time")
    plt.ylabel("Infected ratio")
    plt.title("Infected Over time")
    if (infected_below is not None and infected_above is not None):
        plt.savefig(f"img/{name}_comparison.png")
    else:
        plt.savefig(f"img/{name}_{gamma}_{beta}.png")
    plt.close()

    return

def compute_threshold(g):
    #Beta/Gamma = Threshold
    #Compute largest eigenvalue of the network??? E' QUALCOSA CHE ABBIAMO O SI RICAVA?
    #Cosa vuol dire in pratica epidemia avviene e epidemia non avviene?
    #Poi fai equazione Beta/Gamma = 1/Lambda
    #L_g  = g.laplacian(normalized = True)
    matrix = g.get_adjacency()
    np_matrix = np.array(matrix.data)
    eigen_values =  np.linalg.eigvals(np_matrix)
    max_eigen = max(eigen_values)

    return round(1/max_eigen,4)
