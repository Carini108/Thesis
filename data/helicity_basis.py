#!python
import internal.lhe_parser as lhe_parser
print('imported')
import sys
import math
import numpy as np

# Note: returns angle in radians
def angle(v, w): return np.arccos(v.dot(w)/(np.linalg.norm(v)*np.linalg.norm(w)))

# open and load file of events
eventfile = sys.argv[1]
F=lhe_parser.EventFile(eventfile)

#################################################################

# create vectors for storage

top_momenta = []
antilep_momenta = []
antitop_momenta = []
lep_momenta = []

# boosting and helicity basis

t_antit_sum = []
p = []
k = []
r = []
n = []
antik = []
antir = []
antin = []

# entanglement matrices




# parse the file event by event
for iev, event in enumerate(F):

    # find particles in the event file (event is a list)
    for part in event:

        ################    top & anti-top    ################
        # top
        if part.pid == 6:
            top_momenta.append(lhe_parser.FourMomentum(part))
        # antitop
        if part.pid == -6:
            antitop_momenta.append(lhe_parser.FourMomentum(part))
        ################    anti-lepton & lepton    ################
        # anti-lepton
        bool_antilep = part.pid==-11 or part.pid==-13 or part.pid==-15
        if bool_antilep:
            antilep_momenta.append(lhe_parser.FourMomentum(part))
        # lepton
        bool_lep = part.pid==11 or part.pid==13 or part.pid==15
        if bool_lep:
            lep_momenta.append(lhe_parser.FourMomentum(part))

    #############################################################
    # BOOST IN THE C.O.M.F. OF tt~
    # comf
    t_antit_sum.append(top_momenta[iev])
    t_antit_sum[iev]+=antitop_momenta[iev]
    # boosts
    top_momenta[iev] = top_momenta[iev].boost_to_restframe(t_antit_sum[iev])
    antitop_momenta[iev] = antitop_momenta[iev].boost_to_restframe(t_antit_sum[iev])
    antilep_momenta[iev] = antilep_momenta[iev].boost_to_restframe(t_antit_sum[iev])
    lep_momenta[iev] = lep_momenta[iev].boost_to_restframe(t_antit_sum[iev])
    # DEFINE HELICITY BASIS
    ## beam direction
    p.append(np.array([t_antit_sum[iev].px, t_antit_sum[iev].py, t_antit_sum[iev].pz]))
    p[iev] = p[iev]/math.sqrt(np.inner(p[iev], p[iev])) # normalization
    ## top direction
    k.append(np.array([top_momenta[iev].px, top_momenta[iev].py, top_momenta[iev].pz]))
    k[iev] = k[iev]/np.linalg.norm(k[iev]) # normalization
    ## on the k,p plane
    theta = angle(p[iev],k[iev])
    r.append((p[iev]-k[iev]*np.cos(theta))/np.sin(theta))
    ## orthogonal to k and r
    n.append(np.cross(k[iev],r[iev]))

    #############################################################
    # BOOST IN THE C.O.M.F. OF SINGLE t / SINGLE t~
    antilep_momenta[iev] = antilep_momenta[iev].boost_to_restframe(top_momenta[iev])
    lep_momenta[iev] = lep_momenta[iev].boost_to_restframe(antitop_momenta[iev])
    
for i in range(len(top_momenta)):
    print(top_momenta[i])
    if i==0: break # stops at 1st element
