#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

# LIBRARIES
# managing files interactively
import sys
from alive_progress import alive_bar
# parser for Les Houches Events file
import internal.lhe_parser as lhe_parser
# statistics and vectors
import math
import statistics as stat
import numpy as np
import matplotlib.pyplot as plt
print('Libraries properly imported')

# FUNCTIONS
# functions that return the cosine and the angle between vectors (in radians)
def cosine(v, w): return (np.inner(v,w)/(np.linalg.norm(v)*np.linalg.norm(w)))
def angle(v, w): return np.arccos(cosine(v,w))

# MAIN
def main():

    #############################################################
    # 1        variables initialization and loading of event file 
    #############################################################

    # 1.1 open file

    # from command line 
    eventfile = sys.argv[1]
    n_events = int(sys.argv[2])
    # skim events with too high energy
    E_thres = 400.0 # energy threshold in GeV
    # read event file
    F = lhe_parser.EventFile(eventfile)

    # 1.2 lists initialization

    # initialize 4-momenta lists for storage						4-vectors
    top_momenta = []
    antitop_momenta = []
    antilep_momenta = []
    lep_momenta = []
    # vectors for boosting momenta and constructing helicity basis
    t_antit_sum = [] # sum of the momenta of top and antitop quark
    p = [] # beam direction											3-vectors	
    # helicity (orthonormal) basis in the top rest frame 
    k = [] # top direction
    r = [] # on (p,k) plane and orthogonal to k
    n = [] # orthogonal to (p,k) plane
    # initialize x_ii 												scalars
    x_kk = []
    x_rr = []
    x_nn = []

    # 2D histograms
    THETA = []
    MASS = []

    #############################################################
    # 2                             parse the file event by event
    #############################################################
    
    print(f"Processing the events...")
    with alive_bar(n_events, force_tty=True) as bar: # to show the overall progress (must need n_events in advance)
        for iev, event in enumerate(F):
            """cuts = True
            for part in event:
                if abs(part.pid)==5:
                    if lhe_parser.FourMomentum(part).pt2() < 25.**2:
                        cuts=False
                    if abs(lhe_parser.FourMomentum(part).pseudorapidity()) > 2.5:
                        cuts=False"""
                    
            #############################################################
            # 3 find parts. in ev. and save 4-mom ('event' is a list)
            #############################################################

            bool_top = False 
            bool_antitop = False
            bool_antilep = False
            bool_lep = False
            for part in event:
                # top
                if part.pid == 6: 
                    bool_top = True
                    top_momenta.append(lhe_parser.FourMomentum(part))
                # antitop
                if part.pid == -6:
                    bool_antitop = True
                    antitop_momenta.append(lhe_parser.FourMomentum(part))
                # anti-lepton
                if part.pid==-11 or part.pid==-13 or part.pid==-15:
                    bool_antilep = True
                    antilep_momenta.append(lhe_parser.FourMomentum(part))
                # lepton
                if part.pid==11 or part.pid==13 or part.pid==15:
                    bool_lep = True
                    lep_momenta.append(lhe_parser.FourMomentum(part))
            # check whether an event lacks one of the particles above
            if bool_top==False or bool_antitop==False or bool_antilep==False or bool_lep==False:
                print(f"Particle missing in event number {iev}!")

            #############################################################
            # 4               perform boosts and construct helicity basis
            #############################################################

            # 4.1 boost to the c.o.m.f.of tt~

            # c.o.m. momentum
            t_antit_sum.append(top_momenta[iev]+antitop_momenta[iev])
            # perform Lorentz boosts on all momenta
            top_momenta[iev] = top_momenta[iev].boost_to_restframe(t_antit_sum[iev])
            antitop_momenta[iev] = antitop_momenta[iev].boost_to_restframe(t_antit_sum[iev])
            antilep_momenta[iev] = antilep_momenta[iev].boost_to_restframe(t_antit_sum[iev])
            lep_momenta[iev] = lep_momenta[iev].boost_to_restframe(t_antit_sum[iev])

            # 4.2 define the helicity basis

            # beam direction
            p.append(np.array([0,0,1]))
            # top direction
            k.append(np.array([top_momenta[iev].px, top_momenta[iev].py, top_momenta[iev].pz]))
            k[iev] = k[iev]/np.linalg.norm(k[iev]) # normalization
            # r on the (k,p) plane, orthogonal to k
            theta = angle(p[iev],k[iev])
            if theta>0.5*np.pi:
                p[iev] = -p[iev]
                theta = np.pi-theta
            if theta<0 or theta>np.pi: # check
                print(f'error theta = {theta}')
            r.append((p[iev]-k[iev]*cosine(p[iev],k[iev]))/np.sin(theta))
            # orthogonal to k and r
            n.append(np.cross(k[iev],r[iev]))

            # check basis is orthonormal
            if abs(np.inner(k[iev],r[iev]))>10**(-6):
                print('kr', np.inner(k[iev],r[iev])) 
            if abs(np.inner(k[iev],n[iev]))>10**(-6):
                print('kn', np.inner(k[iev],n[iev])) 
            if abs(np.inner(n[iev],r[iev]))>10**(-6):
                print('nr', np.inner(n[iev],r[iev]))
            if abs(np.inner(n[iev],p[iev]))>10**(-6):
                print('np', np.inner(n[iev],p[iev]))

            #############################################################
            # 5                         prepare to calculate the C matrix
            #############################################################

            # 5.1 k boost from c.o.m.f. to rest fram of single t / t~

            antilep_momenta[iev] = antilep_momenta[iev].boost_to_restframe(top_momenta[iev])
            lep_momenta[iev] = lep_momenta[iev].boost_to_restframe(antitop_momenta[iev])      

            # 5.2 get the spatial component to project on the helicity basis
            
            antilep_3momentum = np.array([antilep_momenta[iev].px, antilep_momenta[iev].py, antilep_momenta[iev].pz])
            lep_3momentum = np.array([lep_momenta[iev].px, lep_momenta[iev].py, lep_momenta[iev].pz])
            
            # 5.3 calculate the x_ii
            
            if t_antit_sum[iev]*t_antit_sum[iev]>E_thres**2: # or theta>(0.5*np.pi)*0.9:
                continue
            MASS.append(np.sqrt(t_antit_sum[iev]*t_antit_sum[iev]))
            THETA.append(theta/np.pi)
            x_kk.append(cosine(antilep_3momentum, k[iev])*cosine(lep_3momentum, k[iev]))
            x_rr.append(cosine(antilep_3momentum, r[iev])*cosine(lep_3momentum, r[iev]))
            x_nn.append(cosine(antilep_3momentum, n[iev])*cosine(lep_3momentum, n[iev]))
            
            bar() # progress bar

    #############################################################
    # 6          print the results with statistical uncertainties
    #############################################################

    print(f"{len(x_kk)} events found below E<{E_thres} threshold") # signal end of the for loop

    # determine the matrix elements for entanglement detection
    print(f"From the kinematic calculations it follows that:")

    C_kk = -9.*np.mean(x_kk)
    err_kk = 9.*stat.stdev(x_kk)/math.sqrt(len(x_kk)-1)
    print(f"C_kk = {C_kk} +- {err_kk}")

    C_rr = -9.*np.mean(x_rr)
    err_rr = 9.*stat.stdev(x_rr)/math.sqrt(len(x_rr)-1)
    print(f"C_rr = {C_rr} +- {err_rr}")

    C_nn = -9.*np.mean(x_nn)
    err_nn = 9.*stat.stdev(x_nn)/math.sqrt(len(x_nn)-1)
    print(f"C_nn = {C_nn} +- {err_nn}")

    # find entanglement
    print(f"C_kk + C_rr = {C_kk+C_rr}") # at threshold, this quantity should be negative
    C = abs(C_kk+C_rr)-C_nn
    err_C = math.sqrt(err_kk**2+err_nn**2+err_rr**2)
    print(f"C = {C} +- {err_C}")

    #############################################################
    # 7                            plot histograms of x_ii values
    #############################################################

    plt.figure(figsize=(13, 4.5))
    
    # histogram for x_kk
    plt.subplot(1, 3, 1)
    plt.hist(np.array(x_kk), bins=6, color='b', alpha=0.7, edgecolor='black', density='True')
    plt.title('x_kk')
    plt.grid()
    plt.xlabel('values')
    plt.ylabel('frequency')

    # histogram for x_rr
    plt.subplot(1, 3, 2)
    plt.hist(np.array(x_rr), bins=6, color='g', alpha=0.7, edgecolor='black', density='True')
    plt.title('x_rr')
    plt.grid()
    plt.xlabel('values')
    plt.ylabel('frequency')

    # histogram for x_nn
    plt.subplot(1, 3, 3)
    plt.hist(np.array(x_nn), bins=6, color='r', alpha=0.7, edgecolor='black', density='True')
    plt.title('x_nn')
    plt.grid()
    plt.xlabel('values')
    plt.ylabel('frequency')

    # Adjust layout and show plot
    plt.tight_layout()
    plt.savefig('histograms.png')
    plt.savefig('histograms.pdf')
    plt.show()

    ############################
    # plot 2D histogram using pcolorù

    fig2 = plt.figure()
    plt.hist2d(THETA, MASS, bins=(10,int((E_thres-300.0)/100)), range=([0.0,0.5],[300.0,E_thres]))
    plt.xlabel('theta/pi')
    plt.ylabel('m')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Counts')
    plt.savefig('2Dhistograms.pdf') 

if __name__ == "__main__":
    main()


    # istogrammi coseni piatti