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

    # 1.2 variable initialization

    # cosines                                    
    cos_k = []
    cos_antik = []
    cos_r = []
    cos_antir = []
    cos_n = []
    cos_antin = []
    # initialize x_ii 												
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
                    top_momenta = lhe_parser.FourMomentum(part)
                # antitop
                if part.pid == -6:
                    bool_antitop = True
                    antitop_momenta = lhe_parser.FourMomentum(part)
                # anti-lepton
                if part.pid==-11 or part.pid==-13 or part.pid==-15:
                    bool_antilep = True
                    antilep_momenta = lhe_parser.FourMomentum(part)
                # lepton
                if part.pid==11 or part.pid==13 or part.pid==15:
                    bool_lep = True
                    lep_momenta = lhe_parser.FourMomentum(part)
            # check whether an event lacks one of the particles above
            if bool_top==False or bool_antitop==False or bool_antilep==False or bool_lep==False:
                print(f"Particle missing in event number {iev}!")

            #############################################################
            # 4               perform boosts and construct helicity basis
            #############################################################

            # 4.1 boost to the c.o.m.f.of tt~

            # c.o.m. momentum
            t_antit_sum = top_momenta+antitop_momenta

            # perform Lorentz boosts on all momenta
            top_momenta = top_momenta.boost_to_restframe(t_antit_sum)
            antitop_momenta = antitop_momenta.boost_to_restframe(t_antit_sum)
            antilep_momenta = antilep_momenta.boost_to_restframe(t_antit_sum)
            lep_momenta = lep_momenta.boost_to_restframe(t_antit_sum)

            # 4.2 define the helicity basis

            # beam direction
            p = np.array([0,0,1])
            # top direction
            k = np.array([top_momenta.px, top_momenta.py, top_momenta.pz])
            k = k/np.linalg.norm(k) # normalization
            # r on the (k,p) plane, orthogonal to k
            theta = angle(p,k)
            # only want acute thetas
            if theta>0.5*np.pi: 
                p = -p
                theta = np.pi-theta
            if theta<0 or theta>np.pi: # check
                print(f'error theta = {theta}')
            r = (p-k*cosine(p,k))/np.sin(theta)
            # orthogonal to k and r
            n = np.cross(k,r)

            # check basis is orthonormal
            if abs(np.inner(k,r))>10**(-6):
                print('kr not orthogonal', np.inner(k,r)) 
            if abs(np.inner(k,n))>10**(-6):
                print('kn not orthogonal', np.inner(k,n)) 
            if abs(np.inner(n,r))>10**(-6):
                print('nr not orthogonal', np.inner(n,r))
            if abs(np.inner(n,p))>10**(-6):
                print('np not orthogonal', np.inner(n,p))

            #############################################################
            # 5                         prepare to calculate the C matrix
            #############################################################

            if t_antit_sum*t_antit_sum>E_thres**2: # or theta>(0.5*np.pi)*0.9: # resolve in energy and angle
                continue

            # 5.1 k boost from c.o.m.f. to rest fram of single t / t~

            antilep_momenta = antilep_momenta.boost_to_restframe(top_momenta)
            lep_momenta = lep_momenta.boost_to_restframe(antitop_momenta)      

            # 5.2 get the spatial component to project on the helicity basis
            
            antilep_3momentum = np.array([antilep_momenta.px, antilep_momenta.py, antilep_momenta.pz])
            lep_3momentum = np.array([lep_momenta.px, lep_momenta.py, lep_momenta.pz])
            
            # 5.3 calculate the kinematic quantities

            # for 2d hist
            MASS.append(np.sqrt(t_antit_sum*t_antit_sum))
            THETA.append(theta/np.pi)
            # for hists of cosines
            cos_k.append(cosine(antilep_3momentum, k))
            cos_antik.append(cosine(lep_3momentum, k))
            cos_r.append(cosine(antilep_3momentum, r))
            cos_antir.append(cosine(lep_3momentum, r))
            cos_n.append(cosine(antilep_3momentum, n))
            cos_antin.append(cosine(lep_3momentum, n))
            
            bar() # progress bar

        # signal end of the for loop
        events_under_threshold = len(cos_k)
        print(f"{events_under_threshold} events found below E<{E_thres} threshold") 
        for j in range(events_under_threshold):
        # for hists of x_iis
            x_kk.append(cos_k[j]*cos_antik[j])
            x_rr.append(cos_r[j]*cos_antir[j])
            x_nn.append(cos_n[j]*cos_antin[j])

    #############################################################
    # 6          print the results with statistical uncertainties
    #############################################################

    # determine the matrix elements for entanglement detection
    print(f"From the kinematic calculations it follows that:")

    C_kk = -9.*np.mean(x_kk)
    err_kk = 9.*stat.stdev(x_kk)/math.sqrt(len(x_kk)-1) # error propagation ==> times |-9|
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
    x_kk.clear()
    plt.title('x_kk')
    plt.grid()
    plt.xlabel('values')
    plt.ylabel('frequency')

    # histogram for x_rr
    plt.subplot(1, 3, 2)
    plt.hist(np.array(x_rr), bins=6, color='g', alpha=0.7, edgecolor='black', density='True')
    x_rr.clear()
    plt.title('x_rr')
    plt.grid()
    plt.xlabel('values')
    plt.ylabel('frequency')

    # histogram for x_nn
    plt.subplot(1, 3, 3)
    plt.hist(np.array(x_nn), bins=6, color='r', alpha=0.7, edgecolor='black', density='True')
    x_nn.clear()
    plt.title('x_rr')
    plt.grid()
    plt.xlabel('values')
    plt.ylabel('frequency')

    # Adjust layout and save plots
    plt.tight_layout()
    plt.savefig(f'x_histos_Ethr{E_thres}.pdf')

    #############################################################
    # 8                           plot 2D histogram using pcolorù
    #############################################################

    fig2 = plt.figure()
    plt.hist2d(THETA, MASS, bins=(10,int((E_thres-300.0)/100)), range=([0.0,0.5],[300.0,E_thres]))
    plt.xlabel('theta/pi')
    plt.ylabel('m')
    plt.title(f'$E_thr$={E_thres}')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Counts')
    plt.savefig(f'2Dhistograms_Ethr{E_thres}.pdf') 

    #############################################################
    # 9                            plot histograms of the cosines
    #############################################################

    plt.figure(figsize=(13, 4.5))
    
    # histogram for cos_k
    plt.subplot(2, 3, 1)
    plt.hist(np.array(cos_k), bins=6, color='b', alpha=0.7, edgecolor='black', density='True')
    cos_k.clear()
    plt.title('cos_k')
    plt.grid()
    plt.xlabel('values')
    plt.ylabel('frequency')

    # histogram for cos_r
    plt.subplot(2, 3, 2)
    plt.hist(np.array(cos_r), bins=6, color='g', alpha=0.7, edgecolor='black', density='True')
    cos_r.clear()
    plt.title('cos_r')
    plt.grid()
    plt.xlabel('values')
    plt.ylabel('frequency')

    # histogram for cos_n
    plt.subplot(2, 3, 3)
    plt.hist(np.array(cos_n), bins=6, color='r', alpha=0.7, edgecolor='black', density='True')
    cos_n.clear()
    plt.grid()
    plt.xlabel('values')
    plt.ylabel('frequency')

    # histogram for cos_antik
    plt.subplot(2, 3, 4)
    plt.hist(np.array(cos_antik), bins=6, color='b', alpha=0.7, edgecolor='black', density='True')
    cos_antik.clear()
    plt.title('cos_antik')
    plt.grid()
    plt.xlabel('values')
    plt.ylabel('frequency')

    # histogram for cos_antir
    plt.subplot(2, 3, 5)
    plt.hist(np.array(cos_antir), bins=6, color='g', alpha=0.7, edgecolor='black', density='True')
    cos_antir.clear()
    plt.title('cos_antir')
    plt.grid()
    plt.xlabel('values')
    plt.ylabel('frequency')

    # histogram for cos_antin
    plt.subplot(2, 3, 6)
    plt.hist(np.array(cos_antin), bins=6, color='r', alpha=0.7, edgecolor='black', density='True')
    cos_antin.clear()
    plt.grid()
    plt.xlabel('values')
    plt.ylabel('frequency')

    # adjust layout and save plots
    plt.tight_layout()
    plt.savefig(f'cos_histos_Ethr{E_thres}.pdf')    

if __name__ == "__main__":
    main()
