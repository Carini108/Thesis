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
def WikiBoostToRF(four_momentum, four_momentum_RF_target):
    """
    https://en.wikipedia.org/wiki/Lorentz_transformation#Proper_transformations
    
        Apply a general Lorentz boost to the given four-momentum vector.

    Parameters:
        four_momentum (np.array): The input four-momentum vector [E, px, py, pz].
        v (np.array): The velocity vector [vx, vy, vz].

    Returns:
        np.array: The boosted four-momentum vector.
    """
    beta_x = four_momentum_RF_target[1]/four_momentum_RF_target[0]
    beta_y = four_momentum_RF_target[2]/four_momentum_RF_target[0]
    beta_z = four_momentum_RF_target[3]/four_momentum_RF_target[0]
    beta = np.sqrt(beta_x**2 + beta_y**2 + beta_z**2)
    gamma = 1 / np.sqrt(1 - beta**2)

    # Construct the Lorentz boost matrix
    Lambda = np.array(
        [[
            gamma, 
          -gamma * beta_x, 
          -gamma * beta_y, 
          -gamma * beta_z
         ],
         [
             -gamma * beta_x, 
             1 + (gamma - 1) * (beta_x**2 / beta**2),
             (gamma - 1) * (beta_x * beta_y / beta**2),
             (gamma - 1) * (beta_x * beta_z / beta**2)
         ],
         [
             -gamma * beta_y, 
             (gamma - 1) * (beta_y * beta_x / beta**2),
             1 + (gamma - 1) * (beta_y**2 / beta**2),
             (gamma - 1) * (beta_y * beta_z / beta**2)
         ],
         [
             -gamma * beta_z, 
             (gamma - 1) * (beta_z * beta_x / beta**2),
             (gamma - 1) * (beta_z * beta_y / beta**2),
             1 + (gamma - 1) * (beta_z**2 / beta**2)
         ]])

    # Apply the Lorentz boost
    boosted_four_momentum = np.dot(Lambda, four_momentum)

    return boosted_four_momentum


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
    # observable D
    N_cos_plus = 0
    N_cos_minus = 0
    COS_PHI = []
    D = []

    #############################################################
    # 2                             parse the file event by event
    #############################################################
    
    print(f"Processing the events...")
    with alive_bar(n_events, force_tty=True) as bar: # to show the overall progress (must need n_events in advance)
        for iev, event in enumerate(F):
                    
            #############################################################
            # 3 find parts. in ev. and save 4-mom ('event' is a list)
            #############################################################

            bool_top = False 
            bool_bottom = False 
            bool_antilep = False
            bool_neutr = False

            bool_antitop = False
            bool_antibottom = False            
            bool_lep = False
            bool_antineutr = False

            for part in event:
                
                # top
                if part.pid == 6: 
                    bool_top = True
                    top_momenta = lhe_parser.FourMomentum(part)
                # bottom
                if part.pid==5:
                    bool_bottom = True
                    bottom_momenta = lhe_parser.FourMomentum(part)
                # anti-lepton
                if part.pid==-11 or part.pid==-13 or part.pid==-15:
                    bool_antilep = True
                    antilep_momenta = lhe_parser.FourMomentum(part)
                # neutrino
                if part.pid==12 or part.pid==14 or part.pid==16:
                    bool_neutr = True
                    neutr_momenta = lhe_parser.FourMomentum(part)

                # antitop
                if part.pid == -6:
                    bool_antitop = True
                    antitop_momenta = lhe_parser.FourMomentum(part)
                # anti-bottom
                if part.pid==-5:
                    bool_antibottom = True
                    antibottom_momenta = lhe_parser.FourMomentum(part)
                # lepton
                if part.pid==11 or part.pid==13 or part.pid==15:
                    bool_lep = True
                    lep_momenta = lhe_parser.FourMomentum(part)
                # anti-neutrino
                if part.pid==-12 or part.pid==-14 or part.pid==-16:
                    bool_antineutr = True
                    antineutr_momenta = lhe_parser.FourMomentum(part)


            # check whether an event lacks one of the particles above
            if bool_top==False or bool_antitop==False or bool_antilep==False or bool_lep==False or bool_antineutr==False or bool_neutr==False or bool_antibottom==False or bool_bottom==False:
                print(f"Particle missing in event number {iev}!")
            if any([abs(v) > 1e-6 for v in (top_momenta -(bottom_momenta + antilep_momenta + neutr_momenta)).to_list()]):
                print("error top")
            if any([abs(v) > 1e-6 for v in (antitop_momenta -(antibottom_momenta + lep_momenta + antineutr_momenta)).to_list()]):
                print("error antitop")

            #############################################################
            # 4               perform boosts and construct helicity basis
            #############################################################

            # 4.1 boost to the c.o.m.f.of tt~

            # c.o.m. momentum
            t_antit_sum = top_momenta+antitop_momenta

            # convert all 4 momenta to np.arrays
            # t tbar pair
            t_antit_4M = np.array([t_antit_sum.E,t_antit_sum.px,t_antit_sum.py,t_antit_sum.pz])
            # top group
            top_4M = np.array([top_momenta.E,top_momenta.px,top_momenta.py,top_momenta.pz])
            bottom_4M = np.array([bottom_momenta.E,bottom_momenta.px,bottom_momenta.py,bottom_momenta.pz])
            antilep_4M = np.array([antilep_momenta.E,antilep_momenta.px,antilep_momenta.py,antilep_momenta.pz])
            neutr_4M = np.array([neutr_momenta.E,neutr_momenta.px,neutr_momenta.py,neutr_momenta.pz])
            # ANTItop group
            antitop_4M = np.array([antitop_momenta.E,antitop_momenta.px,antitop_momenta.py,antitop_momenta.pz])
            antibottom_4M = np.array([antibottom_momenta.E,antibottom_momenta.px,antibottom_momenta.py,antibottom_momenta.pz])
            lep_4M = np.array([lep_momenta.E,lep_momenta.px,lep_momenta.py,lep_momenta.pz])
            antineutr_4M = np.array([antineutr_momenta.E,antineutr_momenta.px,antineutr_momenta.py,antineutr_momenta.pz])

            # perform Lorentz boosts on all momenta to the rest frame of t anti t
            # top group
            top_4M = WikiBoostToRF(top_4M, t_antit_4M)
            bottom_4M = WikiBoostToRF(bottom_4M, t_antit_4M)
            antilep_4M = WikiBoostToRF(antilep_4M, t_antit_4M)
            neutr_4M = WikiBoostToRF(neutr_4M, t_antit_4M)
            # ANTItop group
            antitop_4M = WikiBoostToRF(antitop_4M, t_antit_4M)
            antibottom_4M = WikiBoostToRF(antibottom_4M, t_antit_4M)
            lep_4M = WikiBoostToRF(lep_4M, t_antit_4M)
            antineutr_4M = WikiBoostToRF(antineutr_4M, t_antit_4M)
            
            #print(lep_4M)
            #print(antitop_4M)

            if any([abs(v) > 1e-4 for v in (top_4M -(bottom_4M + antilep_4M + neutr_4M))]):
                print("error 1 top")
            if any([abs(v) > 1e-4 for v in (antitop_4M -(antibottom_4M + lep_4M + antineutr_4M))]):
                print("error 1 antitop")
            

            # 4.2 define the helicity basis

            # beam direction
            p = np.array([0.0,0.0,1.0])
            # top direction
            k = np.array([top_4M[1], top_4M[2], top_4M[3]])
            k = k/np.linalg.norm(k) # normalization
            # r on the (k,p) plane, orthogonal to k
            theta = angle(p,k)
            # only want acute thetas
            if theta>0.5*np.pi: 
                k = -k # invert k instead of p
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

            inv_m_squared = t_antit_sum*t_antit_sum
            beta_tantit = t_antit_sum.pz/t_antit_sum.E
            if inv_m_squared<345**2 or inv_m_squared>400**2 or beta_tantit>0.9: # or theta>(0.5*np.pi)*0.9: # resolve in energy and angle
                continue

            # 5.1 k boost from c.o.m.f. to rest fram of single t / t~

            # top group
            bottom_4M = WikiBoostToRF( bottom_4M , top_4M )
            antilep_4M = WikiBoostToRF( antilep_4M , top_4M )
            neutr_4M = WikiBoostToRF( neutr_4M , top_4M )            
            # ANTItop group
            antibottom_4M = WikiBoostToRF( antibottom_4M , antitop_4M )
            lep_4M = WikiBoostToRF( lep_4M , antitop_4M )
            antineutr_4M = WikiBoostToRF( antineutr_4M , antitop_4M )

            #print(lep_4M)

            # 5.2 get the spatial component to project on the helicity basis
            
            antilep_3momentum = np.array([antilep_4M[1], antilep_4M[2], antilep_4M[3]])
            lep_3momentum = np.array([lep_4M[1], lep_4M[2], lep_4M[3]])

            if any([abs(v) > 1e-6 for v in ((bottom_4M + antilep_4M + neutr_4M))[1:]]):
                print("error 2 top")
            if any([abs(v) > 1e-6 for v in ((antibottom_4M + lep_4M + antineutr_4M))[1:]]):
                print("error 2 antitop")
            
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
            # calculate D
            cos_phi = cosine(antilep_3momentum,lep_3momentum)
            COS_PHI.append(cos_phi)
            D.append(-3.*cos_phi)
            # calculate Ns for asymmetry
            if cos_phi>=0:
                N_cos_plus += 1
            if cos_phi<0:
                N_cos_minus += 1

            #exit(0)

            bar() # progress bar

        # signal end of the for loop
        events_under_threshold = len(cos_k)

        print(f"{events_under_threshold} events found below E<{E_thres} threshold") 
        for j in range(events_under_threshold):
        # for hists of x_ii
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

    print("OBSERVABLE C")
    print(f"C_kk + C_rr = {C_kk+C_rr}") # at threshold, this quantity should be negative
    C = abs(C_kk+C_rr)-C_nn
    err_C = math.sqrt(err_kk**2+err_nn**2+err_rr**2)
    print(f"C = {C} +- {err_C}")

    print("OBSERVABLE D")
    mean_D = np.mean(D)
    err_D = stat.stdev(D)/math.sqrt(len(D)-1)
    print(f"D = {mean_D} +- {err_D}")

    print("ASYMMETRY OBSERVABLE A_D")
    N_tot = N_cos_plus+N_cos_minus
    A_D = (N_cos_plus-N_cos_minus)/(N_tot)
    print(f"A_D = {A_D} with N_tot = {N_tot}")
    print(f"which implies D = {-2.0*A_D}")

    # Writing results to a file
    output_filename = f"kinematic_quantities_Ethr{E_thres}.txt"
    with open(output_filename, "w") as f:
        f.write(f"From the kinematic calculations it follows that:\n")
        f.write(f"OBSERVABLE C\n")
        f.write(f"C_kk = {C_kk} +- {err_kk}\n")
        f.write(f"C_rr = {C_rr} +- {err_rr}\n")
        f.write(f"C_nn = {C_nn} +- {err_nn}\n")
        f.write(f"C_kk + C_rr = {C_kk+C_rr}\n")
        f.write(f"C = {C} +- {err_C}\n")
        f.write("OBSERVABLE D\n")
        f.write(f"D = {mean_D} +- {err_D}\n")
        f.write("ASYMMETRY OBSERVABLE A_D\n")
        f.write(f"A_D = {A_D} with N_tot = {N_tot}\n")
        f.write(f"which implies D = {-2.0*A_D}\n")

    #############################################################
    # 7                            plot histograms of x_ii values
    #############################################################

    plt.figure(figsize=(13, 4.5))
    
    # histogram for x_kk
    plt.subplot(1, 3, 1)
    plt.hist(np.array(x_kk), bins=6, color='b', alpha=0.7, histtype='step', stacked=False, fill=False, density='True')
    x_kk.clear()
    plt.title('x_kk')
    plt.grid()
    plt.xlabel('values')
    plt.ylabel('frequency')

    # histogram for x_rr
    plt.subplot(1, 3, 2)
    plt.hist(np.array(x_rr), bins=6, color='g', alpha=0.7, histtype='step', stacked=False, fill=False, density='True')
    x_rr.clear()
    plt.title('x_rr')
    plt.grid()
    plt.xlabel('values')
    plt.ylabel('frequency')

    # histogram for x_nn
    plt.subplot(1, 3, 3)
    plt.hist(np.array(x_nn), bins=6, color='r', alpha=0.7, histtype='step', stacked=False, fill=False, density='True')
    x_nn.clear()
    plt.title('x_nn')
    plt.grid()
    plt.xlabel('values')
    plt.ylabel('frequency')

    # Adjust layout and save plots
    plt.tight_layout()
    plt.savefig(f'x_histos_Ethr{E_thres}_nev{n_events}.pdf')

    #############################################################
    # 8                           plot 2D histogram using pcolor
    #############################################################

    plt.figure(figsize=(10, 10)) 

    plt.hist2d(THETA, MASS, bins=(10,int((E_thres-300.0)/100)), range=([0.0,0.5],[300.0,E_thres]))
    plt.xlabel('theta / pi$')
    plt.ylabel('m')
    plt.title(f'E threshold ={ E_thres}')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Counts')
    plt.savefig(f'2Dhistograms_Ethr{E_thres}_nev{n_events}.pdf') 

    #############################################################
    # 9                            plot histograms of the cosines
    #############################################################

    plt.figure(figsize=(13, 4.5))
    
    # histogram for cos_k
    plt.subplot(2, 3, 1)
    plt.hist(np.array(cos_k), bins=6, color='b', alpha=0.7, histtype='step', stacked=False, fill=False, density='True')
    cos_k.clear()
    plt.title('cos( theta_k )')
    plt.grid()
    plt.xlabel('values')
    plt.ylabel('frequency')

    # histogram for cos_r
    plt.subplot(2, 3, 2)
    plt.hist(np.array(cos_r), bins=6, color='g', alpha=0.7, histtype='step', stacked=False, fill=False, density='True')
    cos_r.clear()
    plt.title('cos ( theta_r )')
    plt.grid()
    plt.xlabel('values')
    plt.ylabel('frequency')

    # histogram for cos_n
    plt.subplot(2, 3, 3)
    plt.hist(np.array(cos_n), bins=6, color='r', alpha=0.7, histtype='step', stacked=False, fill=False, density='True')
    cos_n.clear()
    plt.title('cos ( theta_n )')
    plt.grid()
    plt.xlabel('values')
    plt.ylabel('frequency')

    # histogram for cos_antik
    plt.subplot(2, 3, 4)
    plt.hist(np.array(cos_antik), bins=6, color='b', alpha=0.7, histtype='step', stacked=False, fill=False, density='True')
    cos_antik.clear()
    plt.title('cos ( antitheta_k )')
    plt.grid()
    plt.xlabel('values')
    plt.ylabel('frequency')

    # histogram for cos_antir
    plt.subplot(2, 3, 5)
    plt.hist(np.array(cos_antir), bins=6, color='g', alpha=0.7, histtype='step', stacked=False, fill=False, density='True')
    cos_antir.clear()
    plt.title('cos ( antitheta_r )')
    plt.grid()
    plt.xlabel('values')
    plt.ylabel('frequency')

    # histogram for cos_antin
    plt.subplot(2, 3, 6)
    plt.hist(np.array(cos_antin), bins=6, color='r', alpha=0.7, histtype='step', stacked=False, fill=False, density='True')
    cos_antin.clear()
    plt.title('cos ( antitheta_n )')
    plt.grid()
    plt.xlabel('values')
    plt.ylabel('frequency')

    # adjust layout and save plots
    plt.tight_layout()
    plt.savefig(f'cos_histos_Ethr{E_thres}_nev{n_events}.pdf')  

    #############################################################
    # 10                            plot histograms of D
    #############################################################

    # plot the histogram for D
    plt.figure()
    plt.hist(D, bins=6, color='blue', alpha=0.7, histtype='step', stacked=False, fill=False, density='True')
    plt.title('Histogram of Observable D')
    plt.xlabel('Values of D')
    plt.ylabel('Frequency')
    plt.grid(True)

    # adjust layout and save plots
    plt.savefig(f'D_histogram_Ethr{E_thres}_nev{n_events}.pdf')  

    #############################################################
    # 11                             plot histogram of COS_PHI
    #############################################################

    # plot the histogram for COS_PHI
    plt.figure()
    plt.hist(COS_PHI, bins=6, color='purple', alpha=0.7, histtype='step', stacked=False, fill=False, density='True')
    plt.title('Histogram of COS_PHI')
    plt.xlabel('Values of COS_PHI')
    plt.ylabel('Frequency')
    plt.grid(True)

    # adjust layout and save plot
    plt.savefig(f'COS_PHI_histogram_Ethr{E_thres}_nev{n_events}.pdf')

if __name__ == "__main__":
    main()
