#!python
import internal.lhe_parser as lhe_parser
print('imported')
import sys
import math

# open and load file of events
eventfile = sys.argv[1]
F=lhe_parser.EventFile(eventfile)

#################################################################

# create vectors for storage
top_momenta = []
antitop_momenta = []
wplus_momenta = []
wminus_momenta = []
b_momenta = []
antib_momenta = []

# parse the file event by event
for iev, event in enumerate(F):
    # find particles in the event file (event is a list)
    for part in event:
        ################    TOP    ################
        # top
        if part.pid == 6:
            top_momenta.append(lhe_parser.FourMomentum(part))
        # antitop
        if part.pid == -6:
            antitop_momenta.append(lhe_parser.FourMomentum(part))
        ################    W    ################
        # W+
        if part.pid == 24:
            wplus_momenta.append(lhe_parser.FourMomentum(part))
        # W-
        if part.pid == -24:
            wminus_momenta.append(lhe_parser.FourMomentum(part))
        ################    b    ################
        # b
        if part.pid == 5:
            b_momenta.append(lhe_parser.FourMomentum(part))
        # antib
        if part.pid == -5:
            antib_momenta.append(lhe_parser.FourMomentum(part))
    
for i in range(len(top_momenta)):
    # find the top and antitop momenta
    print('\nThe momenta of the byproducts should add up to the mother particle momentum:')
    print(top_momenta[i].pz-wplus_momenta[i].pz-b_momenta[i].pz)
    print(antitop_momenta[i].pz-wminus_momenta[i].pz-antib_momenta[i].pz)
    # Lorentz boosts
    # p.boost_to_restframe(q) boosts p to the system where q is zero
    print('\nThe top momentum is zero in the top rest frame:')
    top_comf = top_momenta[i].boost_to_restframe(top_momenta[i]) 
    print(top_comf) # expected to be 0
    # invariant mass
    print('\nIn any frame, p^2=m^2 (invariant mass):')
    print(math.sqrt(top_momenta[i]*top_momenta[i])) 
    # in the top frame, W+ and b are expected to have 
    bottom_comf = b_momenta[i].boost_to_restframe(top_momenta[i])
    wplus_comf = wplus_momenta[i].boost_to_restframe(top_momenta[i])
    print('\nb and W+ momenta in the top rest frame:')
    print(bottom_comf)
    print(wplus_comf)
    break # stops at 1
