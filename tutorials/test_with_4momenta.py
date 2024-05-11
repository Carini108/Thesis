#!python
import internal.lhe_parser as lhe_parser
print('imported')
import sys
import math

# open and load file of events (input file name from command line!)
eventfile = sys.argv[1]
F=lhe_parser.EventFile(eventfile)

# initialize lists for storage
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
        # top
        if part.pid == 6:
            top_momenta.append(lhe_parser.FourMomentum(part))
        # antitop
        if part.pid == -6:
            antitop_momenta.append(lhe_parser.FourMomentum(part))
        # W+
        if part.pid == 24:
            wplus_momenta.append(lhe_parser.FourMomentum(part))
        # W-
        if part.pid == -24:
            wminus_momenta.append(lhe_parser.FourMomentum(part))
        # b
        if part.pid == 5:
            b_momenta.append(lhe_parser.FourMomentum(part))
        # antib
        if part.pid == -5:
            antib_momenta.append(lhe_parser.FourMomentum(part))
            
# perform some operations on the 4-momenta lists
for i in range(len(top_momenta)):
    
    # find the top and antitop momenta and check if they are zero (within a threshold due to trucation in number representation)
    print('\nThe momenta of the byproducts should add up to the mother particle momentum:')
    print(top_momenta[i].pz-wplus_momenta[i].pz-b_momenta[i].pz)
    print(antitop_momenta[i].pz-wminus_momenta[i].pz-antib_momenta[i].pz)
    
    # tests with Lorentz boosts
    # p.boost_to_restframe(q) boosts p to the system where q is zero
    print('\nThe top momentum is zero in the top rest frame:')
    top_comf = top_momenta[i].boost_to_restframe(top_momenta[i]) 
    print(top_comf) # expected to be 0
    
    # calculate the invariant mass
    print('\nIn any frame, p^2=m^2 (invariant mass):')
    print(math.sqrt(top_momenta[i]*top_momenta[i])) 
    
    # in the top frame, W+ and b are expected to have equal and opposite spatial parts (they come from a 2-body decay)
    bottom_comf = b_momenta[i].boost_to_restframe(top_momenta[i])
    wplus_comf = wplus_momenta[i].boost_to_restframe(top_momenta[i])
    print('\nb and W+ momenta in the top rest frame:')
    print(bottom_comf)
    print(wplus_comf)
    
    break # stops at 1
