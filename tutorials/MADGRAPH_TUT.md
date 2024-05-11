# Fast guide to `madgraph` and reading events

## [Madgraph homepage](http://madgraph.phys.ucl.ac.be/index.html)

### [Running madgraph](https://www.niu.edu/spmartin/madgraph/)
### [Tutorial](https://www.niu.edu/spmartin/madgraph/madtutor.html)
### [Syntax](https://www.niu.edu/spmartin/madgraph/madsyntax.html)

## LHE event file

 * In every directory `RESULTS` produced by a simulation of `madgraph`, there is also a parser, called `RESULTS/bin/internal/lhe_parser.py`, which is a python library tailored for `.lhe` files.  
  Each event is a list of particles, identified by the first number in each row (see `number_scheme_montecarlorpp.pdf`), each of which has defined properties, like four-momentum and mass. For example, 6 stands for a top quark $t$ and -13 is an anti-muon $\mu^+$. The full legend (taken from the embedded parser) for a single line in an event is:

  <head>
    
      #    (?P<pid>-?\d+)\s+           #PID
       #    (?P<status>-?\d+)\s+            #status (1 for output particle)
        #    (?P<mother1>-?\d+)\s+       #mother
        #    (?P<mother2>-?\d+)\s+       #mother
        #    (?P<color1>[+-e.\d]*)\s+    #color1
        #    (?P<color2>[+-e.\d]*)\s+    #color2
        #    (?P<px>[+-e.\d]*)\s+        #px
        #    (?P<py>[+-e.\d]*)\s+        #py
        #    (?P<pz>[+-e.\d]*)\s+        #pz
        #    (?P<E>[+-e.\d]*)\s+         #E
        #    (?P<mass>[+-e.\d]*)\s+      #mass
        #    (?P<vtim>[+-e.\d]*)\s+      #displace vertex
        #    (?P<helicity>[+-e.\d]*)\s*      #helicity
        #    ($|(?P<comment>\#[\d|D]*))  #comment/end of string
        #    ''',66) #verbose+ignore case

  </head>

  * A Les Houches Events formatted file will contain many blocks that look somewhat like this:
  
  <head>
         
               <event> 
               #particles12      1=d +2.2752800e+01 2.18028500e+02 7.54677100e-03 1.13631900e-01
                     21 -1    0    0  503  502 +0.0000000000e+00 +0.0000000000e+00 +1.0553638102e+02 1.0553638102e+02 0.0000000000e+00 0.0000e+00 -1.0000e+00  
                     21 -1    0    0  501  503 -0.0000000000e+00 -0.0000000000e+00 -6.5866818640e+02 6.5866818640e+02 0.0000000000e+00 0.0000e+00 1.0000e+00  
                      6  2    1    2  501    0 +1.2066017443e+01 -1.3269389022e+02 -4.9880426120e+02 5.4437313907e+02 1.7257800054e+02 0.0000e+00 0.0000e+00  
                      24 2    3    3    0    0 -2.8691059218e+01 -1.4415203949e+02 -4.4749300375e+02 4.7768458215e+02 7.9558684337e+01 0.0000e+00 0.0000e+00  
                     -6  2    1    2    0  502 -1.2066017443e+01 +1.3269389022e+02 -5.4327544172e+01 2.1983142835e+02 1.6619602148e+02 0.0000e+00 0.0000e+00  
                    -24  2    5    5    0    0 -6.7814483670e+00 +1.0577339091e+02 +2.6094356564e+01 1.3472906603e+02 7.8974727120e+01 0.0000e+00 0.0000e+00  
                    -13  1    4    4    0    0 +9.8450477594e-01 -5.2962748136e-01 -1.0632386129e+02 1.0632973824e+02 0.0000000000e+00 0.0000e+00 1.0000e+00  
                     14  1    4    4    0    0 -2.9675563994e+01 -1.4362241201e+02 -3.4116914246e+02 3.7135484391e+02 0.0000000000e+00 0.0000e+00 -1.0000e+00  
                      5  1    3    3  501    0 +4.0757076661e+01 +1.1458149269e+01 -5.1311257456e+01 6.6688556922e+01 4.7000000000e+00 0.0000e+00 -1.0000e+00  
                     13  1    6    6    0    0 +1.3042126651e+01 -5.4560017190e+00 -8.5906056214e+00 1.6542778705e+01 0.0000000000e+00 0.0000e+00 -1.0000e+00  
                    -14  1    6    6    0    0 -1.9823575018e+01 +1.1122939263e+02 +3.4684962185e+01 1.1818628733e+02 0.0000000000e+00 0.0000e+00 1.0000e+00  
                     -5  1    5    5    0  502 -5.2845690758e+00 +2.6920499316e+01 -8.0421900735e+01 8.5102362315e+01 4.7000000000e+00 0.0000e+00 1.0000e+00  
              <mgrwt>                       
              <rscale>  2 0.21676356E+03</rscale>                
              <asrwt>0</asrwt>                
              <pdfrwt beam="1">  1       21 0.16236366E-01 0.21802849E+03</pdfrwt>                
              <pdfrwt beam="2">  1       21 0.10133357E+00 0.21802849E+03</pdfrwt>                
              <totfact> 0.26611858E+04</totfact>                
              </mgrwt>  

  </head>

  * In the present directory there is a short test that reads an event file and does some operations with the 4-momenta of the particles.
    For execution, paste the file to `RESULTS/bin`, go to `RESULTS` and type in the terminal:  
    `python ./bin/test_with_4momenta.py ./Events/run_01/unweighted_events.lhe`
  

