__includes ["Network.nls" "Quarantine.nls" "Avoidance.nls" "RiskCulture.nls" "Calibrage.nls" ]

globals 
[ 
  t 
  Reference-Node
  max_I/S_node
  min_I/S_node
  duration
  MaxI
  TimeofMaxI
  TimeToReachAllNodes
  possible-colors 
  QCostTotal 
  ACostTotal
  RCCostTotal
  TotalCost_TrafficLoss_Emissions
  transparency
  mobilepopulationtheo
  mobilepopratio
  maxmobilepopratio
  TimeofMaxMobileRatio
  TotalPop
  ; pour la calibration
  tempS ;tmpS
tempI
tempR
calibS  ; currentS
calibI
calibR
calibKS  ;currentKS
calibKI
calibKR
TotalS
TotalI
TotalR
nxtKS
nxtKI
nxtKR
nxtS
nxtI
nxtR
Total-Population 
Total-Mobile-Population
;Calibration-Step
Calibration-Duration
Epidemy-Duration
Mobile-Population-By-Step
 ]

; Initialization of the number of S, I and R in nodes
; This value is randomly fixed between 0 and (S,I,R)Init

to initialize-your-state
  initialize-population
;  initialize-mobility-rate

;  set gsInit gs
;  set giInit gi
;  set grInit gr
end

to initialize-population
  ifelse One-Infected? [one-infected-guy]  ;If One-Infected? switch is on, we ask only one node of the network to have 1 infected (all the other nodes have 0 infected)
  [
   ifelse Initial-Population = "Uniform"        
   [
     set currentS SInit 
     set currentI IInit 
     set currentR Rinit  
     set nextS currentS  
     set nextI currentI  
     set nextR currentR 
   ]
   [                                      
     set currentS max (list 1 (random-normal SInit (SInit * pop-SD))) ; population standard deviation
     set currentI max (list 1 (random-normal IInit (IInit * pop-SD)))
     set currentR max (list 1 (random-normal RInit (RInit * pop-SD)))
     set nextS currentS  
     set nextI currentI  
     set nextR currentR 
   ]
  ]
end

to one-infected-guy  ; to have no infected persons on the nodes
  ifelse Initial-Population = "Uniform"        
   [
     set currentS SInit 
     set currentI 0
     set currentR Rinit  
     set nextS currentS  
     set nextI currentI  
     set nextR currentR 
   ]
   [                                      
     set currentS max (list 1 (random-normal SInit (SInit * pop-SD))) ; population standard deviation
     set currentI 0
     set currentR max (list 1 (random-normal RInit (RInit * pop-SD)))
     set nextS currentS  
     set nextI currentI  
     set nextR currentR 
   ]
end

to infect-one-guy  ; to infect one person on all the network
  if One-Infected? 
  [  
    ask one-of nodes 
    [
      set currentI 1 
      set nextI currentI
      set shape "star"
    ]
  ]
end

to update-node-infection-state
  if any? nodes with [currenti > 0]
  [
    set TimeToReachAllNodes t 
  ask nodes with [currenti > 0]
  [
   set infected? 1 
  ]
  ]
end

to compute-timeToReachAllNodes
  ask nodes with [currenti > 0 and infected? = 0]
  [
    set TimeToReachAllNodes t 
    set infected? 1
  ]
end


;to initialize-mobility-rate
;   ifelse Mobility-Rate = "Uniform"        
;   [
;     set gs gsInit
;     set gi giInit
;     set gr grInit
;   ]
;   [                                      
;     set gs max (list 0.0001 (random-normal gsInit (gsInit * mob-SD))) ; mobility stadard deviation
;     set gi max (list 0.0001 (random-normal giInit (giInit * mob-SD)))
;     set gr max (list 0.0001 (random-normal grInit (grInit * mob-SD)))
;   ]
;end

to setup
  clear-all
  fix-random-seed
  define-colors
  setup-network
  reset-ticks
    do-calibrate
  ask nodes [
  set initial-gs gs
  set initial-gi gi
  set initial-gr gr
  ]
  reset-ticks
  infect-one-guy ; if One-Infected? is on then only one person is infected during the setup
  update-node-infection-state
  set TotalPop sum [currentS + currentI + currentR] of nodes
end

;make sure colors of nodes and edges are visible
to define-colors
  set possible-colors [1 3 5 7 9 11 13 15 17 19 21 23 25 27]
  
end


; random-seed fixed to be sure the network topology remains the same at each "setup"
to fix-random-seed
;  if FixRandomSeed? [random-seed 47822] 
if FixRandomSeed? [random-seed 47512]
end
; verifier que ça n'a pas d'ifluence sur l'ordre de tirage des noeuds dans tous les appels ask-nodes

to go
  tick
  if S1_Quarantine? [update-quarantine-state compute-QCostOutNodes compute-QCostTotal]
  if S2_Avoidance? [update-avoidance-state compute-ACostOutNodes compute-ACostTotal ]
  if S3_RiskCulture? [update-riskculture-state compute-RCCost compute-RCCostTotal]
  set t t + deltaT
  RKstepNetwork
  Compute-TotalCost
  display-nodes
  FindMaxI
  compute-mobilepopulation
  compute-timeToReachAllNodes
  Compute-Duration
  let epsilon 1.0
  if (ticks - TimeofMaxI) > 0 and ( sum [currentI] of nodes ) < epsilon
  [
;    set duration ticks
    stop
  ]
   
end

to Compute-TotalCost
  set TotalCost_TrafficLoss_Emissions QCostTotal  + ACostTotal + RCCostTotal
end


; Node procedure needed to compute RK4 : first step of the RK4 global method
; SIR model

to calcK [coef]
  set tmpS currentS + coef * currentKS
  set tmpI currentI + coef * currentKI
  set tmpR currentR + coef * currentKR
  set nextKS (- beta / (tmpS + tmpI + tmpR)) * tmpI * tmpS
  set nextKI (beta / (tmpS + tmpI + tmpR)) * tmpI * tmpS - alpha * tmpI
  set nextKR alpha * tmpI
end

; Node procedure needed to compute RK4 : second step of the RK4 global method

to coupling 
 
 ; Outflow : sum of mij for each out-link of the nodes  
 ; this outflow is multiplied by the population of the nodes later
 
  let sumOutS 0
  let sumOutI 0
  let sumOutR 0
  
    ask my-out-links [
    set sumOutS sumOutS + mij
    set sumOutI sumOutI + mij
    set sumOutR sumOutR + mij
    ]
    
  ; Inflow : sum of mij * gs * tmpS, mij * gi * tmpI, mij * gr * tmpR for each in-link
  ; where tmpS, tmpI and tmpR are the number of S, I, R of the source node during the RK4 procedure 
    
  let sumInS 0
  let sumInI 0
  let sumInR 0
  
  ask in-link-neighbors [
    let mi [Mij] of out-link-to myself
    set sumInS sumInS + gs * tmpS * mi
    set sumInI sumInI + gi * tmpI * mi
    set sumInR sumInR + gr * tmpR * mi
  ]
  
  ; Update of the number of S, I and R in each node, taking into account outflows and inflows
   
  set nextKS nextKS - sumOutS * gs * tmpS + sumInS 
  set nextKI nextKI - sumOutI * gi * tmpI + sumInI
  set nextKR nextKR - sumOutR * gr * tmpR + sumInR  
end

; Last step of RK4 method

to stepK [coef] 
  set nextS nextS + coef * nextKS
  set nextI nextI + coef * nextKI
  set nextR nextR + coef * nextKR
  set currentKS nextKS 
  set currentKI nextKI
  set currentKR nextKR
end

; All steps of RK4 methods to compute the new values of S,I and R in nodes

to RKstepNetwork
  (foreach (list 0 (1 / 2) (1 / 2) 1) (list (1 / 6) (1 / 3) (1 / 3) (1 / 6)) [
    ask nodes [calcK ?1 * deltaT]
    ask nodes [coupling]
    ask nodes [stepK ?2 * deltaT]
  ])
  ask nodes [
    set currentS nextS
    set currentI nextI
    set currentR nextR
  ]
end

to FindMaxI
  if MaxI <= sum [currentI] of nodes 
  [
    set MaxI  sum [currentI] of nodes     
    set TimeofMaxI ticks
  ]

end


;to-report maxITheoMin
;  let const alpha / beta
;  report  ( NB_Nodes * ( ( Iinit - ( 0.01 * Iinit ) ) + ( Sinit - ( 0.01 * Iinit ) ) - const + const * ln( const ) - const * ln( ( Sinit - ( 0.01 * Iinit ) ) ) ) ) / ( ( ( Sinit + Iinit + Rinit ) - ( 0.01 * ( Sinit + Iinit + Rinit ) ) ) * NB_Nodes )
;end

to-report maxITheo
  let N1 (SInit + IInit + RInit )
  let const ( alpha * N1 ) / beta
  report NB_Nodes * ( Iinit  + Sinit  - const + const * ln( const ) - const * ln( Sinit ) )
end


to compute-mobilepopulation
  
  ; population mobile totale theorique sans stratégies
  set mobilepopulationtheo mobilepopulationtheo + sum [initial-gs * currentS + initial-gi * currentI + initial-gr * currentR ] of nodes
 
  ; ratio entre le cout total et le population mobile theorique
  set mobilepopratio TotalCost_TrafficLoss_Emissions / mobilepopulationtheo
  
  ; calcul du taux de perte de trafic max
  if mobilepopratio != 0
  [
    if maxmobilepopratio <= mobilepopratio
    [
      set maxmobilepopratio mobilepopratio
      set TimeofMaxMobileRatio t
    ]
  ]
end

to Compute-Duration
      set duration ticks
end
@#$#@#$#@
GRAPHICS-WINDOW
522
10
915
424
10
10
18.24
1
10
1
1
1
0
0
0
1
-10
10
-10
10
1
1
1
ticks
30.0

INPUTBOX
17
182
108
242
Alpha
0.2
1
0
Number

INPUTBOX
113
182
206
242
Beta
0.5
1
0
Number

BUTTON
14
18
80
51
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
104
18
167
51
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
15
285
120
345
SInit
40000
1
0
Number

INPUTBOX
128
285
233
345
IInit
0
1
0
Number

INPUTBOX
238
284
345
344
RInit
0
1
0
Number

TEXTBOX
17
161
167
179
SIR Parameters Values
12
0.0
1

TEXTBOX
20
259
222
289
Initial Conditions of nodes
12
0.0
1

TEXTBOX
13
440
215
470
Mobility Rate of nodes by step
12
0.0
1

MONITOR
565
500
648
545
S
[CurrentS] of reference-node
2
1
11

MONITOR
673
500
754
545
I
[currentI] of reference-node
2
1
11

MONITOR
774
501
854
546
R
[currentR] of reference-node
2
1
11

INPUTBOX
16
78
90
138
NB_Nodes
10
1
0
Number

MONITOR
620
445
784
490
Population Reference Node
[currentS + currentI + currentR]  of reference-node
2
1
11

MONITOR
642
585
744
630
Population
sum [currentS + currentI + currentR] of nodes
2
1
11

CHOOSER
253
78
391
123
Link-Weights
Link-Weights
"Uniform" "Random"
0

INPUTBOX
255
11
353
71
deltaT
0.01
1
0
Number

TEXTBOX
11
57
161
75
Network Parameters
12
0.0
1

TEXTBOX
663
429
813
447
Upper Node
12
0.0
1

TEXTBOX
638
567
788
585
Total Population
12
0.0
1

MONITOR
555
642
637
687
S glob
sum [currentS] of nodes
2
1
11

MONITOR
662
644
746
689
I glob
sum [currentI] of nodes
2
1
11

MONITOR
763
644
844
689
R glob
sum [currentR] of nodes
2
1
11

SWITCH
19
567
173
600
S1_Quarantine?
S1_Quarantine?
1
1
-1000

INPUTBOX
181
553
300
613
Quarantine_Threshold
0.04
1
0
Number

TEXTBOX
21
543
171
561
Strategies
12
0.0
1

SWITCH
16
626
171
659
S2_Avoidance?
S2_Avoidance?
0
1
-1000

INPUTBOX
181
617
300
677
Avoidance_Threshold
0.035
1
0
Number

SWITCH
19
693
170
726
S3_RiskCulture?
S3_RiskCulture?
1
1
-1000

INPUTBOX
182
682
347
742
RiskCulture_Probability
0.955
1
0
Number

MONITOR
1288
332
1481
377
Duration of Epidemic (ticks)
duration
2
1
11

MONITOR
925
332
1108
377
Maximum number of infected
MaxI
2
1
11

MONITOR
1115
332
1284
377
Time Maximum Infected
TimeofMaxI
2
1
11

MONITOR
514
804
762
849
Total Risk Culture Cost
RCCostTotal
2
1
11

CHOOSER
370
246
514
291
Initial-Population
Initial-Population
"Uniform" "Random"
0

TEXTBOX
528
747
678
801
*************\nRisk Culture \n*************
15
0.0
1

TEXTBOX
264
749
414
803
***********\nAvoidance\n***********
15
0.0
1

TEXTBOX
12
747
162
801
************\nQuarantine\n************
15
0.0
1

INPUTBOX
377
294
507
354
pop-SD
0.1
1
0
Number

PLOT
924
11
1478
316
Global Evolution
ticks
S, I, R Global
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"S-Global" 1.0 0 -12087248 true "" "plotxy t / deltaT sum [currentS] of nodes"
"I-Global" 1.0 0 -2674135 true "" "plotxy t / deltaT sum [currentI] of nodes"
"R-Global" 1.0 0 -14454117 true "" "plotxy t / deltaT sum [currentR] of nodes"

MONITOR
791
802
999
847
Total Cost
TotalCost_TrafficLoss_Emissions
2
1
11

SWITCH
364
126
517
159
fixrandomseed?
fixrandomseed?
1
1
-1000

SWITCH
107
353
263
386
One-Infected?
One-Infected?
0
1
-1000

SWITCH
304
594
468
627
Redistribution?
Redistribution?
0
1
-1000

SLIDER
237
163
360
196
random-degree
random-degree
0.1
1
0.4
0.1
1
NIL
HORIZONTAL

MONITOR
1590
549
1839
594
Mobile Population Theoretical
mobilepopulationtheo
2
1
11

MONITOR
1572
85
1820
130
Lost Mobile Population Ratio
mobilepopratio * 100
2
1
11

MONITOR
1572
134
1819
179
Max Lost Mobile Population Ratio
maxmobilepopratio * 100
2
1
11

MONITOR
1591
603
1840
648
Time of Max Mobile Population Ratio
TimeofMaxMobileRatio / deltaT
2
1
11

PLOT
1489
201
1894
496
Traffic Loss Rate
ticks
Max Mobile Pop Ratio
0.0
10.0
0.0
100.0
true
false
"" ""
PENS
"default" 1.0 0 -5298144 true "" "if S1_quarantine? or S2_Avoidance? or S3_RiskCulture? \n[plotxy t / deltaT mobilepopratio * 100]"

TEXTBOX
806
744
956
798
******************\n Total Traffic Loss\n******************
15
0.0
1

TEXTBOX
1622
25
1772
70
*****************\nTraffic Loss Rate\n*****************
12
0.0
1

MONITOR
1120
380
1295
425
Time To Reach All Nodes
TimeToReachAllNodes / deltaT
17
1
11

MONITOR
1304
380
1490
425
Proportion of infected nodes (%)
count nodes with [infected? = 1] / nb_nodes * 100
17
1
11

MONITOR
10
809
219
854
Total Quarantine Cost
QCostTotal
2
1
11

MONITOR
253
805
489
850
Total Avoidance Cost
ACostTotal
2
1
11

MONITOR
928
387
1105
432
Theoric Max I 
maxITheo
2
1
11

CHOOSER
97
79
249
124
Network-Topology
Network-Topology
"Complete" "Chain" "Regular" "Random" "Smallworld"
0

MONITOR
421
12
504
57
NIL
count links
17
1
11

SLIDER
365
163
511
196
Regular-Net-Order
Regular-Net-Order
2
NB_Nodes - 1
2
1
1
NIL
HORIZONTAL

SLIDER
97
126
358
159
SmallWorld_Probability
SmallWorld_Probability
0
1
0.3
0.1
1
NIL
HORIZONTAL

MONITOR
223
203
315
248
R0
beta / alpha
2
1
11

MONITOR
930
434
1064
479
Max Infected prop
MaxITheo / ( sum [currentS + currentI + currentR] of nodes )
2
1
11

MONITOR
12
465
171
510
gs
[gs] of node 0
17
1
11

MONITOR
178
466
335
511
gi
[gi] of node 0
17
1
11

MONITOR
340
468
505
513
gr
[gr] of node 0
17
1
11

INPUTBOX
276
401
437
461
mobility-rate
0.12
1
0
Number

PLOT
1126
475
1470
700
Calibration
ticks
calibS, calibI, calibR
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"calibS" 1.0 0 -14439633 true "" "plotxy ticks calibS"
"calibI" 1.0 0 -14070903 true "" "plotxy ticks calibI"
"caliR" 1.0 0 -2674135 true "" "plotxy ticks calibR"

MONITOR
934
554
1118
599
NIL
calibS
17
1
11

MONITOR
933
602
1119
647
NIL
calibI
17
1
11

MONITOR
934
650
1120
695
NIL
calibR
17
1
11

MONITOR
1184
719
1440
764
NIL
Calibration-Duration
17
1
11

MONITOR
1185
817
1440
862
NIL
Mobile-Population-By-Step
17
1
11

MONITOR
1183
769
1440
814
NIL
Total-Mobile-Population
17
1
11

TEXTBOX
1029
501
1179
546
**************\nCALIBRATION\n**************
12
0.0
1

@#$#@#$#@
Runge-kutta 4th Order Method for the simulation software NetLogo
================================================================

Description
-----------

The program RK4 (*rk4.nls*) is a library for the simulation software NetLogo [1]. This library is performing an iterative method for the approximation of solutions of ordinary differential equations: the Runge-Kutta 4th Order Method [2]. This library works with any number of equations and any number of variables. It is delivered with a demonstration (*rk4.nlogo*) showing the Butterfly attractor [3, 4] studied by Edward Lorenz.

Resources
---------

[Continue reading ](http://flow.chasset.net/netlogo-rk4/)

Author
------

Pierre-Olivier Chasset

Version
-------

1.0.1

License
-------

This program is released under the [GNU Affero General public license](http://www.gnu.org/licenses/agpl.html).

How to cite?
------------

Chasset P.O. (2013). *RK4: Runge-Kutta 4th Order Method for the simulation software NetLogo*. Software, http://flow.chasset.net/netlogo-rk4/.

References
----------

[1] Wilensky U. (1999). NetLogo. Evanston, IL: Center for Connected Learning and Computer-Based Modeling, Northwestern University. Software, http://ccl.northwestern.edu/netlogo/.

[2] Runge C. and Kutta M.W. (around 1900)

[3] Lorenz E.N. (1963). Deterministic nonperiodic flow. Journal of the Atmospheric Sciences 20 (2): 130–141.
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.1.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="Impact enlevement 1 infected par noeud" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>MaxITheo</metric>
    <metric>MaxI</metric>
    <metric>duration / deltaT</metric>
    <metric>TimeofMaxI / deltaT</metric>
    <metric>MaxITheo / TotalPop</metric>
    <metric>MaxI / TotalPop</metric>
    <enumeratedValueSet variable="SInit">
      <value value="100"/>
      <value value="1000"/>
      <value value="10000"/>
      <value value="100000"/>
      <value value="1000000"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="GR" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>MaxI</metric>
    <metric>TimeofMaxI</metric>
    <metric>duration</metric>
    <enumeratedValueSet variable="Regular-Net-Order">
      <value value="4"/>
      <value value="20"/>
      <value value="50"/>
      <value value="70"/>
      <value value="99"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="SM" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>MaxI</metric>
    <metric>TimeofMaxI</metric>
    <metric>duration</metric>
    <enumeratedValueSet variable="SmallWorld_Probability">
      <value value="0.2"/>
      <value value="0.4"/>
      <value value="0.6"/>
      <value value="0.8"/>
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Strategies-10Noeuds" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>MaxI</metric>
    <metric>TimeofMaxI</metric>
    <metric>duration</metric>
    <metric>TotalCost_TrafficLoss_Emissions</metric>
    <metric>RCCostTotal</metric>
    <metric>ACostTotal</metric>
    <metric>QCostTotal</metric>
    <metric>mobilepopratio * 100</metric>
    <metric>maxmobilepopratio * 100</metric>
    <enumeratedValueSet variable="Redistribution?">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <steppedValueSet variable="Quarantine_Threshold" first="0" step="0.2" last="1"/>
    <steppedValueSet variable="RiskCulture_Probability" first="0" step="0.2" last="1"/>
    <steppedValueSet variable="Avoidance_Threshold" first="0" step="0.2" last="1"/>
  </experiment>
  <experiment name="Strategies-30Noeuds-quarantaine" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>MaxI</metric>
    <metric>TimeofMaxI</metric>
    <metric>duration</metric>
    <metric>TotalCost_TrafficLoss_Emissions</metric>
    <metric>RCCostTotal</metric>
    <metric>ACostTotal</metric>
    <metric>QCostTotal</metric>
    <metric>mobilepopratio * 100</metric>
    <metric>maxmobilepopratio * 100</metric>
    <enumeratedValueSet variable="Redistribution?">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <steppedValueSet variable="Quarantine_Threshold" first="2.4E-5" step="5.0E-7" last="2.5E-5"/>
  </experiment>
  <experiment name="Strategies-30Noeuds-évitement" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>MaxI</metric>
    <metric>TimeofMaxI</metric>
    <metric>duration</metric>
    <metric>TotalCost_TrafficLoss_Emissions</metric>
    <metric>RCCostTotal</metric>
    <metric>ACostTotal</metric>
    <metric>QCostTotal</metric>
    <metric>mobilepopratio * 100</metric>
    <metric>maxmobilepopratio * 100</metric>
    <enumeratedValueSet variable="Redistribution?">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <steppedValueSet variable="Avoidance_Threshold" first="0" step="5.0E-6" last="5.0E-5"/>
  </experiment>
  <experiment name="Strategies-30Noeuds-RiskCulture" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>MaxI</metric>
    <metric>TimeofMaxI</metric>
    <metric>duration</metric>
    <metric>TotalCost_TrafficLoss_Emissions</metric>
    <metric>RCCostTotal</metric>
    <metric>ACostTotal</metric>
    <metric>QCostTotal</metric>
    <metric>mobilepopratio * 100</metric>
    <metric>maxmobilepopratio * 100</metric>
    <enumeratedValueSet variable="Redistribution?">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <steppedValueSet variable="RiskCulture_Probability" first="0" step="0.2" last="1"/>
  </experiment>
  <experiment name="Quarantaine-redibfalse" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>MaxI</metric>
    <metric>TimeofMaxI</metric>
    <metric>duration</metric>
    <metric>TotalCost_TrafficLoss_Emissions</metric>
    <metric>QCostTotal</metric>
    <metric>TimeToReachAllNodes / deltaT</metric>
    <enumeratedValueSet variable="NB_Nodes">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SInit">
      <value value="40000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="One-Infected?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Alpha">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Beta">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="S1_Quarantine?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="S2_Avoidance?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="S3_RiskCulture?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="redistribution?">
      <value value="false"/>
    </enumeratedValueSet>
    <steppedValueSet variable="Quarantine_Threshold" first="0" step="0.005" last="0.05"/>
  </experiment>
  <experiment name="Avoidance-RedibFalse" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>MaxI</metric>
    <metric>TimeofMaxI</metric>
    <metric>duration</metric>
    <metric>TotalCost_TrafficLoss_Emissions</metric>
    <metric>ACostTotal</metric>
    <metric>TimeToReachAllNodes / deltaT</metric>
    <enumeratedValueSet variable="NB_Nodes">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SInit">
      <value value="40000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="One-Infected?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Alpha">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Beta">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="S1_Quarantine?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="S2_Avoidance?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="S3_RiskCulture?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="redistribution?">
      <value value="false"/>
    </enumeratedValueSet>
    <steppedValueSet variable="Avoidance_Threshold" first="0" step="0.005" last="0.05"/>
  </experiment>
  <experiment name="RiskCulture" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>MaxI</metric>
    <metric>TimeofMaxI</metric>
    <metric>duration</metric>
    <metric>TotalCost_TrafficLoss_Emissions</metric>
    <metric>RCCostTotal</metric>
    <metric>TimeToReachAllNodes / deltaT</metric>
    <enumeratedValueSet variable="NB_Nodes">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SInit">
      <value value="40000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="One-Infected?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Alpha">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Beta">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="S1_Quarantine?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="S2_Avoidance?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="S3_RiskCulture?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="redistribution?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="RiskCulture_Probability" first="0.005" step="0.05" last="1"/>
  </experiment>
  <experiment name="Quarantaine-redibtrue" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>MaxI</metric>
    <metric>TimeofMaxI</metric>
    <metric>duration</metric>
    <metric>TotalCost_TrafficLoss_Emissions</metric>
    <metric>QCostTotal</metric>
    <metric>TimeToReachAllNodes / deltaT</metric>
    <enumeratedValueSet variable="NB_Nodes">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SInit">
      <value value="40000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="One-Infected?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Alpha">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Beta">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="S1_Quarantine?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="S2_Avoidance?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="S3_RiskCulture?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="redistribution?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="Quarantine_Threshold" first="0" step="0.005" last="0.05"/>
  </experiment>
  <experiment name="Avoidance-RedibTrue" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>MaxI</metric>
    <metric>TimeofMaxI</metric>
    <metric>duration</metric>
    <metric>TotalCost_TrafficLoss_Emissions</metric>
    <metric>ACostTotal</metric>
    <metric>TimeToReachAllNodes / deltaT</metric>
    <enumeratedValueSet variable="NB_Nodes">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SInit">
      <value value="40000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="One-Infected?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Alpha">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Beta">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="S1_Quarantine?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="S2_Avoidance?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="S3_RiskCulture?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="redistribution?">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="Avoidance_Threshold" first="0" step="0.005" last="0.05"/>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

links arn
1.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Polygon -7500403 true true 150 150 105 210 195 210

@#$#@#$#@
0
@#$#@#$#@
