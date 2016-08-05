#edit 160805 2025
#Def werkt voo beide tabellen
# including def function_list(LLSAdepot, LLSAdata): , def function_ALLocate(LLSAdepot, LLSAdata):
#setting: 15,5,8,8,0.9 >>>>> count:  14 Availability result:  0.9000 Investment:  88.0 Tot_EBO:  0.1798 >>>>> 41111 21111
#15,5,8,8,0.95, maar begrensd op Tot_EBO<0.7 en niet Availability. 
#count:  9 Availability result:  0.5667 Investment:  51.0 Tot_EBO:  0.6942 en 31111 20000
# including def function_list(LLSAdepot, LLSAdata): , def function_ALLocate(LLSAdepot, LLSAdata):
#1005_3 ALLeen def funcie depot: def function_depot(LLSAdepot,LLSAdata):
#Inventory_1005_6_combined_ALLocateOK ALLocate werkt op basis van list
#Inventory_1005_6_combined_ALLocateOK_2 list Ok voor depotlevel en s[base_id,x]
#Inventory_1005_6_combined_ALLocateOK_2 LLSAdata[0,depotlevel,s[base_id,x],x,P.EBO] = LLSAdepot[depotlevel,x,P_depot.EBO] weer teruggezet naar regel 115
#Inventory_1005_9_lrunames setting: 15,5,8,8,0.9 >>>>> count:  14 Availability result:  0.9000 Investment:  88.0 Tot_EBO:  0.1798 >>>>> 41111 21111
#Inventory_1005_9_lrunames 
#Inventory_1005_9toch aanpassing function_depot in ALLocation voor base_id =0 en base_id>0
#setting: 15,5,8,8,0.9 >>>>> count:  14 Availability result:  0.9000 Investment:  88.0 Tot_EBO:  0.1798 >>>>> 41111 21111
#Inventory_1005_10_150706_1801 30,5,12,12,0.999 geen problemen Error program Availability not increasing" 73322 42222
#1152 3168 11 11 11 11 11 en count:  29 Availability result:  0.9969 Investment:  181.0 Tot_EBO:  0.0044
#test08_1_werkt count:  14 Availability result:  0.8100 Investment:  88.0 Tot_EBO:  0.1798 41111 21111
# waarom s[base_id,x] in function_depot(LLSAdepot,LLSAdata,base_id,depotlevel,s[base_id,x],x): s[base_id,x] verwijderd.
#function_depot(LLSAdepot,LLSAdata,base_id,s[0,x],x)   verplaatst in: function_list_depot
#vervang woord depotlevel door s[0,x]
#150717 test14-met availability ALLe parameters bij return verwijderd.
#150717 nog steeds probleem negative waarden voor DeltaEBO
#test17 error met nan 150,5,18,18,0LLSAdata[base_id,s[0,x],s[base_id,x],x,P.EBO] > LRUdata[x,P_lru.Qty_LRU]: avail.99999
#150723_1147. Lijkt volledig te werken met 15,5,8,8,0.9
#150725 0.23 error 150,5,18,18,0.99999 in availability niet increasing voor A = 0.9953 bij depot wisseling. Lijkt in de bases berekening te zitten als gevolg van nbinom.pmf
#150725 " LLSAdata[base_id,s[0,x],s[base_id,x],x,P.VBO] in function-base gedeactiveerd.
#150726 ingebouwd vergelijk var en mu voor binom. Als eerste 3 digits gelijk dan poisson ipv nbinom.  
#150729 " ALLE BASES POISSON
#1507290 Aangebracht toets tbv van availability LLSAdata[base_id,s[0,x],s[base_id,x],x,P.EBO] > LRUdata[x,P_lru.Qty_LRU]: 
#150801  function_base(LLSAdepot,LLSAdata,base_id,s,x) teveel in function_sparestart. Uitgevinkt.
#150801 Fr depot calculatie van _list_depot naar _depot verschoven. Daardoor overbodig in _sparestart]
#150802 Verwijderd list_base
#150806 Investment van LRUdata naar depot omdat anders base_id = 0 moet worden toegekend.
#Thijs. 
#Louis problemen . Negatieve waarden Delta EBO voor oversized depot of baselevel.
#150826 Toegevoegd MTTR_MSI
#150826 Years omgezet in hours

import time
import numpy as np
from scipy.spatial import ConvexHull
from scipy.stats import poisson
from scipy.stats import nbinom
from scipy.integrate import quad
from numpy import ones,zeros,unravel_index,empty,nan,isnan,subtract,ma,diff,newaxis
import matplotlib.pyplot as plt
import pylab as pl
import sys 
import pandas as pd
from enum import IntEnum
from prettytable import PrettyTable, ALL            
import csv
from numpy import genfromtxt
import pprint
import itertools

start_time = time.time()
#
class P_lru(IntEnum):   
    FR = 2 
    MTTR_MSI =3
    MTTR_base = 7
    Qty_LRU = 8
    RHO = 9
    MLDT = 10
    sparebase = 19
class P_lrudepot(IntEnum): 
    MTTR_depot = 1 
    sparedepot = 2
    Investment = 3 
class P_depot(IntEnum):
    FR = 1
    P_s = 4
    ROS = 5
    EBO = 6
    Qty_LRU = 8  
    Mu = 13
    VBO = 14
    Var = 15
class P(IntEnum):
    DeltaEBO= 0
    P_s = 4
    ROS = 5
    EBO = 6
    Mu = 13
    VBO = 14
    Var = 15
    r_neg = 16
    p_neg =17
class S(IntEnum):
    max_count = 0
    max_bases = 1
    max_stocklevel_depot = 2
    max_stocklevel_bases = 3
    max_availability = 4
    MPMT = 5
    MTBPM = 6
    MBaseLDT = 7
    QTY_subsys = 8


#LRUdata = genfromtxt('./data/metric01_new.csv', delimiter=',', skip_header=2)
contents = []

def to_float(x):
    try:
        x = float(x)
    except ValueError:
        x = 0
    return x

#https://docs.python.org/2/library/csv.html
with open('./data/LRU26.csv', 'rb') as f:
    reader = csv.reader(f, delimiter=",")
    for row in reader:
        contents.append(row)
        
#verwijder de eerste twee regels
contents = contents[2:]
# QTY LRU > 1. waarom?
LRUdata = np.array(map(lambda x: map(lambda a: to_float(a), [0] + x[1:]),contents))
Setdata = genfromtxt('./data/set_metric26.csv', delimiter=',')

AvailabilityTarget = float(Setdata[S.max_availability])

max_count = int(Setdata[S.max_count])
#minimum bases = 2, nl.  depot 0 en base 1, dus : max_bases = 2
#minimum stocklevel = 1 
max_stocklevel_depot = int(Setdata[S.max_stocklevel_depot])
max_stocklevel_bases = int(Setdata[S.max_stocklevel_bases])
QTY_subsys = int(Setdata[S.QTY_subsys])
MPMT = float(Setdata[S.MPMT])
MTBPM = float(Setdata[S.MTBPM])
MBaseLDT = float(Setdata[S.MBaseLDT])



#minimum  ma_lru = 1
max_lru_base_lines=LRUdata.shape[0]
max_param=LRUdata.shape[1]

depot_contents = []

with open('./data/depot26.csv', 'rb') as f:
    reader = csv.reader(f, delimiter=",")
    for row in reader:
        depot_contents.append(row)
        
#verwijder de eerste twee regels
depot_contents = depot_contents[2:]

#maak een lijst met lru names door de eerste kolom van de regels te nemen
unique_lru_names = map(lambda x: x[0], depot_contents)
max_lru = len(unique_lru_names)
LRUdepot = np.array(map(lambda x: map(lambda a: to_float(a), [0] + x[1:]),depot_contents))
max_param_depot=LRUdepot.shape[1]

#maak een dictionary
lrus_base_ids = {}

for unique_lru_name in unique_lru_names:
    base_id_rows = []
    for rownr, rowval in enumerate(contents):
        if rowval[0] == unique_lru_name:
            base_id_rows.append(rownr)
    lrus_base_ids[unique_lru_name] = base_id_rows
max_bases = (max(map(lambda x: len(x), lrus_base_ids.values())))
highest_base_id = max(map(lambda x: int(x[1]), contents))

print 'Startwaarden'
print 'max_count: ',max_count,'max_bases: ',max_bases,'max_stocklevel_depot:',max_stocklevel_depot,\
'max_stocklevel_bases:',max_stocklevel_bases,'QTY LRUs:',max_lru, 'AvailabilityTarget',AvailabilityTarget

LRUdata = ones((max_lru_base_lines,max_param)) * LRUdata

LRUdepot = ones((max_lru,max_param_depot)) * LRUdepot
LSAdata = ones((max_bases + 1,max_stocklevel_depot,max_stocklevel_bases,max_lru,max_param)) * nan
print LSAdata.shape
LSAdepot = ones((max_stocklevel_depot,max_lru,max_param)) * nan
Totdepot_BO =zeros((max_lru))  
Totbases_BO =zeros((max_lru))
Totbases_FR =zeros((max_lru))
Investment = zeros((max_count))

DeltaCost = zeros((max_bases + 1,max_lru))
Tot_EBO = zeros((max_count))
Tottotbases_BO = zeros((max_count))
Tottotbases_FR = zeros((max_count))

a = zeros((max_lru))
b = zeros((max_lru))
c = zeros((max_lru))
d = zeros((max_lru))
Tot_spares = 0
s = zeros((max_bases + 1,max_lru))
Availability_tot = ones((max_lru))
Availability_sys = ones((max_count))
Availability_base = ones((max_bases + 1,max_lru))
Unavailability_base= zeros((max_bases + 1,max_lru))


base_id=0
depotlevel=0
count = 0
max_s = 0
x = 0

function_base_count = 0
function_depot_count = 0


def get_lru_name(lru_index):
    return unique_lru_names[lru_index]

def get_lru_baserownr_by_baseindex(lru_index, base_index):
    base_index = base_index - 1
    lru_name = get_lru_name(lru_index)
    base_id_rows = lrus_base_ids[lru_name]
    
    return base_id_rows[base_index] 

def get_lru_base_name(lru_index, base_index):
    rownr = get_lru_baserownr_by_baseindex(lru_index, base_index)
    return contents[rownr][1]

def get_lru_baserownr_by_basename(lru_index, base_name):
    lru_name = get_lru_name(lru_index)
    base_id_rows = lrus_base_ids[lru_name]
    
    for base_id_row in base_id_rows:
        if contents[base_id_row][1] == str(base_name):
            return base_id_row
        
def get_lru_baseindex_by_basename(lru_index, base_name):
    lru_name = get_lru_name(lru_index)
    base_id_rows = lrus_base_ids[lru_name]
    
    for index, base_id_row in enumerate(base_id_rows):
        if contents[base_id_row][1] == str(base_name):
            return index + 1
 
def get_lru_data(lru_index, base_index, param):  
    return LRUdata[get_lru_baserownr_by_baseindex(lru_index, base_index), param]   

# function A. depot formules
def function_depot(LLSAdepot,LLSAdata,s,x):
#     print ' ' 
    
    global function_depot_count
    function_depot_count = function_depot_count+1
#     print 'start function_depot'
#     print 'LRU id: ',x
    
    
    ##Mu in pipeline at depot. EQ 13
    LLSAdepot[s[0,x],x,P_depot.Mu] = LLSAdepot[s[0,x],x,P_depot.FR] * LRUdepot[x,P_lrudepot.MTTR_depot]     
    # The considered stocklevel
    if s[0,x] == 0:            
        #Demand probability P (depot)
        LLSAdepot[s[0,x],x,P_depot.P_s] = poisson.pmf(s[0,x], LLSAdepot[s[0,x],x,P_depot.Mu])
        # ROS (depot)
        LLSAdepot[s[0,x],x,P_depot.ROS] = 1  
        # EBO (depot)       
        LLSAdepot[s[0,x],x,P_depot.EBO] = LLSAdepot[s[0,x],x,P_depot.Mu]
        # VBO (depot)
        LLSAdepot[s[0,x],x,P_depot.VBO] = LLSAdepot[s[0,x],x,P_depot.EBO]
#         print ''
#         print 'A0: Depot------------------------------'  
#         print 'A0: LRU id: ',x ,'s[0,x]: ',s[0,x]
#         print 'A0: Demand rate FR at depot: ','%.3f' % LLSAdepot[s[0,x],x,P_depot.FR]
#         print 'A0: #Average_Qty_LRU_ in pipeline at depot. EQ 17at depot: ','%.3f' % LLSAdepot[s[0,x],x,P_depot.Mu]
#         print 'A0: Demand probability P at depot: ','%.3f' % LLSAdepot[s[0,x],x,P_depot.P_s]
#         print 'A0: ROS (Risk Out of Stock) at depot: ','%.3f' % LLSAdepot[s[0,x],x,P_depot.ROS]
#         print 'A0: EBO (Expected Back Order) at depot: ','%.3f' % LLSAdepot[s[0,x],x,P_depot.EBO]
#         print 'A0: VBO (Expected Variance of Back Order) at depot: ','%.3f' % LLSAdepot[s[0,x],x,P_depot.VBO]
#         print ''  
#    
    else:   
        #ROS (depot)
        LLSAdepot[s[0,x],x,P_depot.ROS] = LLSAdepot[s[0,x]-1,x,P_depot.ROS]\
            - LLSAdepot[s[0,x]-1,x,P_depot.P_s] 
        # VBO (depot)
        LLSAdepot[s[0,x],x,P_depot.VBO] = LLSAdepot[s[0,x]-1,x,P_depot.VBO]\
           - (LLSAdepot[s[0,x]-1,x,P_depot.EBO] - LLSAdepot[s[0,x],x,P_depot.ROS])\
           - LLSAdepot[s[0,x]-1,x,P_depot.EBO] - (LLSAdepot[s[0,x]-1,x,P_depot.EBO]\
           - LLSAdepot[s[0,x],x,P_depot.ROS]) **2 + LLSAdepot[s[0,x]-1,x,P_depot.EBO]**2 
        # EBO (depot)
        LLSAdepot[s[0,x],x,P_depot.EBO] = LLSAdepot[s[0,x]-1,x,P_depot.EBO]\
            - LLSAdepot[s[0,x],x,P_depot.ROS] 
        #Demand probability P at depot
        LLSAdepot[s[0,x],x,P_depot.P_s] = poisson.pmf(s[0,x], LLSAdepot[s[0,x],x,P_depot.Mu])
#         print ''
#         print 'A1: Depot------------------------------'  
#         print 'A1: LRU id: ',x ,'s[0,x]: ',s[0,x]
#         print 'A1: Demand rate FR at depot: ','%.3f' % LLSAdepot[s[0,x],x,P_depot.FR]
#         print 'A1: #Average_Qty_LRU_ in pipeline at depot. EQ 17at depot: ','%.3f' % LLSAdepot[s[0,x],x,P_depot.Mu]
#         print 'A1: Demand probability P at depot: ','%.3f' % LLSAdepot[s[0,x],x,P_depot.P_s]
#         print 'A1: ROS (Risk Out of Stock) at depot: ','%.3f' % LLSAdepot[s[0,x],x,P_depot.ROS]
#         print 'A1: EBO (Expected Back Order) at depot: ','%.3f' % LLSAdepot[s[0,x],x,P_depot.EBO]
#         print 'A1: VBO (Expected Variance of Back Order) at depot: ','%.3f' % LLSAdepot[s[0,x],x,P_depot.VBO]
#         print ''   
    

def function_base(LLSAdepot, LLSAdata,base_id,s,x): 
#     print ' '

    global function_base_count
    function_base_count = function_base_count+1
#     print 'start function base'
#     print 'LRU id: ',x ,'base_id:',base_id
    #LLSAdata[0,s[0,x],s[base_id,x],x,P.EBO] = LLSAdepot[s[0,x],x,P_depot.EBO]

    #Mu in pipeline at base. EQ 17        
    
    LLSAdata[base_id,s[0,x],s[base_id,x],x,P.Mu]=\
    get_lru_data(x,base_id,P_lru.Qty_LRU) * get_lru_data(x,base_id,P_lru.FR) * (get_lru_data(x,base_id,P_lru.RHO) * get_lru_data(x,base_id,P_lru.MTTR_base)+\
    (1-get_lru_data(x,base_id,P_lru.RHO))*(get_lru_data(x,base_id,P_lru.MLDT) + LLSAdepot[s[0,x],x,P.EBO] / LLSAdepot[s[0,x],x,P_depot.FR]))           
    if s[base_id,x] == 0:                 
        # P_s(s[base_id,x] = 0)
        LLSAdata[base_id,s[0,x],s[base_id,x],x,P.P_s] = poisson.pmf(s[base_id,x],LLSAdata[base_id,s[0,x],s[base_id,x],x,P.Mu])
        # ROS(s[base_id,x] = 0)
        LLSAdata[base_id,s[0,x],s[base_id,x],x,P.ROS] = 1 
        # EBO (s[base_id,x] = 0)  =  Mu in pipeline at base. EQ 17 
        LLSAdata[base_id,s[0,x],s[base_id,x],x,P.EBO] = LLSAdata[base_id,s[0,x],s[base_id,x],x,P.Mu]  
        # VBO (baselelev = 0) = EBO (s[base_id,x] = 0)
        LLSAdata[base_id,s[0,x],s[base_id,x],x,P.VBO] = LLSAdata[base_id,s[0,x],s[base_id,x],x,P.EBO] 
        # VAR (baselelev = 0) EQ20
        LLSAdata[base_id,s[0,x],s[base_id,x],x,P.Var] = \
        LLSAdata[base_id,s[0,x],s[base_id,x],x,P.Mu] \
        - (get_lru_data(x,base_id,P_lru.FR)*get_lru_data(x,base_id,P_lru.Qty_LRU)\
           * (1-get_lru_data(x,base_id,P_lru.RHO))/LLSAdepot[s[0,x],x,P_depot.FR])**2
        
#         print 'B0: baselevel ','%d' % s[base_id,x],'depotlevel',s[0,x]
#         print 'B0: P.Mu','%.3f' % LLSAdata[base_id,s[0,x],s[base_id,x],x,P.Mu]
# #         print 'Var','%.3f' % LLSAdata[base_id,s[0,x],s[base_id,x],x,P.Var]
#         print 'B0: P.EBO baselevel = 0  ',LLSAdata[base_id,s[0,x],s[base_id,x],x,P.EBO] 
#            
#         print 'B0: P.ROS baselevel = 0  ',LLSAdata[base_id,s[0,x],s[base_id,x],x,P.ROS]
#         print 'B0: P.P_s baselevel = 0  ',LLSAdata[base_id,s[0,x],s[base_id,x],x,P.P_s]    
#         print " "             
    else:                                                               
        LLSAdata[base_id,s[0,x],s[base_id,x],x,P.Var] = \
        LLSAdata[base_id,s[0,x],s[base_id,x],x,P.Mu] \
        - (get_lru_data(x,base_id,P_lru.FR)*get_lru_data(x,base_id,P_lru.Qty_LRU) * (1-get_lru_data(x,base_id,P_lru.RHO))/LLSAdepot[s[0,x],x,P_depot.FR])**2\
        * (LLSAdepot[s[0,x],x,P.EBO] - LLSAdepot[s[0,x],x,P.VBO])   
        # r_neg EQ18
        LLSAdata[base_id,s[0,x],s[base_id,x],x,P.r_neg] = LLSAdata[base_id,s[0,x],s[base_id,x],x,P.Mu]**2\
        /(LLSAdata[base_id,s[0,x],s[base_id,x],x,P.Var]-LLSAdata[base_id,s[0,x],s[base_id,x],x,P.Mu])
        # p_neg EQ19
        LLSAdata[base_id,s[0,x],s[base_id,x],x,P.p_neg] = LLSAdata[base_id,s[0,x],s[base_id,x],x,P.Mu]/\
        LLSAdata[base_id,s[0,x],s[base_id,x],x,P.Var] 
#             print 'voor ROS' , LLSAdata[base_id,s[0,x],s[base_id,x]-1,x,P.ROS], 'P_s', LLSAdata[base_id,s[0,x],s[base_id,x]-1,x,P.P_s] 
        LLSAdata[base_id,s[0,x],s[base_id,x],x,P.ROS] = LLSAdata[base_id,s[0,x],s[base_id,x]-1,x,P.ROS] - LLSAdata[base_id,s[0,x],s[base_id,x]-1,x,P.P_s]  
#             print ' na ROS' , LLSAdata[base_id,s[0,x],s[base_id,x]-1,x,P.ROS], 'P_s', LLSAdata[base_id,s[0,x],s[base_id,x]-1,x,P.P_s]            
        LLSAdata[base_id,s[0,x],s[base_id,x],x,P.VBO] = LLSAdata[base_id,s[0,x],s[base_id,x]-1,x,P.VBO] - LLSAdata[base_id,s[0,x],s[base_id,x],x,P.ROS]\
                    - LLSAdata[base_id,s[0,x],s[base_id,x]-1,x,P.P_s] - LLSAdata[base_id,s[0,x],s[base_id,x]-1,x,P.EBO]\
                    - (LLSAdata[base_id,s[0,x],s[base_id,x],x,P.ROS] - LLSAdata[base_id,s[0,x],s[base_id,x]-1,x,P.P_s]) **2\
                     + LLSAdata[base_id,s[0,x],s[base_id,x]-1,x,P.EBO]**2  
        LLSAdata[base_id,s[0,x],s[base_id,x],x,P.EBO] = LLSAdata[base_id,s[0,x],s[base_id,x]-1,x,P.EBO] - LLSAdata[base_id,s[0,x],s[base_id,x],x,P.ROS]

       
        if s[0,x] == 0:
#             print 'poisson'
#             LLSAdata[base_id,s[0,x],s[base_id,x],x,P.P_s] = nbinom.pmf(s[base_id,x],LLSAdata[base_id,s[0,x],s[base_id,x],x,P.r_neg],LLSAdata[base_id,s[0,x],s[base_id,x],x,P.p_neg])         
            LLSAdata[base_id,s[0,x],s[base_id,x],x,P.P_s] = poisson.pmf(s[base_id,x], LLSAdata[base_id,s[0,x],s[base_id,x],x,P.Mu])
             
        else:
#             LLSAdata[base_id,s[0,x],s[base_id,x],x,P.P_s] = nbinom.pmf(s[base_id,x],LLSAdata[base_id,s[0,x],s[base_id,x],x,P.r_neg],LLSAdata[base_id,s[0,x],s[base_id,x],x,P.p_neg])         
            LLSAdata[base_id,s[0,x],s[base_id,x],x,P.P_s] = poisson.pmf(s[base_id,x], LLSAdata[base_id,s[0,x],s[base_id,x],x,P.Mu])
#             print 'nbinom'
#         
#             print' '
#         print 'B1: baselevel','%d' % s[base_id,x],' at depotlevel','%d' % s[0,x],'for base_id_',base_id
#         print 'B1: P.Mu','%.5f' % LLSAdata[base_id,s[0,x],s[base_id,x],x,P.Mu]
#         print 'Var','%.5f' % LLSAdata[base_id,s[0,x],s[base_id,x],x,P.Var]
#         print 'p_neg: ','%.3f' % LLSAdata[base_id,s[0,x],s[base_id,x],x,P.p_neg]
#         print 'r_neg: ','%.3f' % LLSAdata[base_id,s[0,x],s[base_id,x],x,P.r_neg]
#         print 'B1: P.EBO baselevel  ',LLSAdata[base_id,s[0,x],s[base_id,x],x,P.EBO]
#         print 'B1: P.ROS baselevel',LLSAdata[base_id,s[0,x],s[base_id,x],x,P.ROS]
#         print 'B1: P.P_s baselevel',LLSAdata[base_id,s[0,x],s[base_id,x],x,P.P_s]
# #         print 'B1: baselevel -1: ' , s[base_id,x]-1, ' at depotlevel','%d' % s[0,x]
#                
#         print 'B1: P.EBO baselevel-1',LLSAdata[base_id,s[0,x],s[base_id,x]-1,x,P.EBO]
#         print 'B1: P.ROS baselevel-1]',LLSAdata[base_id,s[0,x],s[base_id,x]-1,x,P.ROS]  
#         print 'B1: P.P_s baselevel-1]',LLSAdata[base_id,s[0,x],s[base_id,x]-1,x,P.P_s] 
#                           
#         print 'B1 end'
    return

def function_sparestart(LLSAdepot, LLSAdata): 
#     print "Start function-sparestart"

    w = PrettyTable(["LRU_id","s[0]","s[1]","s[2]","s[3]","s[4]","Delta_BO","TOT_BO","Availability","P.EBO"])  
    w.max_width = 10
    w.padding_width = 1 # One space between column edges and contents (default)
    w.hrules = ALL
    
    base_id=1
    
    for x in range(max_lru):
        lru_max_bases = len(lrus_base_ids[unique_lru_names[x]])
        # Tbv gaten doe en test met depot 0 en base 1 in sparestart. Zichtbaar in ALLocate neg. waarden voor deltaEBO . Daarom +1. zie ook einde script.
        for s[0,x] in range (0,int(LRUdepot[x,P_lrudepot.sparedepot])+1):
            Totbases_BO[x] = 0
            LLSAdepot[s[0,x],x,P_depot.FR] = 0
           
            for base_id in range (1,lru_max_bases+1):
                #L1: Demand rate FR at depot EQ12
                LLSAdepot[s[0,x],x,P_depot.FR] += get_lru_data(x, base_id,P_lru.Qty_LRU)*get_lru_data(x, base_id,P_lru.FR)*\
                (1-get_lru_data(x, base_id,P_lru.RHO))         
            #L2: Vaststellen wat startwaarde is voor ieder depotlevel = s[0,x] en voor baselevel = 0, iedere x
            
            for base_id in range (1,lru_max_bases+1): 
                #baselevel =0
               
                s[base_id,x]=0 
                function_depot(LLSAdepot,LLSAdata,s,x)
                function_base(LLSAdepot,LLSAdata,base_id,s,x)
                Totbases_BO[x] += LLSAdata[base_id,s[0,x],s[base_id,x],x,P.EBO]
            Totdepot_BO[x] =  Totbases_BO[x] + LLSAdepot[s[0,x],x,P_depot.EBO]
#             print 'LRU_id', x, 's[0,x]', s[0,x].
#             print 'Totbases_BO', Totbases_BO[x]
#             print 'Totdepot_BO',Totdepot_BO[x]
# 
#             print 'P.EBO',LLSAdata[base_id,s[0,x],s[base_id,x],x,P.EBO]
#             print 'P_depot',LLSAdepot[s[0,x],x,P_depot.EBO]
#             print ''
            
            #L3:  Tottot_BO voor ieder depotlevel s[0,x] =(P.sparedepot,x) en voor baselevel (P.sparebase) op basis van de voorgaande startwaarde.
            # bijzondere volgorde van baselevel en  base_id om juiste lijst outputformaat  te verkrijgen
            for base_id in range (1,lru_max_bases+1):
                for s[base_id,x]  in range (0,int(get_lru_data(x, base_id,P_lru.sparebase))+1):   
                    function_base(LLSAdepot,LLSAdata,base_id,s,x)      
                    if s[base_id,x] == 0:  # startwaarden voor ieder depotlevel = s[0,x] en voor baselevel = 0, iedere x
                        DeltaEBO[base_id,x]=0
                        Tottot_BO = Totdepot_BO[x]
                    else:
                        #The decrease in EBOfor baselevel = s[base_id,x]
                        DeltaEBO[base_id,x]=LLSAdata[base_id,s[0,x],s[base_id,x]-1,x,P.EBO]-LLSAdata[base_id,s[0,x],s[base_id,x],x,P.EBO]
                    Tottot_BO = Tottot_BO - DeltaEBO[base_id,x]
                         
#                     print 'function_sparestart'
#                  
#                     print 'LRU id: ',x ,'base_id:',base_id,
#                     print 'baselevel','%d' % s[base_id,x],' at depotlevel','%d' % s[0,x]
#                     print 'LLSAdata[base_id,s[0,x],s[base_id,x],x,P.EBO]','%.3f' % LLSAdata[base_id,s[0,x],s[base_id,x],x,P.EBO]
#                     w.add_row([x,'%d' % s[0,x],'%d' % s[1,x],'%d' % s[2,x],'%d' % s[3,x],'%d' % s[4,x],\
#                                '%.3f' % DeltaEBO[base_id,x],'%.3f' % Tottot_BO,'%.3f' % Availability_base[base_id,x],\
#                                '%.3f' % LLSAdata[base_id,s[0,x],s[base_id,x],x,P.EBO]])  
        
        # Tbv gaten doe en test met depot 0 en base 1 in sparestart. Zichtbaar in ALLocate neg. waarden voor deltaEBO
#         s[0,x] = s[0,x]-1    
#         for base_id in range (0,lru_max_bases+1 + 1): 
#             print s[base_id,x]  
#             print'yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy'                
       
#     print w 
#     print 'end function_sparestart'
    return Totbases_BO
#function_sparestart(LSAdepot, LSAdata)     

#function D. vaststellen Tot EBO with s[base_id,x] voor ALLe bases bij ieder depotlevel

#     print 'E1: count', count,'LRU id: ',x ,'base_id:',base_id
#     print 'E1: baselevel = ' ,'%d' % s[base_id,x]
#     print ''
#     print '-----Start function_list_base'
    s[base_id,x] = s[base_id,x] +1
#     print 'E1d: base_id:',base_id,'baselevel +1:', '%d' % s[base_id,x]
    function_base(LLSAdepot,LLSAdata,base_id,s,x)   
    d[x] = LLSAdata[base_id,s[0,x],s[base_id,x],x,P.EBO]

    s[base_id,x] = s[base_id,x] -1
#     print 'E1d: dddd','%.3f' % d[x]
#     print '-----E1: End function_list_base'  
#     print '-----E1: Start function_list_base'
#     print 'E1c: base_id',base_id,'baselevel:', '%d' % s[base_id,x] 
    function_base(LLSAdepot,LLSAdata,base_id,s,x)   
    c[x] = LLSAdata[base_id,s[0,x],s[base_id,x],x,P.EBO]
#     print 'E1c: cccc','%.3f' % c[x] 
#     print '-----E1: End function_list_base'
    DeltaEBO[base_id,x]  = (c[x] - d[x])
    DeltaCost[base_id,x] = DeltaEBO[base_id,x]/(LRUdepot[x,P_lrudepot.Investment]*get_lru_data(x,base_id,P_lru.Qty_LRU))                     

#     print 'E1: DeltaEBO[base_id,x] c[x] - d[x]','%.3f' % DeltaEBO[base_id,x]
# #initieer depot>>> depotlevel = 0, bases  >>> baselevel = 0
# def function_init(LLSAdepot,LLSAdata):
#     for x in range(max_lru): 
#         lru_max_bases = len(lrus_base_ids[unique_lru_names[x]])  
#         for base_id in range (0,lru_max_bases+1 + 1):              
#             s[base_id,x] = 0 
# #             print 'E: Start'  
#             function_list_depot(LSAdepot,LSAdata,base_id,s,x,lru_max_bases)          

#stock result table
def function_stockresult(LSAdepot,LSAdata,base_id,s,x): 
    ww = PrettyTable(["LRU_id"] + map(lambda x: "Base_id_%d" % x, range(highest_base_id + 1)))
      
    for lru_id in range(max_lru):   
        data = [unique_lru_names[lru_id], int(s[0, lru_id])]
        for base_id in range(1, highest_base_id + 1):
            base_index = get_lru_baseindex_by_basename(lru_id, base_id)
            if base_index is not None:
                data.append(int(s[base_index, lru_id]))
            else:
                data.append(None)    
        ww.add_row(data) 
      
    LRUdepot[x,P_lrudepot.Investment]   
    ww.add_column("Tot_Spares", map(lambda x: "%d" % x, s.sum(0)))  
    ww.max_width = 30
    ww.padding_width = 1 # One space between column edges and contents (default)
    ww.hrules = ALL
    print ww

def function_detailedresult(LLSAdepot,LLSAdata,base_id,s,x,lru_max_bases):      
    print '-------------------------------------------------------------------------'
    TotStock = 0
    FR_sys = 0
    MTTR_sys = 0
    Unavail_sys = 0
    BO_singlesubsys_totMSIx_siBase = zeros((max_bases + 1,max_lru))
    Ainh_singlesubsys_totMSIx_siBase = ones((max_bases + 1,max_lru))
    Asup_singlesubsys_totMSIx_siBase = ones((max_bases + 1,max_lru))
    Aop_singlesubsys_totMSIx_siBase = ones((max_bases + 1,max_lru))
    QTY_singlesubsys_totMSIx_siBase = zeros((max_bases + 1,max_lru))
    MTBF_singlesubsys_totMSIx_siBase = zeros((max_bases + 1,max_lru))
    FR_singlesubsys_totMSIx_siBase = zeros((max_bases + 1,max_lru))
    MWT_singlesubsys_totMSIx_siBase = zeros((max_bases + 1,max_lru))
    MTTR_singlesubsys_totMSIx_siBase = zeros((max_bases + 1,max_lru))
    ROS_singlesubsys_totMSIx_siBase = zeros((max_bases + 1,max_lru))
   
    '---------------------------------------------------'
    
    BO_singlesubsys_siBase = zeros((max_bases + 1))
    FR_singlesubsys_siBase = zeros((max_bases + 1))
    
    MTBF_singlesubsys_siBase = zeros((max_bases + 1))
    UnAinh_singlesubsys_siBase = zeros((max_bases + 1))
    MTTR_singlesubsys_siBase = zeros((max_bases + 1))
    Ainh_singlesubsys_siBase = ones((max_bases + 1))
    
    fPM_singlesubsys_siBase = zeros((max_bases + 1))
    MTBPM_singlesubsys_siBase = zeros((max_bases + 1))
    MPMT_singlesubsys_siBase = zeros((max_bases + 1))
    
    
    MTBAM_singlesubsys_siBase = zeros((max_bases + 1))
    UnAach_singlesubsys_siBase = zeros((max_bases + 1))
    MAMT_singlesubsys_siBase = zeros((max_bases + 1))
    Aach_singlesubsys_siBase = ones((max_bases + 1))
    
    Asup_singlesubsys_siBase = ones((max_bases + 1))
    MWT_singlesubsys_siBase = zeros((max_bases + 1))
    MTBDemand_singlesubsys_siBase = zeros((max_bases + 1))
    UnSO_singlesubsys_siBase = zeros((max_bases + 1))
    ROS_singlesubsys_siBase = zeros((max_bases + 1))
    
    MBaseLDT_singlesubsys_siBase = zeros((max_bases + 1))
    MTBOM_singlesubsys_siBase = zeros((max_bases + 1))
    UnAop_singlesubsys_siBase = zeros((max_bases + 1))
    Aop_singlesubsys_siBase = ones((max_bases + 1))
    MDT_singlesubsys_siBase = zeros((max_bases + 1)) 
    Aop_singlesubsys_siBase = ones((max_bases + 1))
    
    QTY_singlesubsys_siBase = zeros((max_bases + 1))
          
    Ainh_allMSI_allBase = 1
    Asup_allMSI_allBase = 1
    Aop_allMSI_allBase =1
    
    Ainh_allMSI_allBase = 1
    Asup_allMSI_allBase = 1
    Aop_allMSI_allBase = 1
    BO_allMSI_allBase = 0
    FR_allMSI_allBase = 0
    MWT_allMSI_allBase = 0
      
    print '------------------Data per MSIx per subsystem for each single Base--------------------'
    w = PrettyTable(["Base_id","LRU_id","Spares","BO","FR","MTBF","MWT","MTTR","Ai","Asup","Aop","QTY of MSIx per subsystem"]) 
    w.max_width = 10
    w.padding_width = 1 # One space between column edges and contents (default)
    w.hrules = ALL   
    
    for base_id in range (1,lru_max_bases+1): 
        for x in range(max_lru): 
            lru_max_bases = len(lrus_base_ids[unique_lru_names[x]])  
            
            QTY_singlesubsys_totMSIx_siBase[base_id,x] = get_lru_data(x,base_id,P_lru.Qty_LRU)/QTY_subsys        
            BO_singlesubsys_totMSIx_siBase[base_id,x] = LLSAdata[base_id,s[0,x],s[base_id,x],x,P.EBO]/QTY_subsys
            
            
            FR_singlesubsys_totMSIx_siBase[base_id,x] = get_lru_data(x,base_id,P_lru.FR)*QTY_singlesubsys_totMSIx_siBase[base_id,x]
            MTBF_singlesubsys_totMSIx_siBase[base_id,x] =1/FR_singlesubsys_totMSIx_siBase[base_id,x]      
            MWT_singlesubsys_totMSIx_siBase[base_id,x] = BO_singlesubsys_totMSIx_siBase[base_id,x]/FR_singlesubsys_totMSIx_siBase[base_id,x]
            MTTR_singlesubsys_totMSIx_siBase[base_id,x] = get_lru_data(x,base_id,P_lru.MTTR_MSI)
            Ainh_singlesubsys_totMSIx_siBase[base_id,x]  = MTBF_singlesubsys_totMSIx_siBase[base_id,x]/(MTBF_singlesubsys_totMSIx_siBase[base_id,x]\
                                                                                                        +MTTR_singlesubsys_totMSIx_siBase[base_id,x])
            Asup_singlesubsys_totMSIx_siBase[base_id,x] = 1-BO_singlesubsys_totMSIx_siBase[base_id,x]
            
            'moet anders berekend worden via UnA'
            Aop_singlesubsys_totMSIx_siBase[base_id,x] = Ainh_singlesubsys_totMSIx_siBase[base_id,x]*Asup_singlesubsys_totMSIx_siBase[base_id,x]
            
            w.add_row([base_id,x,s[base_id,x],\
                        '%.5f' % BO_singlesubsys_totMSIx_siBase[base_id,x],\
                        '%.7f' % FR_singlesubsys_totMSIx_siBase[base_id,x],\
                        '%.1f' % MTBF_singlesubsys_totMSIx_siBase[base_id,x],\
                        '%.1f' % MWT_singlesubsys_totMSIx_siBase[base_id,x],\
                        '%.3f' % get_lru_data(x,base_id,P_lru.MTTR_MSI),\
                        '%.5f' % Ainh_singlesubsys_totMSIx_siBase[base_id,x],\
                        '%.5f' % Asup_singlesubsys_totMSIx_siBase[base_id,x],\
                        '%.5f' % Aop_singlesubsys_totMSIx_siBase[base_id,x],\
                        QTY_singlesubsys_totMSIx_siBase[base_id,x]])
            
    print w
            
    print '------------------allMSI_siBase figures for set of all x and total quantity MSIx--------------------'
    w = PrettyTable(["Base_id","BO","FR","MTBF","MWT","MTTR","MAMT","MBaseLDT","MTBPM","MPMT","MDT","Ainh","Aach","Asup","Aop","TotLRU per sub"]) 
    w.max_width = 10
    w.padding_width = 1 # One space between column edges and contents (default)
    w.hrules = ALL          
    
    for base_id in range (1,lru_max_bases+1): 
        for x in range(max_lru):       
            lru_max_bases = len(lrus_base_ids[unique_lru_names[x]]) 
            
            BO_singlesubsys_siBase[base_id] += BO_singlesubsys_totMSIx_siBase[base_id,x]
#             ROS_singlesubsys_siBase[base_id] = LLSAdata[base_id,s[0,x],s[base_id,x],x,P.ROS]
            'Base Corrective Maintenance, Inherent Availability'
            FR_singlesubsys_siBase[base_id] += FR_singlesubsys_totMSIx_siBase[base_id,x]  
            MTBF_singlesubsys_siBase[base_id] = 1/FR_singlesubsys_siBase[base_id]
            UnAinh_singlesubsys_siBase[base_id] += FR_singlesubsys_totMSIx_siBase[base_id,x] * MTTR_singlesubsys_totMSIx_siBase[base_id,x]
            MTTR_singlesubsys_siBase[base_id] = UnAinh_singlesubsys_siBase[base_id]/FR_singlesubsys_siBase[base_id]
            Ainh_singlesubsys_siBase[base_id] = MTBF_singlesubsys_siBase[base_id]/(MTBF_singlesubsys_siBase[base_id]+MTTR_singlesubsys_siBase[base_id])
            
            'Base Preventive Maintenance'
            MTBPM_singlesubsys_siBase[base_id] = MTBPM
            fPM_singlesubsys_siBase[base_id] = 1/MTBPM_singlesubsys_siBase[base_id]
            MPMT_singlesubsys_siBase[base_id] = MPMT 
            
            'Base Achieved Maintenance, Achieved Availability'
            MTBAM_singlesubsys_siBase[base_id] =  1/(FR_singlesubsys_siBase[base_id]+fPM_singlesubsys_siBase[base_id])
            UnAach_singlesubsys_siBase[base_id] = UnAinh_singlesubsys_siBase[base_id] + fPM_singlesubsys_siBase[base_id]* MPMT_singlesubsys_siBase[base_id]
            MAMT_singlesubsys_siBase[base_id] = UnAach_singlesubsys_siBase[base_id]/(FR_singlesubsys_siBase[base_id]+ fPM_singlesubsys_siBase[base_id])
            Aach_singlesubsys_siBase[base_id] = MTBAM_singlesubsys_siBase[base_id]/(MTBAM_singlesubsys_siBase[base_id] + MAMT_singlesubsys_siBase[base_id])
            
            'Base Supply'  
            MWT_singlesubsys_siBase[base_id] = BO_singlesubsys_siBase[base_id]/(FR_singlesubsys_siBase[base_id])
            Asup_singlesubsys_siBase[base_id] = 1-BO_singlesubsys_siBase[base_id]
            MTBDemand_singlesubsys_siBase[base_id] = MWT_singlesubsys_siBase[base_id] /(1/Asup_singlesubsys_siBase[base_id]-1)
#            
            
            'Base Operational Availability'  
            MBaseLDT_singlesubsys_siBase[base_id] = MBaseLDT                
            MTBOM_singlesubsys_siBase[base_id] = MTBAM_singlesubsys_siBase[base_id]
            UnAop_singlesubsys_siBase[base_id] = UnAach_singlesubsys_siBase[base_id]+ FR_singlesubsys_siBase[base_id]*\
            (MWT_singlesubsys_siBase[base_id]+MBaseLDT_singlesubsys_siBase[base_id])
                                            
            MDT_singlesubsys_siBase[base_id]  =  UnAop_singlesubsys_siBase[base_id]/(FR_singlesubsys_siBase[base_id]+ fPM_singlesubsys_siBase[base_id])                                           
            Aop_singlesubsys_siBase[base_id]  = MTBOM_singlesubsys_siBase[base_id]/(MTBOM_singlesubsys_siBase[base_id]+MDT_singlesubsys_siBase[base_id])  
                     
            QTY_singlesubsys_siBase[base_id] += QTY_singlesubsys_totMSIx_siBase[base_id,x] 
            TotStock += s[base_id,x]
            
        w.add_row([base_id,\
                    '%.5f' % BO_singlesubsys_siBase[base_id],\
                    '%.7f' % FR_singlesubsys_siBase[base_id],\
                    '%.1f' % MTBF_singlesubsys_siBase[base_id],\
                    '%.3f' % MWT_singlesubsys_siBase[base_id],\

                    '%.3f' % MTTR_singlesubsys_siBase[base_id],\
                  
                    '%.3f' % MAMT_singlesubsys_siBase[base_id],\
                    '%.3f' % MBaseLDT_singlesubsys_siBase[base_id],\
                    '%.3f' % MTBPM_singlesubsys_siBase[base_id],\
                    '%.3f' % MPMT_singlesubsys_siBase[base_id],\
                    '%.3f' % MDT_singlesubsys_siBase[base_id],\
                    '%.5f' % Ainh_singlesubsys_siBase[base_id],\
                    '%.5f' % Aach_singlesubsys_siBase[base_id],\
                    '%.5f' % Asup_singlesubsys_siBase[base_id],\
                    '%.5f' % Aop_singlesubsys_siBase[base_id],
                     QTY_singlesubsys_siBase[base_id]])
              
                                 
        print'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz'
        print 'allMSI_siBase'
        print 'Base_id: ',base_id
        print 'BO_singlesubsys_siBase','%.7f' % BO_singlesubsys_siBase[base_id]
#         print 'ROS_singlesubsys_siBase','%.7f' % ROS_singlesubsys_siBase[base_id]
        print ' '
        print 'Corrective Maintenance'
        print 'FR_singlesubsys_siBase','%.7f' % FR_singlesubsys_siBase[base_id]
        print 'MTBF_singlesubsys_siBase','%.3f' % MTBF_singlesubsys_siBase[base_id]
        print 'UnAinh_singlesubsys_siBase','%.5f' % UnAinh_singlesubsys_siBase[base_id]
        print 'MTTR_singlesubsys_siBase','%.3f' % MTTR_singlesubsys_siBase[base_id]
        print "Ainh_singlesubsys_siBase", '%.5f' % Ainh_singlesubsys_siBase[base_id]
        print "QTY_subsys",QTY_subsys
        print ''
        print 'Preventive Maintenance'
        print 'MTBPM_singlesubsys_siBase','%.3f' % MTBPM_singlesubsys_siBase[base_id]
        print 'MPMT_singlesubsys_siBase','%.3f' % MPMT_singlesubsys_siBase[base_id]
        print ''
        print 'Achieved Maintenance'
        print 'MTBAM_singlesubsys_siBase','%.3f' % MTBAM_singlesubsys_siBase[base_id]
        print 'MAMT_singlesubsys_siBase','%.3f' % MAMT_singlesubsys_siBase[base_id]
        print 'UnAach_singlesubsys_siBase','%.5f' % UnAach_singlesubsys_siBase[base_id]
        print 'Aach_singlesubsys_siBase','%.5f' % Aach_singlesubsys_siBase[base_id]
        print ' '
        print 'Base Supply'
        print 'MTBDemand_singlesubsys_siBase','%.3f' % MTBDemand_singlesubsys_siBase[base_id]
        print "dient gelijkt te zijn aan: MTBF_singlesubsys_siBase "
        print 'MWT_singlesubsys_siBase','%.3f' % MWT_singlesubsys_siBase[base_id]
        print 'Asup_singlesubsys_siBase','%.5f' % Asup_singlesubsys_siBase[base_id]
        print ''
        print 'Base Operational Availability'  
        print 'MTBOM_singlesubsys_siBase', '%.3f' % MTBOM_singlesubsys_siBase[base_id]
        print 'MDT_singlesubsys_siBase','%.3f' % MDT_singlesubsys_siBase[base_id]
        print 'UnAop_singlesubsys_siBase','%.5f' % UnAop_singlesubsys_siBase[base_id]
        print 'Aop_singlesubsys_siBase','%.5f' % Aop_singlesubsys_siBase[base_id]
        
       
        print''
        print 'TotLRU per subsystem',QTY_singlesubsys_siBase[base_id]
    print w 
         
    BO_allMSI_allBase +=BO_singlesubsys_siBase[base_id]
    FR_allMSI_allBase +=FR_singlesubsys_siBase[base_id]
 
    MWT_allMSI_allBase = BO_allMSI_allBase/FR_allMSI_allBase
    Ainh_allMSI_allBase *= Ainh_singlesubsys_siBase[base_id]
    Asup_allMSI_allBase *=  Asup_singlesubsys_siBase[base_id]
    Aop_allMSI_allBase *=   Aop_singlesubsys_siBase[base_id]
    
#     print w 

    TotStock += s[base_id,x]
    
    print 'ooooooooooooooooooooooooooooooooooooooooooooo '
    print "System parameters"
    print 'BO_allMSI_allBase','%.5f' % BO_allMSI_allBase
    print 'FR_allMSI_allBase',FR_allMSI_allBase
    print 'MTBF_allMSI_allBase',1/FR_allMSI_allBase
    print 'MWT_allMSI_allBase','%.2f' % MWT_allMSI_allBase
     
    print 'Ainh_allMSI_allBase','%.5f' % Ainh_allMSI_allBase
    
    print 'Asup_allMSI_allBase','%.5f' % Asup_allMSI_allBase
#     print '                                   AsupBO_singleWT',1-BO_allMSI_allBase/QTY_subsys
    print 'Aop_allMSI_allBase','%.5f' % Aop_allMSI_allBase 
    
    print''
#     print '                                   Average QTY_singlesubsys_siBase', int(BO_allMSI_allBase/(1-Asup_allMSI_allBase))
    print '++++++++++++++++++++++++++++++++++++++++++QTY_subsys',QTY_subsys
    print 'Tot stock: ',TotStock
#     print 's[max_s[0]][max_s[1]]', s[max_s[0]][max_s[1]], 's[0,x]', s[0,x]
#     print 'count', count,'max_s',max_s,'s[max_s]','%d' % s[max_s[0]][max_s[1]]
#     print 'LLSAdata[base_id,s[0,x],s[max_s[0]][max_s[1]],x,P.EBO]', LLSAdata[base_id,s[0,x],s[max_s[0]][max_s[1]],x,P.EBO]    
#     print '-------------------------------------------------------------' 
#     print ''




    function_detailedresult(LSAdepot,LSAdata,base_id,s,x,lru_max_bases)


    function_stockresult(LSAdepot,LSAdata,base_id,s,x)
    print ''
    print 'function depot count',function_depot_count,',function base count; ',function_base_count
    print 'count: ',count
    print 'BO (bases+depot) system: ', '%.4f' % Tot_EBO[count]
    print ''
    print 'Tottotbases_BO[count]','%.5f' % Tottotbases_BO[count]
    print 'Tottotbases_FR[count]','%.5f' % Tottotbases_FR[count]
    print 'MTBF_SYS', 1/Tottotbases_FR[count]
    print 'MWT','%.2f' % (Tottotbases_BO[count]/Tottotbases_FR[count])
    print 'Availability supply system: ', '%.5f' % Availability_sys[count]
    print ''
    print 'Investment system: ', Investment[count]
    print 'QTY',int(Tottotbases_BO[count]/(1-Availability_sys[count]))
    print 'MPMT', MPMT
    print ''
    
 
    return count

# count = function_ALLocate(LSAdepot,LSAdata)


# dim_count = count+1 #edit thijs
#   
#   
# xaxis = zeros((dim_count)) #edit thijs
# yaxis = zeros((dim_count)) #edit thijs
# Invaxis = zeros((dim_count)) #edit thijs
#   
#   
# for ccount in range(0,dim_count): # edit thijs
#     xaxis[ccount] = ccount
#     yaxis[ccount]= Availability_sys[ccount]
#   
#       
# plt.subplot(311)
# plt.xlabel('round')
# plt.ylabel('Availability')
# plt.title('Availability of each round')
#    
# plt.plot (xaxis,yaxis)
# plt.grid(True)
#  
#   
# for ccount in range(0,dim_count): # edit thijs
#     xaxis[ccount] = ccount
#     Invaxis[ccount]= Investment[ccount]
#    
# plt.subplot(312)
# plt.xlabel('round')
# plt.ylabel('Investment')
# plt.title('Investment of each round')
# plt.plot (xaxis,Invaxis)
# plt.grid(True)
#  
# for ccount in range(0,dim_count): # edit thijs
#     xaxis[ccount] = Investment[ccount]
#     Invaxis[ccount]= Availability_sys[ccount]
# plt.subplot(313)
# plt.xlabel('Investment')
# plt.ylabel('Availability')
# plt.title('Availability function of Investment')
# plt.plot (xaxis,Invaxis)
# plt.grid(True)
# plt.show()
#  
# print("--- %s seconds ---" % (time.time() - start_time))


# function allocation
# BO = ones((max_bases + 1,max_stocklevel_depot,max_stocklevel_bases,max_lru,max_param)) * nan
TotStock = 0
BO = zeros((max_stocklevel_bases,max_bases+1,max_stocklevel_depot)) * nan
max_totalstock = max_stocklevel_depot-1
DeltaBOdepot1= zeros((max_totalstock+1))
BOdepotoptalloc2 = np.array([])
BOdepotoptalloc2_base = np.array([])
contentarray =  zeros((max_totalstock+1))
contentdiffarray = zeros((max_totalstock))
buckcontentarray =  zeros((max_totalstock+1))
buckcontentdiffarray = zeros((max_totalstock))
buckcontentarray_base =  zeros((max_totalstock))

w = PrettyTable(["LRU_id","base_id","s[0]","s[1]","s[2]","s[3]","s[4]","s[5]","BO"]) 
w.max_width = 10
w.padding_width = 1 # One space between column edges and contents (default)
w.hrules = ALL

for x in range(max_lru):
    #initialize for each x
    BOdepotoptalloc2 = np.array([])
    BOdepotoptalloc2_base = np.array([])
    BOdepotoptalloc2_baselevel = np.array([])

    
    lru_max_bases = len(lrus_base_ids[unique_lru_names[x]])
    for s[0,x] in range (max_stocklevel_depot): 
    
        # Bepalen demand rate FR at depot
        LSAdepot[s[0,x],x,P_depot.FR] = 0
        for base_id in range (1,lru_max_bases+1): 
            #L1: Demand rate FR at depot 
            LSAdepot[s[0,x],x,P_depot.FR] += get_lru_data(x, base_id, P_lru.Qty_LRU)*get_lru_data(x,base_id,P_lru.FR)*\
            (1-get_lru_data(x,base_id,P_lru.RHO))             
        
        #Bepaal BO for depot at s[0,x] stocklevel for x
        function_depot(LSAdepot,LSAdata,s,x)         
        # Bepaal BO for all bases and all base_stocklevels for x
        for baselevel  in range (0,max_stocklevel_bases):
            for base_id in range (1,lru_max_bases+1):
                s[base_id,x]=baselevel  
                function_base(LSAdepot,LSAdata,base_id,s,x)      
                BO[s[base_id,x],base_id,s[0,x]] = LSAdata[base_id,s[0,x],s[base_id,x],x,P.EBO]
#                     print 'LRU id: ',x ,'base_id:',base_id,
#                     print 'baselevel','%d' % s[base_id,x],' at s[0,x]','%d' % s[0,x],
#                     print 'BO','%.3f' % LSAdata[base_id,s[0,x],s[base_id,x],x,P.EBO]                          
#                     w.add_row([x,base_id, '%d' % s[0,x],'%d' % s[1,x],'%d' % s[2,x],'%d' % s[3,x],'%d' % s[4,x],'%d' % s[5,x],\
#                                '%.3f' % BO[base_id,s[0,x],s[base_id,x],x]])             
#     print w 
        np.set_printoptions(precision=4)
    
        ' Bepaal Table 3.5 en 3.6 for LRU x'
#     print 'BO.shape',BO.shape
#     print '1 BO' 
#     print BO   


    for s[0,x] in range (max_stocklevel_depot): 
        #Define array BO: Expected backorder level for part #1 at the bases, with depot stocklevel
        #equal to zero (* = BO < 0.0001)
        
#         print '+++++++++++++'
        #select out of BO elements for depotlevel = s[0,x] only
        BOdepot = BO[:,:,s[0,x]]
       
#         print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
#         print 'Table 3.3. BO for depot stocklevel depotlevel:',s[0,x],'and x:',x
#         print BOdepot
        #Create an array of zeros with shape(max_stocklevel_bases,lru_max_bases+1)
        Mask0 = zeros((max_stocklevel_bases,lru_max_bases+1))
        #make first column 1. >>>>>>>>> Baselevel0 >>> contents 1
        Mask0[:,0] = 1
        #Create boolean array . First columm (base_id = 0) = True . Rest False
        BoolMask0_idx = (Mask0 ==1)
#         print 'BoolMask0_idx'
#         print BoolMask0_idx
        # BOdepot  with mask for first column >>>>>>base_id =0
        BOdepot = ma.array(BOdepot,mask = BoolMask0_idx)
#         print ''
#         print 'BOdepot  with mask for first column >>>>>>base_id =0'
#         print BOdepot      
        #Define empty array w
        w = np.array([])
        #Table 3.2 Loop through baselevel (row) and base_id columns) to create DeltaBOdepot. First differences
        for baselevel in range (1,max_stocklevel_bases):
            for base_id in range (0,max_bases+1): 
                DeltaBOdepot = np.subtract(BOdepot[baselevel-1,base_id] ,BOdepot[baselevel,base_id])   
#                 print DeltaBOdepot
                w = np.append(w,DeltaBOdepot)
#                 print w
        DeltaBOdepot= w.reshape(max_stocklevel_bases-1,max_bases+1)
#         print DeltaBOdepot
        #'Insert row 0 all 0'
        DeltaBOdepot = np.insert(DeltaBOdepot, 0, 0, axis=0)
#         print DeltaBOdepot
        
        #DeltaBOdepot array indexing for BOdepot = < 0.0001
        BoolBOdepot_idx = (BOdepot < 0.000000001)      
        DeltaBOdepot = ma.array(DeltaBOdepot,mask = BoolBOdepot_idx)
#         print DeltaBOdepot
        
        #Sum all BOdepot for depotlevel = 0
        BOdepotsumdepot0 = BOdepot[0].sum()
#         print 'BOdepotsumdepot0',BOdepotsumdepot0
        BOdepotoptalloc = np.array([])
        BOdepotoptalloc_base = np.array([])
        # Make first base_id = 0 , just dummy
        BOdepotoptalloc_base =[0]
#         print 'BOdepotoptalloc_base'
        BOdepotoptalloc_baselevel = np.array([])
        # Make first base_id = 0 , just dummy
        BOdepotoptalloc_baselevel =[0]
#         print 'BOdepotoptalloc_base'
        BOdepotoptalloc = np.append(BOdepotoptalloc,BOdepotsumdepot0)
#         print 'BOdepotoptalloc',BOdepotoptalloc
#         print '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'

        for totalstock in range (1,max_stocklevel_bases):   
#             print 'stocklevel depot', s[0,x],'totalstock',max_totalstock
            #Find index  maximum DeltaBOdepot values
            DeltaBOdepot_max = unravel_index(DeltaBOdepot.argmax(),DeltaBOdepot.shape)
#             print ''
#             print 'Table 3.2: DeltaBOdepot at base_id: ',DeltaBOdepot_max[1],':',DeltaBOdepot[DeltaBOdepot_max[0],DeltaBOdepot_max[1]]
           
           
            #Table 3.3. Expected backorders for optimal allocation of stock at the bases for part #1
#             given an depot stock of zero. We use a technique called marginal analysis to find the minimum backorder levels
            #as additional parts are added to the base inventory.
            #create list for single totalstock
            BOdepotoptalloc = np.append(BOdepotoptalloc,BOdepotoptalloc[totalstock-1] - DeltaBOdepot[DeltaBOdepot_max[0],DeltaBOdepot_max[1]])  
#             print 'Table 3.3: BOdepotoptalloc',BOdepotoptalloc
            BOdepotoptalloc_base = np.append(BOdepotoptalloc_base,DeltaBOdepot_max[1]) 
            BOdepotoptalloc_baselevel = np.append(BOdepotoptalloc_baselevel,DeltaBOdepot_max[0])
#             print  'Table 3.3: BOdepotoptalloc_base',BOdepotoptalloc_base

            DeltaBOdepot[DeltaBOdepot_max[0],DeltaBOdepot_max[1]]=ma.masked
#             print 'DeltaBOdepot with masked values'
#             print DeltaBOdepot 
        # create array of for all baselevels   
        BOdepotoptalloc2 = np.append(BOdepotoptalloc2,BOdepotoptalloc)
        BOdepotoptalloc2_base = np.append(BOdepotoptalloc2_base,BOdepotoptalloc_base)
        BOdepotoptalloc2_baselevel = np.append(BOdepotoptalloc2_baselevel,BOdepotoptalloc_baselevel)

    BOdepotoptalloc3= BOdepotoptalloc2.reshape(max_totalstock+1,max_stocklevel_depot)
#     print BOdepotoptalloc3.shape
    BOdepotoptalloc3_base= BOdepotoptalloc2_base.reshape(max_totalstock+1,max_stocklevel_depot)
    BOdepotoptalloc3_baselevel= BOdepotoptalloc2_baselevel.reshape(max_totalstock+1,max_stocklevel_depot)
    print ''
    print 'LRU_id =', x
    print 'Table 3.4: BO and location'
    print ''
    print BOdepotoptalloc3
    print 'base'
    print BOdepotoptalloc3_base
    print 'baselevel'
    print BOdepotoptalloc3_baselevel
    #Table 3.4. Expected backorders for optimal allocation of stock at the RSLs for part
    #1 given an RV stock of zero through eight
    #vertical flip table
    dd= np.fliplr(BOdepotoptalloc3)
    dd_base = np.fliplr(BOdepotoptalloc3_base)
#     print dd
    content = zeros(max_totalstock+1)
    content_base = zeros(max_totalstock+1)
    ccontent = zeros((max_lru))
    #looop through diagonals
    for diag in range (max_totalstock,-1,-1):
        #Select diagonal
        ddd =  dd.diagonal(diag)
        ddd_base = dd_base.diagonal(diag)
#         print 'dd.diagonal(diag)',ddd
#         print 'dd_base.diagonal(diag)',ddd_base
#         print 'diag ' ,diag 
        #Give position 
        dddmin = unravel_index(ddd.argmin(),ddd.shape)
        
        content[diag] = ddd[dddmin]
        content_base[diag] = ddd_base[dddmin]
#         content_base[diag] = ddd_base[dddmin]
#         print 'base',content_base[diag],'content =',content[diag]
#         print  'location min = ',dddmin,'content =',content[diag]
#     print '0000000000000000'
#     print 'content'
#     print content
#     print 'base'
#     print content_base
    
    
    print ''
    print 'Table 3.5 or 3.6'
    #slice the whole sequence, with a step of -1, i.e., in reverse. 
    content = content[::-1]
#     print content.shape

    print 'Minimum BO for x:',x
    print content
    
    contentdiff= -np.diff(content)
    print ''
    print 'Minimum BO first difference array for x:',x
    print contentdiff
    
    
    buckcontent = content/(LRUdepot[x,P_lrudepot.Investment]*get_lru_data(x,base_id,P_lru.Qty_LRU))
    print''
    print 'Bigger bang for buck >>> EBO/Investment'
    print buckcontent
    
    buckcontentdiff= -np.diff(buckcontent)
    print ''
    print 'Minimum BO/Buck first difference array for x:',x
    print buckcontentdiff
    
    content_base = content_base[::-1]
    content_base = content_base[1:]
    print ''
    print 'Minimum BO/Buck first difference base array for x:',x
    print content_base
    
    #Stack for all x
    print ''
    
    contentarray = np.vstack([contentarray, content])
    contentdiffarray = np.vstack([contentdiffarray, contentdiff])
    buckcontentarray = np.vstack([buckcontentarray, buckcontent])
    buckcontentdiffarray = np.vstack([buckcontentdiffarray, buckcontentdiff])
    buckcontentarray_base = np.vstack([buckcontentarray_base, content_base])
    
    if x>= 0:
        diag = max_totalstock+1
        xaxis = zeros((diag))
        yaxis = zeros((diag))
#         points = np.array([])
        ppoints  = zeros((0,2))
        print 'pppppppppppppppppppppppppppppppppp'

#         diag = 0
        for diag in range (max_totalstock,-1,-1):
            xaxis[diag] = diag
            yaxis[diag]= content[diag]
#             print 'xaxis',xaxis[diag],
#             print 'yaxis',yaxis[diag]
            points = np.append(xaxis[diag],yaxis[diag])
           
#             print 'points',points
#             print'diag',max_totalstock+1
#             print 'points.shape',points.shape
            ppoints = np.vstack([ppoints,points])
#             print 'ppoints.shape',ppoints.shape
#         print ppoints

        hull = ConvexHull(ppoints)
        print 'hull.points'
        print hull.points
#         print 'hull.simplices'
#         print hull.simplices
        print 'hull.vertices'
        print hull.vertices
#         print 'hull.neighbors'
#         print hull.neighbors
#         print 'hull.equations'
#         print hull.equations
        y = np.sort(hull.vertices)
        print 'y',y
        #mark function points
        plt.plot(ppoints[:,0], ppoints[:,1], 'o')
        for simplex in hull.simplices:
#             print simplex
#             print ppoints[simplex, 0],ppoints[simplex, 1]
            plt.plot(ppoints[simplex, 0], ppoints[simplex, 1])
        plt.grid(True)
        plt.show()
         
#         plt.xlabel('Total stock level LRU')
#         plt.ylabel('BO')
#         plt.title('BO versus Total stock level LRU')
#         plt.plot (xaxis,yaxis)
#         plt.grid(True)
# #         plt.show()
   
   
print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
print 'Summary'
print 
print 'true maximum stocklevel bases', max_stocklevel_bases-1
print 'true maximum stocklevel depot',max_stocklevel_depot-1
print 'true number of bases',max_bases
print 'true maximum total stock level per LRU ',max_totalstock
print '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

print ''
print 'Table 3.5 and 3.6 Minimum BO/Buck all x'

#delete first row
buckcontentarray = np.delete(buckcontentarray, 0, 0)
print buckcontentarray

print ''
print 'Minimum BO/Buck first difference array for all x'
buckcontentdiffarray = np.delete(buckcontentdiffarray, 0, 0)
print buckcontentdiffarray

contentdiffarray = np.delete(contentdiffarray, 0, 0)


print ''
print 'Minimum BO/Buck first difference base array for al x'
buckcontentarray_base = np.delete(buckcontentarray_base, 0, 0)
print buckcontentarray_base
print ''
#All s[base_id,x] = 0
for x in range(max_lru):
    lru_max_bases = len(lrus_base_ids[unique_lru_names[x]])
    for base_id in range (0,lru_max_bases+1):    
        s[base_id,x] = 0

print'Sum BO all parts for Stocklevel = 0:',contentarray[:,0].sum()

Buck_optimal = contentarray[:,0].sum()
for level in range (max_lru* (max_totalstock)):
    print'Total Stock level',level+1,
    buckcontentdiffarray_max = unravel_index(buckcontentdiffarray.argmax(),buckcontentdiffarray.shape)
    x= buckcontentdiffarray_max[0]
    print 'LRU_id =', x,',',
#     print buckcontentdiffarray
    # in plaats van = ma.masked gebruikt =0 omdat ma.masked een fout oplevert
    Buck_optimal -= contentdiffarray[buckcontentdiffarray_max[0],buckcontentdiffarray_max[1]]
    
    print 'BO:','%.4f' % Buck_optimal,',',
    base_id = int(buckcontentarray_base[buckcontentdiffarray_max[0],buckcontentdiffarray_max[1]])
    print 'base_id = ', base_id,'stocklevel',s[base_id,x]
     
    #remove buckcontentdiffarray element for buckcontentdiffarray = max
    buckcontentdiffarray[buckcontentdiffarray_max[0],buckcontentdiffarray_max[1]]=0
#     print buckcontentdiffarray
   
    if Buck_optimal < 0.1:
        print 'BO Target reached: ' ,Buck_optimal
        print s
        break 
    #increment base_id for buckcontentdiffarray = max
    s[base_id,x] +=1
#     print s
    TotStock += 1

    
    

print '-------------------------------------'   
print s
print 'Tot stock: ',TotStock
print 'QTY_subsys: ',QTY_subsys

A =1- Buck_optimal
print 'Availability',A
print x
