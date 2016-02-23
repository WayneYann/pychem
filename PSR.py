# coding: utf-8


"""
A combustor. Two separate stream - one pure methane and the other air, both at
300 K and 60 atm flow into an adiabatic combustor where they mix. We are
interested in the steady-state burning solution. Since at 300 K no reaction
will occur between methane and air, we need to use an 'igniter' to initiate
the chemistry. A simple igniter is a pulsed flow of atomic hydrogen. After the
igniter is turned off, the system approaches the steady burning solution.
"""

import math
import csv

import cantera as ct

# use reaction mechanism GRI-Mech 3.0

gas = ct.Solution('gri30.xml')

# create a reservoir for the fuel inlet, and set to pure methane.
gas.TPX = 300.0, 60*ct.one_atm, 'CH4:1.0'
fuel_in = ct.Reservoir(gas)
fuel_mw = gas.mean_molecular_weight

# use predefined function Air for the air inlet
#air = ct.Solution('air.xml')
#air_in = ct.Reservoir(air)
#air_mw = air.mean_molecular_weight

# use predefined function Air for the air inlet
air = ct.Solution('gri30.xml')
air.TPX = 300.0, 60*ct.one_atm, 'O2:1.0'
air_in = ct.Reservoir(air)
air_mw = air.mean_molecular_weight

# to ignite the fuel/air mixture, we'll introduce a pulse of radicals. The
# steady-state behavior is independent of how we do this, so we'll just use a
# stream of pure atomic hydrogen.
gas.TPX = 300.0, 60*ct.one_atm, 'H:1.0'
igniter = ct.Reservoir(gas)

# create the combustor, and fill it in initially with stoichiometric mixture of CH4 and O2
gas.TPX = 300.0, 60*ct.one_atm, 'CH4:1.0,O2:2'
combustor = ct.IdealGasReactor(gas)
combustor.volume = 1.0

# create a reservoir for the exhaust and initialied with equilibirum composition of the stoichiometric mixture of CH4 and O2
gas.equilibrate("HP")
exhaust = ct.Reservoir(gas)

# lean combustion, phi = 1.0
equiv_ratio = 1.0

# compute fuel and air mass flow rates
factor = 0.1  #0.1
air_mdot = factor * 7.52 * air_mw
fuel_mdot = factor * equiv_ratio * fuel_mw

# create and install the mass flow controllers. Controllers m1 and m2 provide
# constant mass flow rates, and m3 provides a short Gaussian pulse only to
# ignite the mixture
m1 = ct.MassFlowController(fuel_in, combustor, mdot=fuel_mdot)

m2 = ct.MassFlowController(air_in, combustor, mdot=air_mdot)

# The igniter will use a Gaussian time-dependent mass flow rate.
fwhm = 0.2
amplitude = 1.0
t0 = 0.4
igniter_mdot = lambda t: amplitude * math.exp(-(t-t0)**2 * 4 * math.log(2) / fwhm**2)
m3 = ct.MassFlowController(igniter, combustor, mdot=igniter_mdot)

# put a valve on the exhaust line to regulate the pressure
v = ct.Valve(combustor, exhaust, K=1.0)

# the simulation only contains one reactor
sim = ct.ReactorNet([combustor])

# take single steps to 6 s, writing the results to a CSV file for later
# plotting.
tfinal = 6.0
tnow = 0.0
Tprev = combustor.T
tprev = tnow
outfile = open('combustor.csv','wb')
fieldnames = ['time','T','tres']+combustor.thermo.species_names
csvwriter = csv.writer(outfile)
csvwriter.writerow(fieldnames)


while tnow < tfinal:
    tnow = sim.step(tfinal)   #Take a single internal time step toward time *t* [s]. The time after taking the step is returned.
    tres = combustor.mass/v.mdot(tnow)
    Tnow = combustor.T
    if abs(Tnow - Tprev) > 1.0 or tnow-tprev > 2e-2:
        tprev = tnow
        Tprev = Tnow
        csvwriter.writerow([tnow, combustor.T, tres] + list(combustor.thermo.X))
outfile.close()


import numpy as np


selectspecies=["OH","CH4","O2","CO2","H2O","CO"]
selspeX=[]
with open('combustor.csv','rb') as csvfile:
    data = csv.DictReader(csvfile)
    column0 = [row['time'] for row in data]
    csvfile.seek(0)
    columnT = [row['T'] for row in data]
    columnT.pop(0)
    
    for i,item in enumerate(selectspecies):
        csvfile.seek(0)
        column1 = [row[item] for row in data]
        column1.pop(0)
        selspeX.append(column1)
csvfile.close() 


get_ipython().magic(u'matplotlib inline')
from matplotlib import pyplot as plt


for i,item in enumerate(selectspecies):
    plt.plot(column0,selspeX[i],label=item,linewidth=2)
plt.legend(loc=9)
plt.xlim(0.015,1.0)


plt.plot(column0,columnT,label="T",linewidth=2)
plt.legend()
plt.xlim(0.015,1.0)
plt.ylim(0,4500)
