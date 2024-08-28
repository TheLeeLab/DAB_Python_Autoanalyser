from gillespy2.core import (
    Model,
    Species,
    Reaction,
    Parameter
)

import numpy as np
import matplotlib.pyplot as plt
class Two_Step_Dimerization(Model):
    def __init__(self, m0, t, steps, fi, n_c, parameter_values=None):
        self.m0 = m0
        self.t = t
        self.steps = steps
        self.fi = fi
        self.n_c = n_c
        f = []
        f_size = [] # all the aggregates of size i-2
        f_size_real = [] # all the aggregates of real size
        # First call the gillespy2.Model initializer.
        Model.__init__(self, name='Two_Step_Dimerization')
        if not parameter_values:
            print('No parameters passed. Terminating.')
            exit()
        else:
            # Define parameters for the rates of creation and dissociation.
            parameter_list = []
            for i in parameter_values:
                parameter_list.append(Parameter(name=i, expression=parameter_values[i]))
            self.add_parameter(parameter_list)

        # Define variables for the molecular species representing M and D.
        m = Species(name='m', initial_value=self.m0)
        for i in range(size):
            f.append(Species(name = 'f'+str(i), initial_value = self.fi[i]))
            f_size.append('f'+str(i))
            f_size_real.append(i+2)
        self.f_size = f_size
        self.f_size_real = f_size_real
        self.add_species([m])
        for i in range(size):
            self.add_species([f[i]])

        # The list of reactants and products for a Reaction object are each a
        # Python dictionary in which the dictionary keys are Species objects
        # and the values are stoichiometries of the species in the reaction.
        # primary nucleation
        rf_m = Reaction(name="rf_m", rate=parameter_list[2], reactants={m:n_c}, products={f[n_c-2]:1})
        # reverse reaction, rate could not be correct
        rr_m = Reaction(name="rr_m", rate=parameter_list[2], reactants={f[n_c-2]:1}, products={m:n_c})
        a = 0
        rf, rr, rrr = [], [], []
        # elongation & fragmentation
        for i in range(size - 1):
            if i < (self.n_c-2): # when i = n_c-3, length = n_c-1, aggregates not form
                rf.append(Reaction(name="rf"+str(i), rate=parameter_list[1], reactants={f[i]:1, m:1}, products={f[i+1]:1}))
                rr.append(Reaction(name="rr"+str(i), rate=parameter_list[1], reactants={f[i+1]:1}, products={f[i]:1, m:1}))
                print("aggregate not form")
            else: 
                # elongation: growth from f[i] to f[i+1]
                rf.append(Reaction(name="rf"+str(i), rate=parameter_list[0], reactants={f[i]:1, m:1}, products={f[i+1]:1}))
                # fragmentation: monomer dissociation, k_off
                rr.append(Reaction(name="rr"+str(i), rate=parameter_list[5], reactants={f[i+1]:1}, products={f[i]:1, m:1}))
                # fragmentation: non-monomer fragmentation, 
                if (i%2) == 0: 
                    for j in range(int(i/2-1)):
                        rrr.append(Reaction(name="rrr"+str(a), rate=parameter_list[4], reactants={f[i]:1}, products={f[j]:1, f[i-j-2]:1}))
                        a = a+1
                    if i>0:
                        rrr.append(Reaction(name="rrr"+str(a), rate=parameter_list[3], reactants={f[i]:1}, products={f[int(i/2-1)]:2}))
                        a = a+1 
                else:
                    for j in range(int((i-1)/2)):
                        rrr.append(Reaction(name="rrr"+str(a), rate=parameter_list[4], reactants={f[i]:1}, products={f[j]:1, f[i-j-2]:1}))
                        a = a+1

        self.add_reaction([rf_m])
        self.add_reaction([rr_m])
        for i in range(size - 1):
            self.add_reaction([rf[i]])
            self.add_reaction([rr[i]])
        for j in range(len(rrr)):
            self.add_reaction([rrr[j]])

        # Set the timespan for the simulation.
        self.timespan(np.linspace(0, self.t, self.steps))
import argparse

full_description = 'This program simulates the fibril size distribution over time.'
parser = argparse.ArgumentParser(description=full_description)
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-kn', '--nucleation_rate', type=float, help='this is the primary nucleation rate', required=True)
requiredNamed.add_argument('-k+', '--elongation_rate', type=float, help='this is the elongation rate', required=True)
requiredNamed.add_argument('-k-', '--fragmentation_rate', type=float, help='this is the elongation grow rate', required=True)
requiredNamed.add_argument('-ko', '--dissociation_rate', type=float, help='this is the dissociation rate', required=True)
requiredNamed.add_argument('-s', '--steps', type=int, help='this is the simulation steps', required=True)
requiredNamed.add_argument('-m', '--monomers', type=int, help='this is the initial monomer number', required=True)
requiredNamed.add_argument('-nc', '--critical_nucleus_size', type=int, help='this is the critical nucleus size', required=True)
args = parser.parse_args()
for arg in vars(args):
    print('Input variables:\n{}: {}'.format(arg,getattr(args,arg)))
kn = args.nucleation_rate
k_plus = args.elongation_rate
k_minus = args.fragmentation_rate
k_off = args.dissociation_rate
steps = args.steps
m0 = args.monomers
n_c = args.critical_nucleus_size
traj = 5

u = (2*(m0*k_plus-k_off)/k_minus)**(0.5)
kappa_inverse = (2*(m0*k_plus-k_off)*k_minus)**(-0.5)
print("average_length =", u)
print("short_time_limit =",kappa_inverse, "sec")
t = kappa_inverse*30
size = int(u*3)
# commend line: python Trial_of_fibril.py -kn 0.00001 -k+ 0.025 -k- 0.015 -t 1000 -s 100 -m 100 -si 49

import time
# variables definitions
k0 = 0
kp, km, km2, koff = k_plus*2, k_minus, k_minus*2, k_off
fi= [] # fi: initial number aggregates, fi[0]: aggregate of size 2;  fi[1]: aggregate of size 3 
for i in range(size):
    fi.append(0)
initial_condition = {'kp':kp,'k0':k0,'kn':kn,'km':km, 'km2':km2, 'koff':koff}

# time calculation starts
time_start = time.time()
model = Two_Step_Dimerization(m0, t, steps, fi, n_c, parameter_values=initial_condition)
# the output of run() is a list object [{'time': array([,,,]), 'f1':array([,,,])},{...}].
results = model.run(number_of_trajectories=traj) 

M_per_trajecotry=[]
time_per_grajectory = np.array(results[0]['time'])

for i in results:
    M_temp = np.zeros_like(i['time'])
    for name in i.keys(): ## i.key(): 'time', 'm', 'f0', 'f1', 'f2', ... 
        if not name=='time'and not name=='m':
            M_temp = M_temp+i[name]*(int(name.replace('f',''))+2)
    M_per_trajecotry.append(M_temp)

# plot the figure
fig, ax = plt.subplots()

for i in range(traj):
    ax.plot(time_per_grajectory, M_per_trajecotry[i], label = 'traj'+str(i+1))
ax.legend(loc='right')
ax.title.set_text('Concentration')
plt.xlabel('Time(s)')
plt.ylabel('Number of molecules')
fig.savefig('l='+str(size)+'t='+str(t)+'m0='+str(m0)+'Per_traj_concentration.png')

"""
the following part is to calculate the average concentrations of all species, over time
"""

# results[0] is the first entry of list [{'time': array([,,,]), 'f1':array([,,,])},{...}].
# it is the dictionary {'time': array([,,,]), 'f1':array([,,,])}. 
# results[0].keys() returns 'time', 'm', 'f0', 'f1', hence species_names.
species_names = results[0].keys()
# the following two are dictionaries
extracted_results = {}
ensemble_averaged_trajectories = {}

for name in species_names:
    # extracted_results, ensemble_averaged_trajectories are labelled by names
    extracted_results[name] = []
    ensemble_averaged_trajectories[name] = []
for i in results:   # i is the i-th entry of results: {'time': array([,,,]), 'f1':array([,,,])}
    for name in species_names:
        # i[name] is an array object, when it print out, it will show [   ], but not [,,,], which is a list object.
        # extracted_results is {'time': [array([,,,]), array([,,,]),...], 'f1': [array([,,,]), array([,,,]),...] }
        # extracted_results is a dictionary of list of 1-d arrays
        # we can use append function on array or list objects
        extracted_results[name].append(i[name])
	
for name in species_names:
    # converting lists of 1-d arrays into 2-d arrays
	extracted_results[name] = np.array(extracted_results[name])
    # ensemble_average_trajectories is a dictionary of arrays, saving 'time', 'm', etc over time
    # ensemble_average_trajectories  = {'time': array([...]), 'f1': array([...]), ...}
	ensemble_averaged_trajectories[name] = np.mean(extracted_results[name], axis = 0)
#print(ensemble_averaged_trajectories)

# time evolution
all_variables = [] # this is a place to store all variables: 'time', 'm', etc.  
# f_size: all the aggregates
for name in model.f_size:
    all_variables.append(ensemble_averaged_trajectories[name])
all_variables = np.array(all_variables)
#print(all_variables)
all_variables_trans = np.transpose(all_variables)
#print(all_variables_trans)
#print(all_variables.shape)
"""
import os
import imageio
filenames = []
for i in range(len(ensemble_averaged_trajectories['time'])):
    # plot the figure
    fig2, ax = plt.subplots()
    ax.plot(model.f_size_real, all_variables_trans[i], label = 't='+str(i*t/steps))
    ax.legend(loc='best')
    ax.title.set_text('Size_distribution')
    plt.xlabel('Aggregate size')
    plt.ylabel('Number of molecules')
    plt.ylim([-0.05, m0/u/10])
    
    # create file name and append it to a list
    filename = f'{i}.png'
    filenames.append(filename)
    
    # save frame
    plt.savefig(filename)
    plt.close()
# build gif
with imageio.get_writer('l='+str(size)+'t='+str(t)+'m0='+str(m0)+'size_distribution_over_time.gif', mode='I') as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)
        
# Remove files
for filename in set(filenames):
    os.remove(filename)

"""
# total aggregate number over time
P = np.sum(all_variables, axis = 0) # total aggregate number over time

# total aggregate mass over time
size_sequence = np.array([np.linspace(2,len(all_variables)+1,len(all_variables))])
#print(size_sequence.shape)
M = np.matmul(size_sequence, all_variables)
#print(M.shape)
#print(M)
M = np.squeeze(M)
#print(M.shape)
#print(all_variables*np.linspace(2,len(all_variables)+1,len(all_variables)))
#M = np.sum(all_variables*np.linspace(2,len(all_variables)+1,len(all_variables)), axis = 0) # linspace(start value, stop value, intervals)
Total_mass = []
for i in range(len(ensemble_averaged_trajectories['m'])):
    Total_mass.append(M[i]+ensemble_averaged_trajectories['m'][i])
    

# calculate the equilibrium condition
#Keq = []
#for i in range(len(ensemble_averaged_trajectories['time'])):
#    Keq.append(10*ensemble_averaged_trajectories['f2'][i]/ensemble_averaged_trajectories['m'][i]**2)
#Keq = np.array(Keq)
    
# plot the figure
fig, ax = plt.subplots()

ax.plot(ensemble_averaged_trajectories['time'], ensemble_averaged_trajectories['m'], label = 'm')
#ax.plot(ensemble_averaged_trajectories['time'], Keq, label = 'f2/m^2')
ax.plot(ensemble_averaged_trajectories['time'], P, label = 'P(t)')
ax.plot(ensemble_averaged_trajectories['time'], M, label = 'M(t)')
ax.plot(ensemble_averaged_trajectories['time'], M/P, label = 'u(t)')
#ax.plot(ensemble_averaged_trajectories['time'], Total_mass, label = 'Total(t)')
#for i in range(int(size/10)):
#    ax.plot(ensemble_averaged_trajectories['time'], ensemble_averaged_trajectories['f'+str(i*10)], label = 'f'+str(i*10+2))
ax.legend(loc='right')
ax.title.set_text('Average_concentration')
plt.xlabel('Time(s)')
plt.ylabel('Number of molecules')
plt.ylim([-0.05, 2*u])


# Save the figure
fig.savefig('l='+str(size)+'t='+str(t)+'m0='+str(m0)+'Average_concentration.png')
time_end = time.time()
print('time elapsed {} sec'.format(time_end - time_start))