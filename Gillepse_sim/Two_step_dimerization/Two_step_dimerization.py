from gillespy2.core import (
    Model,
    Species,
    Reaction,
    Parameter
)

# variables definitions
k1, k2, k3, k4 = 0.00005, 0.0005, 0, 0
m0, d0, c0 = 50, 50, 0
t, steps = 4000, 500

import numpy as np
import matplotlib.pyplot as plt
class Two_Step_Dimerization(Model):
    def __init__(self, parameter_values=None):
        # First call the gillespy2.Model initializer.
        Model.__init__(self, name='Two_Step_Dimerization')

        # Define parameters for the rates of creation and dissociation.
        k_1 = Parameter(name='k_1', expression=k1)
        k_2 = Parameter(name='k_2', expression=k2)
        k_3 = Parameter(name='k_3', expression=k3)
        k_4 = Parameter(name='k_4', expression=k4)
        self.add_parameter([k_1, k_2, k_3, k_4])

        # Define variables for the molecular species representing M and D.
        m = Species(name='M', initial_value=m0)
        d = Species(name='D',   initial_value=d0)
        c = Species(name='C', initial_value = c0)
        self.add_species([m, d, c])

        # The list of reactants and products for a Reaction object are each a
        # Python dictionary in which the dictionary keys are Species objects
        # and the values are stoichiometries of the species in the reaction.
        r_1 = Reaction(name="r_creation1", rate=k_1, reactants={m:2}, products={d:1})
        r_2 = Reaction(name="r_dissociation2", rate=k_2, reactants={d:1}, products={m:2})
        r_3 = Reaction(name="r_creation3", rate=k_3, reactants={d:2}, products={c:1})
        r_4 = Reaction(name="r_dissociation4", rate=k_4, reactants={c:1}, products={d:2})
        self.add_reaction([r_1, r_2, r_3, r_4])

        # Set the timespan for the simulation.
        self.timespan(np.linspace(0, t, steps))

model = Two_Step_Dimerization()
# the output of run() is a list object [{'time': array([,,,]), 'C':array([,,,])},{...}].
results = model.run(number_of_trajectories=50) 
results.plot(title='Two-step dimerization, (0.005, 0.00005, 0.0005, 0.00005)')
plt.xlabel('Time(s)')
plt.ylabel('Number of molecules')

# Save the figure
plt.savefig("Two-step dimerization, (0.005, 0.00005, 0.0005, 0.00005).png")  

"""
the following part is to calculate the average concentrations of all species, over time
"""

# results[0] is the first entry of list [{'time': array([,,,]), 'C':array([,,,])},{...}].
# it is the dictionary {'time': array([,,,]), 'C':array([,,,])}. 
# results[0].keys() returns 'time', 'M', 'D', 'C', hence species_names.
species_names = results[0].keys()
# the following two are dictionaries
extracted_results = {}
ensemble_averaged_trajectories = {}

for name in species_names:
    # extracted_results, ensemble_averaged_trajectories are labelled by names
    extracted_results[name] = []
    ensemble_averaged_trajectories[name] = []
for i in results:   # i is the i-th entry of results: {'time': array([,,,]), 'C':array([,,,])}
    for name in species_names:
        # i[name] is an array object, when it print out, it will show [   ], but not [,,,], which is a list object.
        # extracted_results is {'time': [array([,,,]), array([,,,]),...], 'C': [array([,,,]), array([,,,]),...] }
        # extracted_results is a dictionary of list of 1-d arrays
        # we can use append function on array or list objects
        extracted_results[name].append(i[name])
	
for name in species_names:
    # converting lists of 1-d arrays into 2-d arrays
	extracted_results[name] = np.array(extracted_results[name])
    # ensemble_average_trajectories is a dictionary of arrays, saving 'time', 'M', etc over time
	ensemble_averaged_trajectories[name] = np.mean(extracted_results[name], axis = 0)
# calculate the equilibrium condition
Keq = []
for i in range(len(ensemble_averaged_trajectories['time'])):
    Keq.append(100*ensemble_averaged_trajectories['D'][i]/ensemble_averaged_trajectories['M'][i]**2)
Keq = np.array(Keq)
    
# plot the figure
fig, ax = plt.subplots()

ax.plot(ensemble_averaged_trajectories['time'], ensemble_averaged_trajectories['M'], label = 'M')
ax.plot(ensemble_averaged_trajectories['time'], ensemble_averaged_trajectories['D'], label = 'D')
#ax.plot(ensemble_averaged_trajectories['time'], ensemble_averaged_trajectories['C'], label = 'C')
ax.plot(ensemble_averaged_trajectories['time'], Keq, label = 'D/M^2')
ax.legend(loc='best')
ax.title.set_text('Average_concentration k1, k2 = 0.00005, 0.0005')
plt.xlabel('Time(s)')
plt.ylabel('Number of molecules')


# Save the figure
fig.savefig('Average_concentration m = 50, d = 50.png')