from gillespy2.core import (
    Model,
    Species,
    Reaction,
    Parameter
)
import numpy
import matplotlib.pyplot as plt
class Dimerization(Model):
    def __init__(self, parameter_values=None):
        # First call the gillespy2.Model initializer.
        Model.__init__(self, name='Dimerization')

        # Define parameters for the rates of creation and dissociation.
        k_c = Parameter(name='k_c', expression=0.5)
        k_d = Parameter(name='k_d', expression=0.5)
        self.add_parameter([k_c, k_d])

        # Define variables for the molecular species representing M and D.
        m = Species(name='M', initial_value=30)
        d = Species(name='D',   initial_value=0)
        self.add_species([m, d])

        # The list of reactants and products for a Reaction object are each a
        # Python dictionary in which the dictionary keys are Species objects
        # and the values are stoichiometries of the species in the reaction.
        r_c = Reaction(name="r_creation", rate=k_c, reactants={m:2}, products={d:1})
        r_d = Reaction(name="r_dissociation", rate=k_d, reactants={d:1}, products={m:2})
        self.add_reaction([r_c, r_d])

        # Set the timespan for the simulation.
        self.timespan(numpy.linspace(0, 100, 101))

model = Dimerization()
results = model.run(number_of_trajectories=10)
results.plot()

# Save the figure
plt.savefig("Dimerization.png")  