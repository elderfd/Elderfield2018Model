# Elderfield *et al.* 2017 Fungicide Model

## What is this?
The package `fungicide_model` provides an implementation of the models described in Elderfield *et al.* 2017 (in submission). These are ODE models of fungicide resistance evolution in septoria leaf blotch of UK winter wheat and in powdery mildew of French grapevine. 

## How do I use it? (A short example)
You can can use the `Parameters` class to define the parameter values used for the simulation, defaults matching those in the paper are provided for all but the fungicide application strategy (`Parameters.strategy`) and the doses to use (`Parameters.highRiskDose` and `Parameters.lowRiskDose`). For example

```python
from fungicide_model import Parameters, SprayStrategies, Model

# Build a parameters object, specifying which model we intend to use
params = Parameters(Model.WheatSeptoria)

# Choose a fungicide application strategy
params.strategy = SprayStrategies.Mixture

# Choose which doses to apply (NB: should be halved under mixture relative to alternation)
params.highRiskDose = 0.5
params.lowRiskDose = 0.5
```

Once the parameter values have been chosen, they are fed into a `Simulation` object which carries out the numerical simulation itself.

```python
from fungicide_model import Simulation

# Provide the chosen parameters to the simulation
sim = Simulation(params)

# Do the calculation
sim.run()
```

Once the simulation has finished (this may take some time depending on the parameter values chosen) the output can then be accessed directly from the object.

```python
from matplotlib import pyplot as plt

fig, ax = plt.subplots()

# Plot out the amount of healthy tissue over time
ax.plot(
    sim.output.index,
    sim.output["S"],
    color = "#006000"
)

# Print out the summary
print("Selection ratio: {}, effective life: {}, lifetime yield: {}".format(
    sim.selectionRatio,
    sim.seasonsTilCriticalLoss,
    sim.yieldTilCriticalLoss
))

plt.show()
```

## What else can it do?
Investigating the `Parameters` object will show the range of other parameter values that can be tweaked. The most interesting is the ability to turn on and off different features of the model. This is done by setting the boolean state of the member variables: `Parameters.densityDependentInfection`, `Parameters.seasonality`, `Parameters.latentPeriod` and `Parameters.fungicideDecay`. Note that changing these members will automatically change the value of the infection rate parameter (`Parameters.beta`) to match the appropriate model.

## Potential improvements in future
This package produced all of the results for the previously mentioned paper and a good portion of my PhD thesis. There are a number of other similar bits of code I've written for related models or extensions of the models included here that will either be added to this package or released separately as and when the publications attached to them are sent for submission.

Given infinite time the following are things I'd like to add or change,

* Update names of members to match those used in the paper (*e.g.* the `seasonality` member variable refers to Phenology from the paper)
* Bring in some of the extensions I've been working on apart from the attached paper (e.g. partial resistance)
* Attach the test suite to the distribution