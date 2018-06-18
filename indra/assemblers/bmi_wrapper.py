import copy
import numpy
from pysb import bng
from pysb.simulator import ScipyOdeSimulator

class BMIModel(object):
    def __init__(self, model):
        self.model = model
        self.dt = 1.0
        self.units = 'seconds'
        self.sim = None
        self.attributes = copy.copy(default_attributes)
        self.species_name_map = {}
        # These attributes are related to the simulation state
        self.time = 0.0
        self.state = None

    # Simulation functions
    def initialize(self, fname=None):
        """Initialize the model for simulation, possibly given a config file.

        Parameters
        ----------
        fname : Optional[str]
            The name of the configuration file to load, optional.
        """
        self.sim = ScipyOdeSimulator(self.model)
        self.state = copy.copy(self.sim.initials)[0]
        for idx, species in enumerate(self.model.species):
            monomer = species.monomer_patterns[0].monomer
            self.species_name_map[monomer.name] = idx

        self.time = 0.0

    def update(self, dt=None):
        """Simulate the model for a given time interval.

        Parameters
        ----------
        dt : Optional[float]
            The time step to simulate, if None, the default built-in time step
            is used.
        """
        dt = dt if dt else self.dt
        tspan = [0, dt]
        # Run simulaton with initials set to current state
        res = self.sim.run(tspan=tspan, initials=self.state)
        # Set the state based on the result here
        self.state  = res.species[-1]
        self.time += dt

    def finalize(self):
        """Finish the simulation and clean up resources as needed."""
        pass

    # Setter functions for state variables
    def set_value(self, var_name, value):
        """Set the value of a given variable to a given value.

        Parameters
        ----------
        var_name : str
            The name of the variable in the model whose value should be set.

        value : float
            The value the variable should be set to
        """
        species_idx = self.species_name_map[var_name]
        self.state[species_idx] = value

    # Getter functions for state
    def get_value(self, var_name):
        """Return the value of a given variable.

        Parameters
        ----------
        var_name : str
            The name of the variable whose value should be returned

        Returns
        -------
        value : float
            The value of the given variable in the current state
        """
        species_idx = self.species_name_map[var_name]
        return self.state[species_idx]

    # Getter functions for basic properties
    def get_attribute(self, att_name):
        """Return the value of a given attribute.

        Atrributes include: model_name, version, author_name, grid_type,
        time_step_type, step_method, time_units

        Parameters
        ----------
        att_name : str
            The name of the attribute whose value should be returned.

        Returns
        -------
        value : str
            The value of the attribute
        """
        return self.attributes.get(att_name)

    def get_input_var_names(self):
        """Return a list of variables names that can be set as input.

        Returns
        -------
        var_names : list[str]
            A list of variable names that can be set from the outside
        """
        pass

    def get_output_var_names(self):
        """Return a list of variables names that can be read as output.

        Returns
        -------
        var_names : list[str]
            A list of variable names that can be read from the outside
        """
        pass

    def get_var_name(self, var_name):
        """Return the internal variable name given an outside variable name.

        Parameters
        ----------
        var_name : str
            The name of the outside variable to map

        Returns
        -------
        internal_var_name : str
            The internal name of the corresponding variable
        """
        return var_name

    def get_var_units(self, var_name):
        """Return the units of a given variable.

        Parameters
        ----------
        var_name : str
            The name of the variable whose units should be returned

        Returns
        -------
        unit : str
            The units of the variable
        """
        return '1'

    def get_var_type(self, var_name):
        """Return the type of a given variable.


        Parameters
        ----------
        var_name : str
            The name of the variable whose type should be returned

        Returns
        -------
        unit : str
            The type of the variable as a string
        """
        return 'float64'

    def get_var_rank(self, var_name):
        """Return the matrix rank of the given variable.

        Parameters
        ----------
        var_name : str
            The name of the variable whose rank should be returned

        Returns
        -------
        rank : int
            The dimensionality of the variable, 0 for scalar, 1 for vector,
            etc.
        """
        return numpy.int16(0)

    def get_start_time(self):
        """Return the initial time point of the model.

        Returns
        -------
        start_time : float
            The initial time point of the model.
        """
        return 0.0

    def get_current_time(self):
        """Return the current time point that the model is at during simulation

        Returns
        -------
        time : float
            The current time point
        """
        return self.time

    def get_time_step(self):
        """Return the time step associated with model simulation.

        Returns
        -------
        dt : float
            The time step for model simulation
        """
        return self.dt

    def get_time_units(self):
        """Return the time units of the model simulation.

        Returns
        -------
        units : str
            The time unit of simulation as a string
        """
        return self.units


default_attributes = {
        'model_name': 'INDRA assembled model',
        'version': '1.0',
        'author_name': 'Benjamin M. Gyori',
        'grid_type': 'none',
        'time_step_type': 'fixed',
        'step_method': 'expliit',
        'time_units': 'seconds'
        }
