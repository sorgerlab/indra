"""This module allows creating a Basic Modeling Interface (BMI) model from
and automatically assembled PySB model. The BMI model can be instantiated
within a simulation workflow system where it is simulated together
with other models."""

import os
import copy
import numpy
import pickle
import textwrap
from lxml import etree
from pysb.bng import generate_equations
from pysb.simulator import ScipyOdeSimulator


class BMIModel(object):
    """This class represents a BMI model wrapping a model assembled by INDRA.

    Parameters
    ----------
    model : pysb.Model
        A PySB model assembled by INDRA to be wrapped in BMI.
    inputs : Optional[list[str]]
        A list of variable names that are considered to be inputs to the model
        meaning that they are read from other models. Note that designating
        a variable as input means that it must be provided by another component
        during the simulation.
    stop_time : int
        The stopping time for this model, controlling the time units up to
        which the model is simulated.
    outside_name_map : dict
        A dictionary mapping outside variables names to inside variable names
        (i.e. ones that are in the wrapped model)
    """
    def __init__(self, model, inputs=None, stop_time=1000,
                 outside_name_map=None):
        self.model = model
        generate_equations(model)

        self.inputs = inputs if inputs else []
        self.stop_time = stop_time
        self.outside_name_map = outside_name_map if outside_name_map else {}

        self.dt = numpy.array(10.0)
        self.units = 'seconds'
        self.sim = None
        self.attributes = copy.copy(default_attributes)
        self.species_name_map = {}
        for idx, species in enumerate(self.model.species):
            monomer = species.monomer_patterns[0].monomer
            self.species_name_map[monomer.name] = idx
        self.input_vars = self._get_input_vars()
        # These attributes are related to the simulation state
        self.state = numpy.array([100.0 for s in self.species_name_map.keys()])
        self.time = numpy.array(0.0)
        self.status = 'start'
        self.time_course = [(self.time, self.state)]
        # EMELI needs a DONE attribute
        self.DONE = False

    def _get_input_vars(self):
        return self.inputs
        # The code below attempts to discover input variables, it is currently
        # inactive but could be made optional later
        # species_is_obj = {s: False for s in self.species_name_map.keys()}
        # for ann in self.model.annotations:
        #     if ann.predicate == 'rule_has_object':
        #         species_is_obj[ann.object] = True
        # # Return all the variables that aren't objects in a rule
        # input_vars = [s for s, tf in species_is_obj.items() if not tf]
        # return input_vars

    # Simulation functions
    def initialize(self, cfg_file=None, mode=None):
        """Initialize the model for simulation, possibly given a config file.

        Parameters
        ----------
        cfg_file : Optional[str]
            The name of the configuration file to load, optional.
        """
        self.sim = ScipyOdeSimulator(self.model)
        self.state = numpy.array(copy.copy(self.sim.initials)[0])
        self.time = numpy.array(0.0)
        self.status = 'initialized'

    def update(self, dt=None):
        """Simulate the model for a given time interval.

        Parameters
        ----------
        dt : Optional[float]
            The time step to simulate, if None, the default built-in time step
            is used.
        """
        # EMELI passes dt = -1 so we need to handle that here
        dt = dt if (dt is not None and dt > 0) else self.dt
        tspan = [0, dt]
        # Run simulaton with initials set to current state
        res = self.sim.run(tspan=tspan, initials=self.state)
        # Set the state based on the result here
        self.state  = res.species[-1]
        self.time += dt
        if self.time > self.stop_time:
            self.DONE = True
        print((self.time, self.state))
        self.time_course.append((self.time.copy(), self.state.copy()))

    def finalize(self):
        """Finish the simulation and clean up resources as needed."""
        self.status = 'finalized'

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
        if var_name in self.outside_name_map:
            var_name = self.outside_name_map[var_name]
            print('%s=%.5f' % (var_name, 1e9*value))
            if var_name == 'Precipitation':
                value = 1e9*value
        species_idx = self.species_name_map[var_name]
        self.state[species_idx] = value

    def set_values(self, var_name, value):
        """Set the value of a given variable to a given value.

        Parameters
        ----------
        var_name : str
            The name of the variable in the model whose value should be set.

        value : float
            The value the variable should be set to
        """
        self.set_value(var_name, value)

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
        if var_name in self.outside_name_map:
            var_name = self.outside_name_map[var_name]
        species_idx = self.species_name_map[var_name]
        return self.state[species_idx]

    def get_values(self, var_name):
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
        return self.get_value(var_name)

    def get_status(self):
        """Return the current status of the model."""
        return self.status

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
        in_vars = copy.copy(self.input_vars)
        for idx, var in enumerate(in_vars):
            if self._map_in_out(var) is not None:
                in_vars[idx] = self._map_in_out(var)
        return in_vars

    def get_output_var_names(self):
        """Return a list of variables names that can be read as output.

        Returns
        -------
        var_names : list[str]
            A list of variable names that can be read from the outside
        """
        # Return all the variables that aren't input variables
        all_vars = list(self.species_name_map.keys())
        output_vars = list(set(all_vars) - set(self.input_vars))
        # Re-map to outside var names if needed
        for idx, var in enumerate(output_vars):
            if self._map_in_out(var) is not None:
                output_vars[idx] = self._map_in_out(var)
        return output_vars

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
        return self._map_out_in(var_name)

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

    def make_repository_component(self):
        """Return an XML string representing this BMI in a workflow.

        This description is required by EMELI to discover and load models.

        Returns
        -------
        xml : str
            String serialized XML representation of the component in the
            model repository.
        """
        component = etree.Element('component')

        comp_name = etree.Element('comp_name')
        comp_name.text = self.model.name
        component.append(comp_name)

        mod_path = etree.Element('module_path')
        mod_path.text = os.getcwd()
        component.append(mod_path)

        mod_name = etree.Element('module_name')
        mod_name.text = self.model.name
        component.append(mod_name)

        class_name = etree.Element('class_name')
        class_name.text = 'model_class'
        component.append(class_name)

        model_name = etree.Element('model_name')
        model_name.text = self.model.name
        component.append(model_name)

        lang = etree.Element('language')
        lang.text = 'python'
        component.append(lang)

        ver = etree.Element('version')
        ver.text = self.get_attribute('version')
        component.append(ver)

        au = etree.Element('author')
        au.text = self.get_attribute('author_name')
        component.append(au)

        hu = etree.Element('help_url')
        hu.text = 'http://github.com/sorgerlab/indra'
        component.append(hu)

        for tag in ('cfg_template', 'time_step_type', 'time_units',
                    'grid_type', 'description', 'comp_type', 'uses_types'):
            elem = etree.Element(tag)
            elem.text = tag
            component.append(elem)

        return etree.tounicode(component, pretty_print=True)

    def export_into_python(self):
        """Write the model into a pickle and create a module that loads it.

        The model basically exports itself as a pickle file and a Python
        file is then written which loads the pickle file. This allows importing
        the model in the simulation workflow.
        """
        pkl_path = self.model.name + '.pkl'
        with open(pkl_path, 'wb') as fh:
            pickle.dump(self, fh, protocol=2)
        py_str = """
        import pickle
        with open('%s', 'rb') as fh:
            model_class = pickle.load(fh)
        """ % os.path.abspath(pkl_path)
        py_str = textwrap.dedent(py_str)
        py_path = self.model.name + '.py'
        with open(py_path, 'w') as fh:
            fh.write(py_str)

    def _map_out_in(self, outside_var_name):
        """Return the internal name of a variable mapped from outside."""
        return self.outside_name_map.get(outside_var_name)

    def _map_in_out(self, inside_var_name):
        """Return the external name of a variable mapped from inside."""
        for out_name, in_name in self.outside_name_map.items():
            if inside_var_name == in_name:
                return out_name
        return None


default_attributes = {
        'model_name': 'indra_model',
        'version': '1.0',
        'author_name': 'Benjamin M. Gyori',
        'grid_type': 'none',
        'time_step_type': 'fixed',
        'step_method': 'explicit',
        'time_units': 'seconds'
        }
