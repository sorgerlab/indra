import json
import logging
import inspect

from .decorators import pipeline_functions, register_pipeline
from indra.statements import get_statement_by_name, Statement


logger = logging.getLogger(__name__)


class AssemblyPipeline():
    """An assembly pipeline that runs the specified steps on a given set of
    statements.

    Ways to initialize and run the pipeline (examples assume you have a list
    of INDRA Statements stored in the `stmts` variable.)

    >>> from indra.statements import *
    >>> map2k1 = Agent('MAP2K1', db_refs={'HGNC': '6840'})
    >>> mapk1 = Agent('MAPK1', db_refs={'HGNC': '6871'})
    >>> braf = Agent('BRAF')
    >>> stmts = [Phosphorylation(map2k1, mapk1, 'T', '185'),
    ...          Phosphorylation(braf, map2k1)]

    1) Provide a JSON file containing the steps, then use the classmethod
    `from_json_file`, and run it with the `run` method on a list of statements.
    This option allows storing pipeline versions in a separate file and
    reproducing the same results. All functions referenced in the JSON file
    have to be registered with the @register_pipeline decorator.

    >>> import os
    >>> path_this = os.path.dirname(os.path.abspath(__file__))
    >>> filename = os.path.abspath(
    ... os.path.join(path_this, '..', 'tests', 'pipeline_test.json'))
    >>> ap = AssemblyPipeline.from_json_file(filename)
    >>> assembled_stmts = ap.run(stmts)

    2) Initialize a pipeline with a list of steps and run it with the `run`
    method on a list of statements. All functions referenced in steps have to
    be registered with the @register_pipeline decorator.

    >>> steps = [
    ...    {"function": "filter_no_hypothesis"},
    ...    {"function": "filter_grounded_only",
    ...     "kwargs": {"score_threshold": 0.8}}
    ... ]
    >>> ap = AssemblyPipeline(steps)
    >>> assembled_stmts = ap.run(stmts)

    3) Initialize an empty pipeline and append/insert the steps one by one.
    Provide a function and its args and kwargs. For arguments that
    require calling a different function, use the RunnableArgument class. All
    functions referenced here have to be either imported and passed as function
    objects or registered with the @register_pipeline decorator and passed as
    function names (strings). The pipeline built this way can be optionally
    saved into a JSON file.

    >>> from indra.tools.assemble_corpus import *
    >>> from indra.ontology.world import load_world_ontology
    >>> from indra.belief.wm_scorer import get_eidos_scorer
    >>> ap = AssemblyPipeline()
    >>> ap.append(filter_no_hypothesis)
    >>> ap.append(filter_grounded_only)
    >>> ap.append(run_preassembly,
    ...           belief_scorer=RunnableArgument(get_eidos_scorer),
    ...           ontology=RunnableArgument(load_world_ontology))
    >>> assembled_stmts = ap.run(stmts)
    >>> ap.to_json_file('filename.json')

    Parameters
    ----------
    steps : list[dict]
        A list of dictionaries representing steps in the pipeline. Each step
        should have a 'function' key and, if appropriate, 'args' and 'kwargs'
        keys. Arguments can be simple values (strings, integers, booleans,
        lists, etc.) or can be functions themselves. In case an argument is a
        function or a result of another function, it should also be
        represented as a dictionary of a similar structure. If a function
        itself is an argument (and not its result), the dictionary should
        contain a key-value pair {'no_run': True}. If an argument is a type
        of a statement, it should be represented as a dictionary {'stmt_type':
        <name of a statement type>}.
    """
    def __init__(self, steps=None):
        # This import is here to avoid circular imports
        # It is enough to import one function to get all registered functions
        from indra.tools.assemble_corpus import filter_grounded_only
        from indra.ontology.bio import bio_ontology
        from indra.preassembler.grounding_mapper.gilda import ground_statements
        from indra.preassembler.custom_preassembly import agent_grounding_matches
        self.steps = steps if steps else []

    @classmethod
    def from_json_file(cls, filename):
        """Create an instance of AssemblyPipeline from a JSON file with
        steps."""
        with open(filename, 'r') as f:
            steps = json.load(f)
        ap = AssemblyPipeline(steps)
        return ap

    def to_json_file(self, filename):
        """Save AssemblyPipeline to a JSON file."""
        with open(filename, 'w') as f:
            json.dump(self.steps, f, indent=1)

    def run(self, statements, **kwargs):
        """Run all steps of the pipeline.

        Parameters
        ----------
        statements : list[indra.statements.Statement]
            A list of INDRA Statements to run the pipeline on.
        **kwargs : kwargs
            It is recommended to define all arguments for the steps functions
            in the steps definition, but it is also possible to provide some
            external objects (if it is not possible to provide them as a step
            argument) as kwargs to the entire pipeline here. One should be
            cautious to avoid kwargs name clashes between multiple functions
            (this value will be provided to all functions that expect an
            argument with the same name). To overwrite this value in other
            functions, provide it explicitly in the corresponding steps kwargs.

        Returns
        -------
        list[indra.statements.Statement]
            The list of INDRA Statements resulting from running the pipeline
            on the list of input Statements.
        """
        logger.info('Running the pipeline')
        for step in self.steps:
            statements = self.run_function(step, statements, **kwargs)
        return statements

    def append(self, func, *args, **kwargs):
        """Append a step to the end of the pipeline.

        Args and kwargs here can be of any type. All functions referenced here
        have to be either imported and passed as function objects or
        registered with @register_pipeline decorator and passed as function
        names (strings). For arguments that require calling a different
        function, use RunnableArgument class.

        Parameters
        ----------
        func : str or function
            A function or the string name of a function to add to the pipeline.
        args : args
            Args that are passed to func when calling it.
        kwargs : kwargs
            Kwargs that are passed to func when calling it.
        """
        if inspect.isfunction(func):
            func_name = func.__name__
            if func_name not in pipeline_functions:
                register_pipeline(func)
        elif isinstance(func, str):
            func_name = func
        else:
            raise TypeError('Should be a function object or a string')
        new_step = self.create_new_step(func_name, *args, **kwargs)
        self.steps.append(new_step)

    def insert(self, ix, func, *args, **kwargs):
        """Insert a step to any position in the pipeline.

        Args and kwargs here can be of any type. All functions referenced here
        have to be either imported and passed as function objects or
        registered with @register_pipeline decorator and passed as function
        names (strings). For arguments that require calling a different
        function, use RunnableArgument class.

        Parameters
        ----------
        func : str or function
            A function or the string name of a function to add to the pipeline.
        args : args
            Args that are passed to func when calling it.
        kwargs : kwargs
            Kwargs that are passed to func when calling it.
        """
        if inspect.isfunction(func):
            func_name = func.__name__
            if func_name not in pipeline_functions:
                register_pipeline(func)
        elif isinstance(func, str):
            func_name = func
        else:
            raise TypeError('Should be a function object or a string')
        new_step = self.create_new_step(func_name, *args, **kwargs)
        self.steps.insert(ix, new_step)

    def create_new_step(self, func_name, *args, **kwargs):
        """Create a dictionary representing a new step in the pipeline.

        Parameters
        ----------
        func_name : str
            The string name of a function to create as a step.
        args : args
            Args that are passed to the function when calling it.
        kwargs : kwargs
            Kwargs that are passed to the function when calling it.

        Returns
        -------
        dict
            A dict structure representing a step in the pipeline.
        """
        assert self.get_function_from_name(func_name)
        new_step = {'function': func_name}
        if args:
            new_step['args'] = [jsonify_arg_input(arg) for arg in args]
        if kwargs:
            new_step['kwargs'] = {
                k: jsonify_arg_input(v) for (k, v) in kwargs.items()}
        return new_step

    @staticmethod
    def get_function_parameters(func_dict):
        """Retrieve a function name and arguments from function dictionary.

        Parameters
        ----------
        func_dict : dict
            A dict structure representing a function and its args and kwargs.

        Returns
        -------
        tuple of str, list and dict
            A tuple with the following elements: the name of the function,
            the args of the function, and the kwargs of the function.
        """
        func_name = func_dict['function']
        args = func_dict.get('args', [])
        kwargs = func_dict.get('kwargs', {})
        return func_name, args, kwargs

    @staticmethod
    def get_function_from_name(name):
        """Return a function object by name if available or raise exception.

        Parameters
        ----------
        name : str
            The name of the function.

        Returns
        -------
        function
            The function that was found based on its name. If not found,
            a NotRegisteredFunctionError is raised.
        """
        if name in pipeline_functions:
            return pipeline_functions[name]
        raise NotRegisteredFunctionError('%s is not registered' % name)

    @staticmethod
    def run_simple_function(func, *args, **kwargs):
        """Run a simple function and return the result.

        Simple here means a function all arguments of which are simple values
        (do not require extra function calls).

        Parameters
        ----------
        func : function
            The function to call.
        args : args
            Args that are passed to the function when calling it.
        kwargs : kwargs
            Kwargs that are passed to the function when calling it.

        Returns
        -------
        object
            Any value that the given function returns.
        """
        statements = kwargs.pop('statements', None)
        if statements is not None:
            return func(statements, *args, **kwargs)
        return func(*args, **kwargs)

    def run_function(self, func_dict, statements=None, **kwargs):
        """Run a given function and return the results.

        For each of the arguments, if it requires an extra
        function call, recursively call the functions until we get a simple
        function.

        Parameters
        ----------
        func_dict : dict
            A dict representing the function to call, its args and kwargs.
        args : args
            Args that are passed to the function when calling it.
        kwargs : kwargs
            Kwargs that are passed to the function when calling it.

        Returns
        -------
        object
            Any value that the given function returns.
        """
        func_name, func_args, func_kwargs = self.get_function_parameters(
            func_dict)
        func = self.get_function_from_name(func_name)
        logger.info('Calling %s' % func_name)
        new_args = []
        new_kwargs = {}
        for arg in func_args:
            arg_value = self.get_argument_value(arg)
            new_args.append(arg_value)
        for k, v in func_kwargs.items():
            kwarg_value = self.get_argument_value(v)
            new_kwargs[k] = kwarg_value
        if statements is not None:
            new_kwargs['statements'] = statements
        if kwargs:
            for k, v in kwargs.items():
                if k not in new_kwargs and k in inspect.getargspec(func).args:
                    new_kwargs[k] = v
        return self.run_simple_function(func, *new_args, **new_kwargs)

    @staticmethod
    def is_function(argument, keyword='function'):
        """Check if an argument should be converted to a specific object type,
        e.g. a function or a statement type.

        Parameters
        ----------
        argument : dict or other object
            The argument is a dict, its keyword entry is checked, and if it is
            there, we return True, otherwise we return False.
        keyword : Optional[str]
            The keyword to check if it's there if the argument is a dict.
            Default: function
        """
        if not isinstance(argument, dict):
            return False
        if keyword not in argument:
            return False
        return True

    def get_argument_value(self, arg_json):
        """Get a value of an argument from its json version."""
        if self.is_function(arg_json, 'function'):
            # Argument is a function
            if arg_json.get('no_run', False):
                value = self.get_function_from_name(arg_json['function'])
            # Argument is a result of a function
            else:
                value = self.run_function(arg_json)
        # Argument is a statement type
        elif self.is_function(arg_json, 'stmt_type'):
            value = get_statement_by_name(arg_json.get('stmt_type'))
        # Argument is a simple value (str, int, boolean, etc.)
        else:
            value = arg_json
        return value

    def __len__(self):
        return len(self.steps)

    def __iter__(self):
        return iter(self.steps)


class NotRegisteredFunctionError(Exception):
    pass


class RunnableArgument():
    """Class representing arguments generated by calling a function.

    RunnableArguments should be used as args or kwargs in AssemblyPipeline
    `append` and `insert` methods.

    Parameters
    ----------
    func : str or function
        A function or a name of a function to be called to generate argument
        value.
    """
    def __init__(self, func, *args, **kwargs):
        if inspect.isfunction(func):
            self.func_name = func.__name__
            if self.func_name not in pipeline_functions:
                register_pipeline(func)
        elif isinstance(func, str):
            self.func_name = func
        else:
            raise TypeError('Should be a function object or a string')
        self.args = args
        self.kwargs = kwargs

    def to_json(self):
        """Jsonify to standard AssemblyPipeline step format."""
        json_dict = {'function': self.func_name}
        new_args = []
        new_kwargs = {}
        for arg in self.args:
            new_args.append(jsonify_arg_input(arg))
        for k, v in self.kwargs.items():
            new_kwargs[k] = jsonify_arg_input(v)
        if new_args:
            json_dict['args'] = new_args
        if new_kwargs:
            json_dict['kwargs'] = new_kwargs
        return json_dict


def jsonify_arg_input(arg):
    """Jsonify user input (in AssemblyPipeline `append` and `insert` methods)
    into a standard step json."""
    if isinstance(arg, RunnableArgument):
        return arg.to_json()
    # If a function object or name of a function is provided, we assume it
    # does not have to be run (function itself is argument).
    if inspect.isfunction(arg):
        func_name = arg.__name__
        if func_name not in pipeline_functions:
            register_pipeline(arg)
        return {'function': func_name, 'no_run': True}
    if isinstance(arg, str) and arg in pipeline_functions:
        return {'function': arg, 'no_run': True}
    # For some functions Statement type has to be argument
    if inspect.isclass(arg) and issubclass(arg, Statement):
        return {'stmt_type': arg.__name__}
    # Argument is a simple value and can be stored as provided
    return arg
