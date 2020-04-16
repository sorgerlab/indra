import types
import json
import logging

from .decorators import pipeline_functions, pipeline
from indra.tools.assemble_corpus import *
from indra.belief.wm_scorer import *
from indra.preassembler.hierarchy_manager import *
from indra.statements import get_statement_by_name


logger = logging.getLogger(__name__)


class AssemblyPipeline():
    """An assembly pipeline that runs the specified steps on a given set of
    statements.

    Ways to initialize and run the pipeline (examples assume you have a list
    of INDRA Statements stored in `stmts` variable.)

    >>> from indra.statements import *
    >>> map2k1 = Agent('MAP2K1', db_refs={'HGNC': '6840'})
    >>> mapk1 = Agent('MAPK1', db_refs={'HGNC': '6871'})
    >>> braf = Agent('BRAF')
    >>> stmts = [Phosphorylation(map2k1, mapk1, 'T', '185'),
    ...          Phosphorylation(braf, map2k1)]

    1) Provide a JSON file containing the steps and use a classmethod
    `from_json_file` and run it with `run` method on a list of statements.
    This option allows to store pipeline versions and reproduce the same
    results. All functions referenced in JSON file have to be registered with
    @pipeline decorator.

    >>> ap = AssemblyPipeline.from_json_file('filename.json')
    >>> assembled_stmts = ap.run(stmts)

    2) Initialize a pipeline with a list of steps and run it with `run` method
    on a list of statements.All functions referenced in steps have to be
    registered with @pipeline decorator.

    >>> steps = [
    ...    {"function": "filter_no_hypothesis"},
    ...    {"function": "filter_grounded_only",
    ...     "kwargs": {"score_threshold": 0.8}}
    ... ]
    >>> ap = AssemblyPipeline(steps)
    >>> assembled_stmts = ap.run(stmts)

    3) Initialize an empty pipeline and append/insert the steps one by one.
    Provide a function and its args and kwargs. For arguments that
    require calling a different function, use RunnableArgument class. All
    functions referenced here have to be either imported and passed as function
    objects or registered with @pipeline decorator and passed as function
    names (strings). The pipeline built this way can be optionally saved into
    a JSON file.

    >>> ap = AssemblyPipeline()
    >>> ap.append(filter_no_hypothesis)
    >>> ap.append(filter_grounded_only, score_threshold=0.8)
    >>> ap.append(run_preassembly,
    ...           belief_scorer=RunnableArgument(get_eidos_scorer),
    ...           hierarchies=RunnableArgument(get_wm_hierarchies))
    >>> assembled_stmts = ap.run(stmts)
    >>> ap.to_json_file('filename.json')

    Parameters
    ----------
    steps : list[dict]
        A list of dictionaries representing steps in the pipeline. Each step
        should have a 'function' key and, if appropriate, 'args' and 'kwargs'
        keys. Arguments can be simple values (strings, integers, booleans,
        lists, etc.) or can be functions themselves. In case an argument is a
        function or a result of another function, it should be also
        represented as a dictionary of a similar structure. If a function
        itself is an argument (and not its result), the dictionary should
        contain a key-value pair {'no_run': True}. If an argument is a type
        of a statement, it should be represented as a dictionary {'stmt_type':
        <name of a statement type>}.
    """
    def __init__(self, steps=None):
        self.steps = steps if steps else []

    @classmethod
    def from_json_file(cls, filename):
        """Create an instance of AssemblyPipeline from a JSON file with steps."""
        with open(filename, 'r') as f:
            steps = json.load(f)
        ap = AssemblyPipeline(steps)
        return ap

    def to_json_file(self, filename):
        """Save AssemblyPipeline to a JSON file."""
        with open(filename, 'w') as f:
            json.dump(self.steps, f, indent=1)

    def run(self, statements):
        """Run all steps of the pipeline."""
        logger.info('Running the pipeline')
        for step in self.steps:
            statements = self.run_function(step, statements)
        return statements

    def append(self, func, *args, **kwargs):
        """Append a step to the end of the pipeline.

        Args and kwargs here can be of any type. All functions referenced here
        have to be either imported and passed as function objects or
        registered with @pipeline decorator and passed as function names
        (strings). For arguments that require calling a different function,
        use RunnableArgument class.
        """
        if isinstance(func, types.FunctionType):
            pipeline(func)
            func_name = func.__name__
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
        registered with @pipeline decorator and passed as function names
        (strings). For arguments that require calling a different function,
        use RunnableArgument class.
        """
        if isinstance(func, types.FunctionType):
            pipeline(func)
            func_name = func.__name__
        elif isinstance(func, str):
            func_name = func
        else:
            raise TypeError('Should be a function object or a string')
        new_step = self.create_new_step(func_name, *args, **kwargs)
        self.steps.insert(ix, new_step)

    def create_new_step(self, func_name, *args, **kwargs):
        """Create a dictionary representing a new step in the pipeline."""
        assert self.get_function_from_name(func_name)
        new_step = {'function': func_name}
        if args:
            new_step['args'] = [jsonify_arg_input(arg) for arg in args]
        if kwargs:
            new_step['kwargs'] = {
                k: jsonify_arg_input(v) for (k, v) in kwargs.items()}
        return new_step

    def get_function_parameters(self, func_dict):
        """Retrieve a function name and arguments from function dictionary."""
        func_name = func_dict['function']
        args = func_dict.get('args', [])
        kwargs = func_dict.get('kwargs', {})
        return func_name, args, kwargs

    def get_function_from_name(self, name):
        """Return a function object by name if available or raise exception."""
        if name in pipeline_functions:
            return pipeline_functions[name]
        raise NotRegisteredFunctionError('%s is not registered' % name)

    def run_simple_function(self, func_name, *args, **kwargs):
        """Run a simple function - simple here means a function all arguments
        of which are simple values (do not require extra function calls).
        """
        func = self.get_function_from_name(func_name)
        if 'statements' in kwargs:
            statements = kwargs['statements']
            del kwargs['statements']
            return func(statements, *args, **kwargs)
        return func(*args, **kwargs)

    def run_function(self, func_dict, statements=None):
        """Run a function. For each of the arguments, if it requires an extra
        function call, recursively call the functions until we get a simple
        function.
        """
        func_name, args, kwargs = self.get_function_parameters(func_dict)
        logger.info('%s is called' % func_name)
        new_args = []
        new_kwargs = {}
        for arg in args:
            arg_value = self.get_argument_value(arg)
            new_args.append(arg_value)
        for k, v in kwargs.items():
            kwarg_value = self.get_argument_value(v)
            new_kwargs[k] = kwarg_value
        if statements:
            new_kwargs['statements'] = statements
        return self.run_simple_function(func_name, *new_args, **new_kwargs)

    def is_function(self, argument, keyword='function'):
        """Check if an argument should be converted to a specific object type,
        e.g. a function or a statement type.
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
        if isinstance(func, types.FunctionType):
            pipeline(func)
            self.func_name = func.__name__
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
    if isinstance(arg, types.FunctionType):
        pipeline(arg)
        return {'function': arg.__name__, 'no_run': True}
    if isinstance(arg, str) and arg in pipeline_functions:
        return {'function': arg, 'no_run': True}
    # For some functions Statement type has to be argument
    if isinstance(arg, Statement):
        return {'stmt_type': arg.__name__}
    # Argument is a simple value and can be stored as provided
    return arg
