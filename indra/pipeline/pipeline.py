import types
import json

from .decorators import pipeline_functions, pipeline
from indra.tools.assemble_corpus import *
from indra.belief.wm_scorer import *
from indra.preassembler.hierarchy_manager import *
from indra.statements import get_statement_by_name


class AssemblyPipeline():
    def __init__(self, steps=None):
        self.steps = steps if steps else []

    def run(self, statements):
        for step in self.steps:
            statements = self.run_function(step, statements)
        return statements

    def append(self, func, *args, **kwargs):
        if isinstance(func_name, types.FunctionType):
            pipeline(func_name)
            func_name = func_name.__name__
        new_step = self.create_new_step(func_name, *args, **kwargs)
        self.steps.append(new_step)

    def insert(self, ix, func_name, *args, **kwargs):
        if isinstance(func_name, types.FunctionType):
            pipeline(func_name)
            func_name = func_name.__name__
        new_step = self.create_new_step(func_name, *args, **kwargs)
        self.steps.insert(ix, new_step)

    def create_new_step(self, func_name, *args, **kwargs):
        self.get_function_from_name(func_name)
        new_step = {'function': func_name}
        if args:
            new_step['args'] = args
        if kwargs:
            new_step['kwargs'] = kwargs
        return new_step

    def get_function_parameters(self, func_dict):
        func_name = func_dict['function']
        args = func_dict.get('args', [])
        kwargs = func_dict.get('kwargs', {})
        return func_name, args, kwargs

    def get_function_from_name(self, name):
        if name in pipeline_functions:
            return pipeline_functions[name]
        raise NotRegisteredFunctionError('%s is not registered' % name)

    def run_simple_function(self, func_name, *args, **kwargs):
        func = self.get_function_from_name(func_name)
        if 'statements' in kwargs:
            statements = kwargs['statements']
            del kwargs['statements']
            return func(statements, *args, **kwargs)
        return func(*args, **kwargs)

    def run_function(self, func_dict, statements=None):
        func_name, args, kwargs = self.get_function_parameters(func_dict)
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
        if not isinstance(argument, dict):
            return False
        if keyword not in argument:
            return False
        return True

    def get_argument_value(self, argument):
        if self.is_function(argument, 'function'):
            # Argument is a function
            if argument.get('no_run', False):
                value = self.get_function_from_name(argument['function'])
            # Argument is a result of a function
            else:
                value = self.run_function(argument)
        # Argument is a statement type
        elif self.is_function(argument, 'stmt_type'):
            value = get_statement_by_name(argument.get('stmt_type'))
        # Argument is a simple value (str, int, boolean, etc.)
        else:
            value = argument
        return value

    @classmethod
    def from_json_file(cls, filename):
        with open(filename, 'r') as f:
            steps = json.load(f)
        ap = AssemblyPipeline(steps)
        return ap


class NotRegisteredFunctionError(Exception):
    pass
