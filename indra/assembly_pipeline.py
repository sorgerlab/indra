import types

from indra.util.decorators import registered_functions, register
from indra.tools.assemble_corpus import *
from indra.belief.wm_scorer import *
from indra.preassembler.hierarchy_manager import *


class AssemblyPipeline():
    def __init__(self, steps=None):
        self.steps = steps if steps else []

    def run(self, statements):
        for step in self.steps:
            self.run_function(step, statements)
        return statements

    def append(self, func, *args, **kwargs):
        if isinstance(func_name, types.FunctionType):
            register(func_name)
            func_name = func_name.__name__
        new_step = self.create_new_step(func_name, *args, **kwargs)
        self.steps.append(new_step)

    def insert(self, ix, func_name, *args, **kwargs):
        if isinstance(func_name, types.FunctionType):
            register(func_name)
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
        if name in registered_functions:
            return registered_functions[name]
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
            if not self.is_function(arg):
                new_args.append(arg)
            elif 'no_run' in arg:
                func = self.get_function_from_name(arg['function'])
                new_args.append(func)
            else:
                result = self.run_function(arg)
                new_args.append(result)
        for k, v in kwargs.items():
            if not self.is_function(v):
                new_kwargs[k] = v
            elif 'no_run' in v:
                func = self.get_function_from_name(v['function'])
                new_kwargs[k] = func
            else:
                result = self.run_function(v)
                new_kwargs[k] = result
        if statements:
            new_kwargs['statements'] = statements
        return self.run_simple_function(func_name, *new_args, **new_kwargs)

    def is_function(self, argument):
        if not isinstance(argument, dict):
            return False
        if 'function' not in argument:
            return False
        return True


class NotRegisteredFunctionError(Exception):
    pass
