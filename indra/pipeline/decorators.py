pipeline_functions = {}


def register_pipeline(function):
    if function.__name__ in pipeline_functions:
        raise ExistingFunctionError(
            '%s is already registered with %s.%s' % (
                function.__name__, function.__module__, function.__name__))
    pipeline_functions[function.__name__] = function
    return function


class ExistingFunctionError(Exception):
    pass
