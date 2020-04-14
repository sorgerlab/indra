pipeline_functions = {}


def pipeline(function):
    if function.__name__ not in pipeline_functions:
        pipeline_functions[function.__name__] = function
    return function
