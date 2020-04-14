registered_functions = {}


def register(function):
    if function.__name__ not in registered_functions:
        registered_functions[function.__name__] = function
    return function
