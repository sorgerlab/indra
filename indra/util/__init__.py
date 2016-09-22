def has_str(obj):
    if type(obj) == str:
        return True
    # Check for an iterable
    if hasattr(obj, '__iter__'):
        for item in obj:
            item_has_str = has_str(item)
            if item_has_str:
                return True
    if hasattr(obj, '__dict__'):
        for item in obj.__dict__.values():
            item_has_str = has_str(item)
            if item_has_str:
                return True
    return False


