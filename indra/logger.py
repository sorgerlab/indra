import logging
logging.basicConfig(format='%(levelname)s: %(name)s - %(message)s', level=logging.DEBUG)

def get_logger(name):
    return logging.getLogger(name)
