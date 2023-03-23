def pytest_configure(config):
    config.addinivalue_line("markers", "webservice: Test using web service")
    config.addinivalue_line("markers", "slow: Test is slow for regular testing")
    config.addinivalue_line("markers", "cron: Test should only run on scheduled test")
    config.addinivalue_line("markers", "nonpublic: Test requires nonpublic environmental variables")
    config.addinivalue_line("markers", "nogha: Test shouldn't be run on GitHub actions")
