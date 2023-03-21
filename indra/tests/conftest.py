def pytest_configure(config):
    config.addinivalue_line("markers", "webservice: Test using web service")
