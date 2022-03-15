# -*- coding: utf-8 -*-

"""Tests for the Enzyme Class client."""

from indra.databases import ec_client


def test_ec_client_loaded():
    """Test the EC client is loaded."""
    assert 'ec-code' == ec_client._client.prefix
    assert ec_client._client.entries
    assert ec_client._client.name_to_id


def test_ec_client_lookup():
    """Test EC name and identifier lookup."""
    name = ec_client.get_name_from_id("1.1.1.1")
    assert "Alcohol dehydrogenase" == name, name

    identifier = ec_client.get_id_from_name("Alcohol dehydrogenase")
    assert "1.1.1.1" == identifier, identifier


if __name__ == '__main__':
    test_ec_client_loaded()
    test_ec_client_lookup()
