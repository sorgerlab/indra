# -*- coding: utf-8 -*-

"""Tests for the MGI client."""

from indra.databases import mgi_client


def test_mgi_client_loaded():
    """Test the MGI client is loaded."""
    assert 'mgi' == mgi_client._client.prefix
    assert mgi_client._client.entries
    assert mgi_client._client.name_to_id


def test_mgi_client_lookup():
    """Test MGI name and identifier lookup."""
    name = mgi_client.get_name_from_id("MGI:97847")
    assert "Raf1" == name, name

    identifier = mgi_client.get_mgi_id_from_mgi_name("Raf1")
    assert "MGI:97847" == identifier, identifier


if __name__ == '__main__':
    test_mgi_client_loaded()
    test_mgi_client_lookup()
