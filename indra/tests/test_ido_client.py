from indra.databases import ido_client


def test_lookup():
    """Test IDO name and identifier lookup."""
    name = ido_client.get_ido_name_from_ido_id("0000403")
    assert "parasite role" == name, name

    identifier = ido_client.get_ido_id_from_ido_name("parasite role")
    assert "0000403" == identifier, identifier


if __name__ == '__main__':
    test_lookup()
