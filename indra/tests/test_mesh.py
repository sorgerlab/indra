from indra.databases import mesh_client


def test_mesh_id_lookup_from_web():
    mesh_id = 'D003094'
    mesh_name = mesh_client.get_mesh_name_from_web(mesh_id)
    assert mesh_name == 'Collagen'


def test_invalid_id():
    mesh_name = mesh_client.get_mesh_name_from_web('34jkgfh')
    assert mesh_name is None


def test_mesh_id_lookup_local():
    mesh_id = 'D005963'
    mesh_name = mesh_client.get_mesh_name(mesh_id, offline=True)
    assert mesh_name == 'Glucosylceramides'


def test_mesh_id_local_missing():
    mesh_id = 'D015242'
    mesh_name = mesh_client.get_mesh_name(mesh_id, offline=True)
    assert mesh_name is None


def test_mesh_id_fallback_to_rest():
    mesh_id = 'D015242'
    mesh_name = mesh_client.get_mesh_name(mesh_id, offline=False)
    assert mesh_name == 'Ofloxacin'



