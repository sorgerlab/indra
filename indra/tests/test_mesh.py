from indra.databases import mesh_client

def test_mesh_id_lookup():
    mesh_id = 'D003094'
    mesh_name = mesh_client.get_mesh_name(mesh_id)
    assert mesh_name == 'Collagen'

def test_invalid_id():
    mesh_name = mesh_client.get_mesh_name('34jkgfh')
    assert mesh_name is None

