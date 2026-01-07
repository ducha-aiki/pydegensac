import importlib

def test_import_pydegensac():
    module = importlib.import_module("pydegensac")
    assert module is not None
