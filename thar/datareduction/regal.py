import pickle
import pathlib
import os

DATABASEDIR = pathlib.Path(__file__).parents[0] / '__STORE__'

class Store:

    def __init__(self):
        if not DATABASEDIR.exists():
            os.mkdir(DATABASEDIR)

    @property
    def keys(self):
        res = []
        for f in os.listdir(DATABASEDIR):
            res.append(self.path_to_key(f))
        return res
    
    def path_to_key(self, p):
        return (str(p.replace('"','/')[:-14]))        

    def key_to_path(self, key):
        fname = str(key).replace('/','"') + "_export.pickle"
        return DATABASEDIR / fname

    def get(self, key):

        try:
            with open(self.key_to_path(key), 'rb') as f:
                o = pickle.load(f)
                return o
        except Exception as ex:
            raise Exception('could not open pickel ', ex)
       
    def store(self, key, ob):
        filepath = self.key_to_path(key)
        try:
            with open(filepath, 'wb') as f:
                pickle.dump(ob, f)
        except Exception as ex:
            raise Exception('could not dump ', ex)