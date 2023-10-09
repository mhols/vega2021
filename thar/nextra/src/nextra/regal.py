import pickle
import pathlib
import os
from nextra import *

DATABASEDIR = pathlib.Path(__file__).parents[0] / '__STORE__'

class Store:

    def __init__(self, storepath=DATABASEDIR):
        
        self.path = pathlib.Path(storepath)
        if not self.path.exists():
            os.mkdir(self.path)
            print(f'new store created at {self.path}')

    def get_path(self):
        return self.path

    def rm(self, key):
        if key in self.keys:
            self.key_to_path(key).unlink()
            print(f"removed {self.key_to_path(key)} from store")
        else:
            print(f"Could not remove: {key} not found in store")

    @property
    def keys(self):
        res = []
        for f in os.listdir(self.path):
            res.append(self.path_to_key(f))
        return res

    def path_to_key(self, p):
        return (str(p.replace('"','/')[:-14]))

    def key_to_path(self, key):
        fname = str(key).replace('/','"') + "_export.pickle"
        return self.path / fname

    def get(self, key):

        try:
            with open(self.key_to_path(key), 'rb') as f:
                o = pickle.load(f)
                return o
        except Exception as ex:
            raise Exception(f'could not open pickel for key {key}, reason: {ex}')

    def store(self, key, ob):
        filepath = self.key_to_path(key)
        try:
            with open(filepath, 'wb') as f:
                pickle.dump(ob, f)
        except Exception as ex:
            raise Exception(f'could not dump {key}, reason {ex}')
