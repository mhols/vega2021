import pickle

DATABASEDIR = './'
STOREDATA = 'stordata.db'
class Store:

    def __init__(self):
        self.items = {}
        self.objects = {}
        try:
            with open(STOREDATA, 'rb') as f:
                self.items = pickle.load(f)
        except:
            print('new store created')

    def get(self, key):

        if key in self.items:
            try:
                return self.objects[key]
            except:
                try:
                    with open (self.items[key], 'rb') as f:
                        ob = pickle.load(f)
                        self.objects[key] = ob
                        return ob 
                except:
                    raise Exception('could not retrieve')

        else:
            raise Exception('no such key')

    def store(self, key, ob):
        filepath =  key +'_extract.pickle'
        try:
            with open(filepath, 'wb') as f:
                pickle.dump(ob, f)
        except Exception as ex:
            raise Exception('could not dump ', ex)
        self.items[key] = filepath
        self.objects[key] = ob

    def save(self):
        for key, o in self.objects.items():
            try:
                with open(self.items[key], 'wb') as f:
                    pickle.dump(o, f)
            except:
                raise Exception('could not save objects')

        try:
            with open(STOREDATA, 'wb') as f:
                pickle.dump(self.items, f)
        except:
            raise Exception('could not save database')
  
    def __setitem__(self, key, ob):
        self.store(key, ob)
    
    def __getitem__(self, key):
        return self.get(key)

    def __delitem__(self, key):
        self.items.pop(key)
        self.objects.pop(key)

    def _hashpath(self, ob):
        return str(hash(str(ob)))+'.pickle'
