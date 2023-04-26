from extract import *
from settings_reference import kwargs

myext_level_1 = Extractor_level_1(kwargs['REFFITSFILE'], **kwargs)
myext_level_1.save_to_store()
del myext_level_1 # no needed any more in memory

# make better lambda map
myext = Extractor(kwargs['REFFITSFILE'], **kwargs)
myext.update()
myext.update()
myext.update()
myext.save_to_store()
