"""
This script is meant to search in the database 
for different substrates for FeS
Author : Pedram Tavadze
Email  : petavazohi@mix.wvu.edu
"""

import pychemia
from pymatgen.io import vasp


dbsettings={'host'  : 'mongo01.systems.wvu.edu', 
            'name'  : 'PyChemiaMasterDB', 
            'user'  : 'guest', 
            'passwd': 'aldo', 
            'ssl'   : True}
pcdb=pychemia.db.get_database(dbsettings)


