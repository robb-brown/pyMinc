try:
	import traceback
	from .pyMinc import *
	from .utilities import *
except:
	traceback.print_exc()
	from pyMinc import *
	from utilities import *
	
#from pyMinctracc import *
