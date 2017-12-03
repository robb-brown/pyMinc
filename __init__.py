try:
	from .pyMinc import *
	from .utilities import *
	try:
		from .pyMinctracc import *
	except:
		pass
except:
	from pyMinc import *
	from utilities import *
	try:
		from pyMinctracc import *
	except:
		pass
	