import sys
sys.path.insert(0, '/var/www/Pse-in-One')

from webserver import app as application
application.debug = True

import os
os.chdir('/var/www/Pse-in-One')