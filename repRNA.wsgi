import sys
sys.path.insert(0, '/var/www/repRNA')

from webserver import app as application
application.debug = True

import os
os.chdir('/var/www/repRNA')