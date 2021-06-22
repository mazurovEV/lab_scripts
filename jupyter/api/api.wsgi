import logging
import sys
logging.basicConfig(stream=sys.stderr)
sys.path.insert(0, '/home/mazurovev/.local/lib/python3.6/site-packages')
sys.path.insert(0, '/home/mazurovev/scripts/api')
from api import app as application
