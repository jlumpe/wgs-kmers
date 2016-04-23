from .database import (Database, CURRENT_DB_VERSION, is_db_directory,
	get_db_root, get_current_db, get_default_db, get_db_version)
from .models import *


open_database = Database.open
