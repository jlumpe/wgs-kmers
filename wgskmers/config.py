import os
from ConfigParser import RawConfigParser

from appdirs import AppDirs


app_dirs = AppDirs('wgskmers', version='0.1')
config_dir = app_dirs.user_config_dir
config_file_path = os.path.join(config_dir, 'config')


def config_exists():
	return os.path.isfile(config_file_path)


def reload_config():
	global config
	config = RawConfigParser()
	config.add_section('databases')
	config.read([config_file_path])

reload_config()


def save_config():
	if not os.path.exists(config_dir):
		os.makedirs(config_dir)
	with open(config_file_path, 'w') as fh:
		config.write(fh)
