from collections import defaultdict
from datetime import datetime
from os import rename
from pickle import dump, load

import os.path

timestamp = datetime.now().strftime("%Y-%m-%d_%H_%M_%S%z")
log_filename = "logfile_{}.pickle".format(timestamp)
logging = True

def write_log(key, data):
    if not logging:
        return
    # TODO: Keep file open?
    with open(log_filename, "ab") as log_file:
        dump((key, data), log_file)

def load_log(log_file, exclude_keys = set()):
    ret = defaultdict(list)
    with open(log_file, "rb") as log_file:
        while True:
            try:
                key, data = load(log_file)
                if key not in exclude_keys:
                    ret[key].append(data)
            except EOFError:
                break
    return dict(ret)

def change_logfile_name(new_filename):
    global log_filename
    if os.path.exists(log_filename):
        os.rename(log_filename, new_filename)
    log_filename = new_filename

def get_logfile_name():
    return log_filename

def stop_logging():
    global logging
    logging = False
def start_logging():
    global logging
    logging = True
