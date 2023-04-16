import os
import json
import urllib.request
import sys

#TODO there are two files with DBRetina_version, remove one.

# Only update this when releasing stable
MAJOR = 2
MINOR = 1
PATCH = 3

PYPI_PACKAGE = "DBRetina"

def is_github_action():
    if "GITHUB_WORKFLOW" in dict(os.environ.items()):
        return True
    else:
        return False    


def get_version():
    
    release_version = f"{MAJOR}.{MINOR}.{PATCH}"
    return release_version