import os
import json
import urllib.request
import sys

#TODO there are two files with DBRetina_version, remove one.

# Only update this when releasing stable
MAJOR = 2
MINOR = 2
PATCH = 3

PYPI_PACKAGE = "DBRetina"

def is_github_action():
    return "GITHUB_WORKFLOW" in dict(os.environ)    


def get_version():
    
    return f"{MAJOR}.{MINOR}.{PATCH}"