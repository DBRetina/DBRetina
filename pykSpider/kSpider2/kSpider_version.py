import os

#TODO there are two files with DBRetina_version, remove one.

# Only update this when releasing stable
MAJOR = 2
MINOR = 2
PATCH = 8

PYPI_PACKAGE = "DBRetina"

def is_github_action():
    return "GITHUB_WORKFLOW" in dict(os.environ)    


def get_version():
    
    return f"{MAJOR}.{MINOR}.{PATCH}"