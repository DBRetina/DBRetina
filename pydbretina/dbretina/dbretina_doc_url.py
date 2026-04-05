def doc_url(command):
    _url =  f"https://dbretina.github.io/DBRetina/usage/dbretina_{command}"
    output =  f"Read more at {_url}"
    output = f"\033[1m{output}\033[0m"

    return output
