import yaml

class ManageDocUrls:
    
    site_root_url = ""
    command_to_url = {}

    def __init__(self, yaml_path):
        all_yaml = self.read_yaml(yaml_path)
        self.site_root_url = all_yaml["site_url"]
        yaml_nav = all_yaml["nav"]
        for nav in yaml_nav:
            if "Getting Started" in nav.keys():
                for nav_item in nav["Getting Started"]:
                    for nav_item_key in nav_item.keys():
                        if nav_item_key == "Usage":
                            for command in nav_item["Usage"]:
                                for k, v in command.items():
                                    self.command_to_url[k.lower()] = v.replace(".md", "")
                

    def read_yaml(self, yaml_path):
        with open(yaml_path, 'r') as file:
            return yaml.safe_load(file)
    
    def get_url(self, command):
        command = command.lower()
        return self.site_root_url + self.command_to_url[command]

def doc_url(command):
    yml = ManageDocUrls("/home/mabuelanin/dib-dev/dbretina/DBRetina/docs/mkdocs.yml")
    return f"Read more at {yml.get_url(command)}"
