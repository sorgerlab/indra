import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
DOCKER_TEMPLATE_DIR = 'dockerfile_fragments'


class DockerTemplateError(Exception):
    pass


def make_dockerfile(template_path):
    # Get the list of templates used
    with open(template_path, 'r') as f:
        template = f.read()

    other_templates = re.findall('\{% ([\w_]+) %\}', template)
    for other_template_name in other_templates:
        other_template_path = os.path.join(HERE, DOCKER_TEMPLATE_DIR,
                                           other_template_name)
        if not os.path.exists(other_template_path):
            raise DockerTemplateError("Template %s does not exist."
                                      % other_template_path)

        other_template_string = make_dockerfile(other_template_path)
        template = template.replace("{%% %s %%}" % other_template_name,
                                    other_template_string)

    return template

