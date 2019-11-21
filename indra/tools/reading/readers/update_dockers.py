import re
import boto3
import logging
import inspect

from os import path
from zipfile import ZipFile
from io import BytesIO

from indra.tools.reading.readers import get_reader_classes

HERE = path.dirname(path.abspath(__file__))
DOCKER_TEMPLATE_DIR = path.join(HERE, 'dockerfile_fragments')


logger = logging.getLogger(__name__)


class DockerTemplateError(Exception):
    pass


def make_dockerfile(template_path):
    # Get the list of templates used
    with open(template_path, 'r') as f:
        template = f.read()

    other_templates = re.findall('\{% ([\w_]+) %\}', template)
    for other_template_name in other_templates:
        other_template_path = path.join(DOCKER_TEMPLATE_DIR,
                                        other_template_name)
        if not path.exists(other_template_path):
            raise DockerTemplateError("Template %s does not exist."
                                      % other_template_path)

        other_template_string = make_dockerfile(other_template_path)
        template = template.replace("{%% %s %%}" % other_template_name,
                                    other_template_string)

    return template


def main():
    for rc in get_reader_classes():
        # Make the dockerfile from the template
        template_path = path.join(path.dirname(inspect.getfile(rc)),
                                  'docker_template')
        if not path.exists(template_path):
            logger.info("%s does not have a dockerfile. Continuing." % rc.name)
            continue
        logger.info("Forming dockerfile for %s." % rc.name)
        dockerfile = make_dockerfile(template_path)
        dockerfile += "\n# Set in-{reader} environment variable\n" \
                      "ENV IN_{reader}_DOCKER=true\n".format(reader=rc.name)

        # Create the buildspec
        build_spec_path = path.join(DOCKER_TEMPLATE_DIR, 'buildspec_fmt.yml')
        with open(build_spec_path, 'r') as f:
            build_spec_fmt = f.read()
        arg_list = re.findall('^ARG[ \t]+(.*?)$', dockerfile, re.MULTILINE)
        args = ' '.join('--build-arg {arg}=${arg}'.format(arg=arg)
                        for arg in arg_list)
        build_spec = build_spec_fmt.format(args=args)

        # Zip up the buildspec and the dockerfile
        zip_output = BytesIO()
        with ZipFile(zip_output, 'w') as zf:
            zf.writestr('Dockerfile', dockerfile)
            zf.writestr('buildspec.yml', build_spec)
        zip_output.seek(0)

        # Put the dockerfile on s3
        s3 = boto3.client('s3')
        s3_key = 'indra-db/{rdr}-dockerfile/{rdr}-autogen.zip'.format(rdr=rc.name.lower())
        logger.info("Writing %s to s3." % s3_key)
        s3.put_object(Bucket='bigmech', Key=s3_key, Body=zip_output)

    return


if __name__ == '__main__':
    main()
