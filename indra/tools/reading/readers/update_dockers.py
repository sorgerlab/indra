import re
import boto3
import logging
import inspect

from os import path, listdir
from io import BytesIO
from zipfile import ZipFile

from indra.tools.reading.readers import get_reader_classes

HERE = path.dirname(path.abspath(__file__))
DOCKER_TEMPLATE_DIR = path.join(HERE, 'dockerfile_fragments')


logger = logging.getLogger(__name__)


class DockerTemplateError(Exception):
    pass


def _make_dockerfile_rec(template_path):
    """Recursive function to generate a dockerfile from templates."""

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

        other_template_string = _make_dockerfile_rec(other_template_path)
        template = template.replace("{%% %s %%}" % other_template_name,
                                    other_template_string)

    return template


def make_zip_package(rc):
    # Make the dockerfile from the template
    dockerfile, arg_list = get_docker_file(rc)
    if dockerfile is None:
        return None, None

    # Create the buildspec
    build_spec_path = path.join(DOCKER_TEMPLATE_DIR, 'buildspec_fmt.yml')
    with open(build_spec_path, 'r') as f:
        build_spec_fmt = f.read()
    args = ' '.join('--build-arg {arg}=${arg}'.format(arg=arg)
                    for arg in arg_list)
    build_spec = build_spec_fmt.format(args=args)

    # Zip up the buildspec and the dockerfile
    zip_output = BytesIO()
    with ZipFile(zip_output, 'w') as zf:
        zf.writestr('Dockerfile', dockerfile)
        zf.writestr('buildspec.yml', build_spec)
    zip_output.seek(0)

    return zip_output, arg_list


def get_docker_file(rc, logging=True):
    template_path = path.join(path.dirname(inspect.getfile(rc)),
                              'docker_template')
    if not path.exists(template_path):
        if logging:
            logger.info("%s does not have a dockerfile. Continuing." % rc.name)
        return None, None
    if logging:
        logger.info("Forming dockerfile for %s." % rc.name)
    dockerfile = _make_dockerfile_rec(template_path)
    dockerfile += "\n# Set in-{reader} environment variable\n" \
                  "ENV IN_{reader}_DOCKER true\n".format(reader=rc.name)
    arg_list = re.findall('^ARG[ \t]+(.*?)$', dockerfile, re.MULTILINE)
    return dockerfile, arg_list


def get_available_readers():
    reader_spec = {}
    for rc in get_reader_classes():
        _, arg_list = get_docker_file(rc, logging=False)
        if arg_list is None:
            continue
        reader_spec[rc.name] = arg_list
    return reader_spec


def print_help():
    reader_spec = get_available_readers()
    msg = ("usage: update_dockers.py [docker build options]\n"
           "                         [--readers {%s}]"
           % (', '.join(reader_spec.keys())))
    msg += ("\n\nUpdate the specialized docker images on ECR via "
            "CodeBuild.\n\n")
    msg += "optional arguments:\n\n"
    msg += "  -h, --help\t\tShow this message and exit.\n"
    msg += ("  --readers\t\tSpecify which readers' images should be "
            "updated.\n"
            "           \t\tIf not entered, all readers will be "
            "updated.\n")

    msg += "\n\navailable reader dockers:"
    for reader_name, arg_list in reader_spec.items():
        msg += '\n\n  ' + reader_name.lower() + '\t\t'
        msg += '\n\t\t'.join(['--%s' % arg.lower() for arg in arg_list])
    print(msg)
    return


def main():
    from sys import argv

    # Provide some help
    if '--help' in argv or '-h' in argv:
        print_help()
        return

    # Allow the user to limit the readers used.
    only_include_readers = []
    if '--readers' in argv:
        next_idx = argv.index('--readers') + 1
        while next_idx < len(argv) and not argv[next_idx].startswith('-'):
            only_include_readers.append(argv[next_idx].upper())
            next_idx += 1
        if not only_include_readers:
            raise ValueError("At least one reader must be specified with "
                             "--readers.")
        logger.info("Updating: %s" % str(only_include_readers))
    else:
        logger.info("Updating all readers.")

    # Get the AWS clients.
    s3 = boto3.client('s3')
    cb = boto3.client('codebuild')

    for rc in get_reader_classes():
        if only_include_readers and rc.name not in only_include_readers:
            logger.info("%s not included. Skipping." % rc.name)
            continue

        # Put the latest dockerfile etc on s3
        zip_output, arg_list = make_zip_package(rc)
        if zip_output is None:
            continue
        s3_key = ('indra-db/{rdr}-dockerfile/{rdr}-autogen.zip'
                  .format(rdr=rc.name.lower()))
        logger.info("Writing %s to s3." % s3_key)
        s3.put_object(Bucket='bigmech', Key=s3_key, Body=zip_output)

        # Trigger the builds. ENV vars are overwritten based on CLI input.
        # ex: --indra_branch=dev will set the env var INDRA_BRANCH=dev
        env_overrides = []
        for arg in arg_list:
            cli_arg = '--%s' % arg.lower()
            if cli_arg in argv:
                cli_idx = argv.index(cli_arg)
                env_overrides.append({'name': arg, 'value': argv[cli_idx + 1],
                                      'type': 'PLAINTEXT'})

        logger.info("Triggering build for %s with env overrides:\n%s"
                    % (rc.name, env_overrides))
        project_name = 'indra_%s_reading_docker' % rc.name.lower()
        try:
            cb.start_build(projectName=project_name,
                           environmentVariablesOverride=env_overrides)
        except cb.exceptions.ResourceNotFoundException:
            logger.error("Project %s does not exist on AWS. Cannot trigger "
                         "build." % project_name)
    return


if __name__ == '__main__':
    main()
