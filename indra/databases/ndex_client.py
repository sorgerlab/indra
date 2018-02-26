from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import io
import json
import time
import requests
import logging
import ndex2.client
from ndex2.niceCXNetwork import NiceCXNetwork

logger = logging.getLogger('ndex_client')

ndex_base_url = 'http://52.37.175.128'

def send_request(ndex_service_url, params, is_json=True, use_get=False):
    """Send a request to the NDEx server.

    Parameters
    ----------
    ndex_service_url : str
        The URL of the service to use for the request.
    params : dict
        A dictionary of parameters to send with the request. Parameter keys
        differ based on the type of request.
    is_json : bool
        True if the response is in json format, otherwise it is assumed to be
        text. Default: False
    use_get : bool
        True if the request needs to use GET instead of POST.

    Returns
    -------
    res : str
        Depending on the type of service and the is_json parameter, this
        function either returns a text string or a json dict.
    """
    if use_get:
        res = requests.get(ndex_service_url, json=params)
    else:
        res = requests.post(ndex_service_url, json=params)
    status = res.status_code
    # If response is immediate, we get 200
    if status == 200:
        if is_json:
            return res.json()
        else:
            return res.text
    # If there is a continuation of the message we get status 300, handled below.
    # Otherwise we return None.
    elif status != 300:
        logger.error('Request returned with code %d' % status)
        return None
    # In case the response is not immediate, a task ID can be used to get
    # the result.
    task_id = res.json().get('task_id')
    logger.info('NDEx task submitted...')
    time_used = 0
    try:
        while status != 200:
            res = requests.get(ndex_base_url + '/task/' + task_id)
            status = res.status_code
            if status != 200:
                time.sleep(5)
                time_used += 5
    except KeyError:
        next
        return None
    logger.info('NDEx task complete.')
    if is_json:
        return res.json()
    else:
        return res.text


def create_network(cx_str, ndex_cred):
    """Creates a new NDEx network of the assembled CX model.

    To upload the assembled CX model to NDEx, you need to have
    a registered account on NDEx (http://ndexbio.org/) and have
    the `ndex` python package installed. The uploaded network
    is private by default.

    Parameters
    ----------
    ndex_cred : dict
        A dictionary with the following entries:
        'user': NDEx user name
        'password': NDEx password

    Returns
    -------
    network_id :  str
        The UUID of the NDEx network that was created by uploading
        the assembled CX model.
    """
    nd = ndex2.client.Ndex2('http://public.ndexbio.org',
                          username=ndex_cred.get('user'),
                          password=ndex_cred.get('password'))
    cx_stream = io.BytesIO(cx_str.encode('utf-8'))
    try:
        logger.info('Uploading network to NDEx.')
        network_uri = nd.save_cx_stream_as_new_network(cx_stream)
    except Exception as e:
        logger.error('Could not upload network to NDEx.')
        logger.error(e)
        return

    network_id = network_uri.rsplit('/')[-1]
    logger.info('The UUID for the uploaded network is: %s' % network_id)
    logger.info('View at: http://ndexbio.org/#/network/%s' % network_id)
    return network_id


def update_network(cx_str, network_id, ndex_cred):
    """Update an existing CX network on NDEx with new CX content.

    Parameters
    ----------
    cx_str : str
        String containing the CX content.
    network_id : str
        UUID of the network on NDEx.
    ndex_cred : dict
        A dictionary with the following entries:
        'user': NDEx user name
        'password': NDEx password
    """
    server = 'http://public.ndexbio.org'
    username = ndex_cred.get('user')
    password = ndex_cred.get('password')
    nd = ndex2.client.Ndex2(server, username, password)

    try:
        logger.info('Getting network summary...')
        summary = nd.get_network_summary(network_id)
    except Exception as e:
        logger.error('Could not get NDEx network summary.')
        logger.error(e)
        return

    # Update network content
    try:
        logger.info('Updating network...')
        cx_stream = io.BytesIO(cx_str.encode('utf-8'))
        nd.update_cx_network(cx_stream, network_id)
    except Exception as e:
        logger.error('Could not update NDEx network.')
        logger.error(e)
        return

    # Update network profile
    ver_str = summary.get('version')
    new_ver = _increment_ndex_ver(ver_str)
    profile = {'name': summary.get('name'),
               'description': summary.get('description'),
               'version': new_ver,
               }
    logger.info('Updating NDEx network (%s) profile to %s',
                network_id, profile)
    profile_retries = 5
    for _ in range(profile_retries):
        try:
            time.sleep(5)
            nd.update_network_profile(network_id, profile)
            break
        except Exception as e:
            logger.error('Could not update NDEx network profile.')
            logger.error(e)

    set_style(network_id, ndex_cred)


def set_style(network_id, ndex_cred):
    # Update network style
    template_uuid = "ea4ea3b7-6903-11e7-961c-0ac135e8bacf"

    server = 'http://public.ndexbio.org'
    username = ndex_cred.get('user')
    password = ndex_cred.get('password')

    source_network = ndex2.create_nice_cx_from_server(username=username,
                                                      password=password,
                                                      uuid=network_id,
                                                      server=server)

    source_network.apply_template(server, template_uuid)

    source_network.update_to(network_id, server=server, username=username,
                             password=password)


def set_provenance(provenance, network_id, ndex_cred):
    server = 'http://public.ndexbio.org'
    username = ndex_cred.get('user')
    password = ndex_cred.get('password')
    nd = ndex2.client.Ndex2(server, username, password)
    try:
        logger.info('Setting network provenance...')
        nd.set_provenance(network_id, provenance)
    except Exception as e:
        logger.error('Could not set network provenance')
        logger.exception(e)


def _increment_ndex_ver(ver_str):
    if not ver_str:
        new_ver = '1.0'
    else:
        major_ver, minor_ver = ver_str.split('.')
        new_minor_ver = str(int(minor_ver) + 1)
        new_ver = major_ver + '.' + new_minor_ver
    return new_ver

