import requests
import json
import time

ndex_base_url = 'http://bel2rdf.bigmech.ndexbio.org'
#ndex_base_url = 'http://52.37.175.128'

def send_request(url_suffix, params):
    res = requests.post(ndex_base_url + url_suffix, data=json.dumps(params))
    res_json = get_result(res)
    return res_json

def get_result(res):
    status = res.status_code
    # If response is immediate, we get 200 
    if status == 200:
        return res.text
    # If there is a continuation of the message
    # we get status 300, handled below. 
    # Otherwise we return None.
    elif status != 300:
        return None
    task_id = res.json()['task_id']
    print 'NDEx task submitted...'
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
    print 'NDEx task complete.'
    return res.text
