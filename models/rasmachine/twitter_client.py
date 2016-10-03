from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import tweepy

def get_oauth(auth_file):
    try:
        fh = open(auth_file, 'rt')
    except IOError:
        print('Could not get Twitter credentials.')
        return None
    lines = [l.strip() for l in fh.readlines()]
    oauth = tweepy.OAuthHandler(lines[0], lines[1])
    oauth.set_access_token(lines[2], lines[3])
    fh.close()
    return oauth

def update_status(msg, auth_file='twitter_cred.txt'):
    twitter_auth = get_oauth(auth_file)
    if twitter_auth is None:
        return
    twitter_api = tweepy.API(twitter_auth)
    twitter_api.update_status(msg)
