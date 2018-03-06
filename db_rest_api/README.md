# INDRA Database RESTful API

One feature of INDRA is the ability to create and maintain a database of INDRA's pre-assembled statements. Because these statements are a valuable resource in and of themselves, a web api was developed to access that knowledge.

## Implementation

Flask was used to write the python api, as you can see in `api.py`. The means of hosting this api are left to the user. We have had greate success with [Zappa](https://github.com/Miserlou/Zappa) and AWS, and would recommend it for a quick and efficient way to get the API up and running.