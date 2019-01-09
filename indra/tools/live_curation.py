from flask import Flask, request, abort, Response
from indra.belief import wm_scorer

scorer = wm_scorer.get_eidos_counts()


app = Flask(__name__)
Compress(app)
CORS(app)


if __name__ == '__main__':
    app.run()
