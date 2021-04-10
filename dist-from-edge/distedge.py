from uwimg import *
from flask import Flask
from flask_cors import CORS, cross_origin
import numpy as np

app = Flask(__name__)
cors = CORS(app, resources={r"/foo": {"origins": "*"}})
app.config['CORS_HEADERS'] = 'Content-Type'

@app.route("/get_result/<x>/<y>", methods=['GET', 'POST'])
@cross_origin(origin='*',headers=['Content-Type','Authorization'])
def get_result(x, y):
  im = load_image("data/dog.jpg")
  dist_from_edge = get_smallest_dist_from_edge(im, int(x), int(y))
  min_dist = get_min_dist_edge(im, int(x), int(y))
  save_image(dist_from_edge, "dist-from-edge-test")
  print("python x", x)
  print("python y", y)
  print("min_dist", min_dist)
  return str(min_dist)

if __name__ == "__main__":
    app.run(debug = True)

'''
print("Test 1: point (20, 100)")
im = load_image("data/dog.jpg")
dist_from_edge = get_smallest_dist_from_edge(im, 20, 100)
save_image(dist_from_edge, "dist-from-edge-1")

print("Test 2: point (160, 400)")
im = load_image("data/dog.jpg")
dist_from_edge = get_smallest_dist_from_edge(im, 160, 400)
save_image(dist_from_edge, "dist-from-edge-2")

print("Test 3: point (510, 120)")
im = load_image("data/dog.jpg")
dist_from_edge = get_smallest_dist_from_edge(im, 510, 120)
save_image(dist_from_edge, "dist-from-edge-3")
'''
