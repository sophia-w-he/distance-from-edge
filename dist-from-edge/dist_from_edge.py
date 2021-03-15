from uwimg import *

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
