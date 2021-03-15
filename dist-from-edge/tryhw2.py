from uwimg import *
im = load_image("data/dog.jpg")
f = make_box_filter(7)
blur = convolve_image(im, f, 1)
save_image(blur, "dog-box7")

im = load_image("data/dog.jpg")
f = make_box_filter(7)
blur = convolve_image(im, f, 1)
thumb = nn_resize(blur, blur.w//7, blur.h//7)
save_image(thumb, "dogthumb")

im = load_image("data/dog.jpg")
f = make_gaussian_filter(2)
blur = convolve_image(im, f, 1)
save_image(blur, "dog-gauss2")

im = load_image("data/dog.jpg")
f = make_gaussian_filter(2)
lfreq = convolve_image(im, f, 1)
hfreq = im - lfreq
reconstruct = lfreq + hfreq
save_image(lfreq, "low-frequency")
save_image(hfreq, "high-frequency")
save_image(reconstruct, "reconstruct")

im = load_image("data/dog.jpg")
res = sobel_image(im)
mag = res[0]
# normalize_image(mag)
# did not work so had to directly apply in c file filter_image.c
save_image(mag, "magnitude")

im = load_image("data/dog.jpg")
dist_from_edge = get_smallest_dist_from_edge(im, 20, 100)
save_image(dist_from_edge, "dist-from-edge")

# EXTRA CREDIT
#im = load_image("figs/salt_petter_building.jpg")
#med = apply_median_filter(im, 3)
#save_image(med, "building-median")
