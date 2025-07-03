def area_spread(labeled_image, mask, step_size = 5): # can be used a svoronoi with step_size = big (image size)
    import numpy as np
    import skimage as sk
    if not np.all(labeled_image == 0):
        while True:
            unmasked_image = sk.segmentation.expand_labels(labeled_image, distance = step_size)
            labeled_image = unmasked_image*mask
            print(labeled_image)
            if np.all((labeled_image != 0) == mask):
                return labeled_image
    else:
        raise TypeError("Labeled image must have labels.")

def voronoi(labeled_image, mask):
    return area_spread(labeled_image, mask, np.inf)
