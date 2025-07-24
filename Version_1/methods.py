def area_spread(labeled_image, mask, step_size = 5): # Spreads area from a labeled image inside of a mask
    import numpy as np
    import skimage as sk
    if not np.all(labeled_image == 0):
        while True:
            unmasked_image = sk.segmentation.expand_labels(labeled_image, distance = step_size)
            labeled_image = unmasked_image*mask
            if np.all((labeled_image != 0) == mask):
                return labeled_image
    else:
        raise TypeError("Labeled image must have labels.")

def voronoi(labeled_image, mask): # Divides a mask into areas based on a labeled image
    from numpy import inf
    return area_spread(labeled_image, mask, inf)

def points_to_labels(points, data, min_value = 0): # Converts a set of points to a labeled image
    import numpy as np
    labels = np.zeros_like(data).astype(int)
    for i,p in enumerate(points.astype(int)):
        labels[p[1],p[0]] = i + 1 + min_value
    return labels

def dendro_to_distance(d, mask, points, distances): #Initiates distance determination based on a dendrogram and datapoints
    import numpy as np
    distance_matrix = np.zeros_like(mask).astype(float)
    confidence_matrix = np.zeros_like(mask).astype(int)
    for s in d.trunk:
        distance_matrix, confidence_matrix = distance_loop(d, s, mask, distance_matrix, confidence_matrix, points, distances)
    return distance_matrix, confidence_matrix

def distance_loop(d, s, mask, distance_matrix, confidence_matrix, points, distances): # Loops through the dendrogram, assignaing distances and certanties to areas
    import numpy as np
    import sklearn as sk
    import skimage as si
    A = 4 # Non-conflicting clouds
    B = 3 # Non-conflicting leaves/branches with conflicts below
    C = 2 # Conflicting leaves
    D = 1 # Conflicting branches (the rest)
    if s.is_branch:
        for c in s.children: # Loops though all substructures such that the rest of the code starts with leaves and works its way down the structure
            distance_matrix, confidence_matrix = distance_loop(d, c, mask, distance_matrix, confidence_matrix, points, distances)
    ind_x = s.indices()[0] # Extracts the x and y pixel coordinates belonging to the structure
    ind_y = s.indices()[1]
    structure_mask = np.zeros_like(mask)
    structure_mask[ind_x,ind_y] = 1 # Creates a mask of the structure
    points_within = []
    for i, p in enumerate(points): # Finds the indices of the distance measurements that are within the structure
        if structure_mask[int(p[1]),int(p[0])]:
            points_within.append(i)
    points_within = np.array(points_within)
    
    if np.abs(np.max(distances[points_within])-np.min(distances[points_within])) < 0.5: # No conflicting distances within the structure
        distance_matrix[ind_x,ind_y] = float(np.mean(distances[points_within])) # Assign distance mean to the entire structure
        confidence_matrix[ind_x,ind_y] = A # Assign confidence level A to structure
        return distance_matrix, confidence_matrix
    else: 
        if s.is_leaf: # Conflicting distances on leaf
            groups = sk.cluster.DBSCAN(eps = 0.5, min_samples = 1).fit(distances[points_within].reshape(-1,1)) # Groups distances that are close enough
            distance_means = distances[points_within]
            for i in range(np.max(groups.labels_)):    
                distance_means = np.where(groups.labels_ == i, np.mean(distances[np.where(groups.labels_==i)]),distance_means) # Assigns group mean to each group
            labels = points_to_labels(points[points_within], structure_mask)
            labeled_image = area_spread(labels, structure_mask) # From points spread distance labels to full structure 
            for i in range(len(distance_means)):
                distance_matrix = np.where(labeled_image == i+1, distance_means[i], distance_matrix) # Converts labels to distances
            confidence_matrix[ind_x, ind_y] = C # Assign confidence level C
            return distance_matrix, confidence_matrix
        elif s not in [d.structure_at(points[i,::-1].astype(int)) for i in points_within]: # Conflicting distances on branch but no new distances compared to branches/leaves above in the structure
            unique_distances = np.unique(distance_matrix)[1:] # Creates a list of unique distances in the structure that have already been assigned them
            groups = sk.cluster.DBSCAN(eps = 0.5, min_samples = 1).fit(unique_distances.reshape(-1,1)) # Gropus the unique distances
            distance_means = unique_distances
            for i in range(np.max(groups.labels_)): # Applies the mean of each group to each group
                distance_means = np.where(groups.labels_ == i, np.mean(unique_distances[np.where(groups.labels_==i)]),distance_means) 
            labeled_image = np.zeros_like(mask).astype(int)
            for i in range(len(unique_distances)):  # Converts distance image to a labeled image
                labeled_image = np.where(distance_matrix == unique_distances[i], i + 1, labeled_image)
            labeled_image = area_spread(labeled_image, structure_mask) # Spreads the areas so that they now cover the entire structure
            for i in range(len(distance_means)): # Converts back from labeled image to one with distances
                distance_matrix = np.where(labeled_image == i+1, distance_means[i], distance_matrix)
            confidence_matrix = np.where(confidence_matrix*structure_mask == A, B, confidence_matrix) # Replaces confidence level A with B in above structure(s)
            confidence_matrix = np.where(confidence_matrix == 0, D*structure_mask, confidence_matrix) # Assigns confidence level D to current structure
            return distance_matrix, confidence_matrix
        else: # Conflicting distances on branch with new distances added that weren't on leaves/branches higher up in structure
            new_indices = np.where(np.array([d.structure_at(points[i,::-1].astype(int)) for i in points_within]) == s)[0] # Finds indices for new distances
            new_points = points[points_within[new_indices]] # Extracts the pixel coordinates of the new points
            new_distances = distances[points_within[new_indices]] # Extracts the distances of the new points
            unique_distances = np.unique(distance_matrix)[1:] # Finds the unique distances in the old points
            labels = points_to_labels(new_points, structure_mask, len(unique_distances)-len(new_distances)) # Converts the new points to a labeled image
            labeled_image_new = si.segmentation.expand_labels(labels, 25)*structure_mask # Expands each point to a circle with a radius of 25 pixels
            all_distances = np.append(unique_distances, new_distances) # Creates a list with all unique distances

            labeled_image = np.zeros_like(mask).astype(int)
            for i in range(len(unique_distances)): # Converts the old distance image to a labeled image
                labeled_image = np.where(distance_matrix == unique_distances[i], i + 1, labeled_image)
            labeled_image = np.where(labeled_image == 0, labeled_image_new, labeled_image) # Fills in previously undetermined distances with the new distance spread that was done
            labeled_image = area_spread(labeled_image, structure_mask) # Does an area spread of the combined labeled image
            for i in range(len(all_distances)): # Converts back to a distance image
                distance_matrix = np.where(labeled_image == i+1, all_distances[i], distance_matrix)
                
            groups = sk.cluster.DBSCAN(eps = 0.5, min_samples = 1).fit(all_distances.reshape(-1,1)) # Groups the distances
            distance_means = all_distances
            for i in range(np.max(groups.labels_)): # Applies the mean of each group to each group
                distance_means = np.where(groups.labels_ == i, np.mean(all_distances[np.where(groups.labels_==i)]),distance_means)
            for i in range(len(unique_distances)): # Combines the areas that are within a group
                distance_matrix = np.where(distance_matrix == all_distances[i], distance_means[i], distance_matrix)
            confidence_matrix = np.where(confidence_matrix*structure_mask == A, B, confidence_matrix) # Replaces confidence level A with B in above structure(s)
            confidence_matrix = np.where(confidence_matrix == 0, D*structure_mask, confidence_matrix) # Assigns confidence level D to current structure
            return distance_matrix, confidence_matrix