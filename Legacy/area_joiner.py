def joiner(data, distances, max_diff, max_var = 0):
    import numpy as np
    groups = [[distances[0]]]
    for d in distances[1:]:
        for g in groups:
            if np.abs(np.mean(g) - d) < max_diff:
                g.append(d)
                continue
            else:
                groups.append([d])
    for i,d in enumerate(distances):
        pass