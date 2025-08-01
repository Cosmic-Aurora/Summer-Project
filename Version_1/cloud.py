class Cloud:
    # Class Properties
    # l0 - galactic longitude coordinate in the bottom left corner
    # b0 - galactic latitude coordinate in the bottom left corner
    # l1 - galactic longitude coordinate in the top right corner
    # b1 - galactic latitude coordinate in the top right corner
    # dl - galactic longitude change for each pixel moved
    # db - galactic latitude change for each pixel moved
    # lc - central pixel galactic longitude value
    # bc - central pixel galactic latitude value
    # xc - central pixel coordinare index along x-axis
    # yc - central pixel coordinate index along y-axis
    # Note: xc and yc will often fall outside the array size

    # delta_x - pixel count along x-axis
    # delta_y - pixel count along y-axis
    # delta_l - width of file in degrees
    # delta_b - height of the file in degrees
    # data - 2D array representing the cloud and its values
    # distances - 2D array where each row represents a distance and the columns are 
    #   (pixel index in x-direction,
    #   pixel index in y-direction,
    #   galactic longitude, 
    #   galactic latitude, 
    #   distance)
    # mean - mean value of the distances withing the mask
    # median - median value of the distances withing the mask
    # var - variance of the distances within the mask
    # n - number of distance measurements withing the mask
    # diff - difference between furthest distance data point and closest within mask (-1 if there is not data)

    def __init__(self, data, distances, hdr, dx, dy):
        import numpy as np
        self.data = data # Saves and calculates some of the cloud properties based on the header, data as well as dx and dy
        self.lc = hdr["CRVAL1"]
        self.dl = hdr["CD1_1"]
        self.xc = hdr["CRPIX1"] - dx

        self.bc = hdr["CRVAL2"]
        self.db = hdr["CD2_2"]
        self.yc = hdr["CRPIX2"] - dy

        self.delta_y, self.delta_x = np.shape(data)
        self.l0 = self.lc - self.dl*self.xc
        self.b0 = self.bc - self.db*self.yc

        self.l1 = self.l0 + self.delta_x*self.dl
        self.b1 = self.b0 + self.delta_y*self.db

        self.delta_l = np.abs(self.l0-self.l1)
        self.delta_b = np.abs(self.b0-self.b1)

        distances_outside = np.where(np.any((distances[:,0] < min(self.l0,self.l1), distances[:,0] > max(self.l0,self.l1), distances[:,1] < min(self.b0,self.b1), distances[:,1] > max(self.b0,self.b1)), axis = 0))
        pruned_distances = np.delete(distances, distances_outside, axis = 0) # Does an initial removal of all distance data points that are outside of the masked rectangle

        l = pruned_distances[:,0]
        b = pruned_distances[:,1]
        x = (self.xc + (l-self.lc)/self.dl).astype(int) # Transforms the galactic coordinates of the distance points to image pixel coordinates
        y = (self.yc + (b-self.bc)/self.db).astype(int)
        
        outside_cloud = np.where(self.data[y, x] == 0) # Checks which remaining distance points that are outside of the actual cloud area
        x_values = np.array([np.delete(x, outside_cloud)]).T # Removes the outside distances from all arrays and matrices containing them
        y_values = np.array([np.delete(y, outside_cloud)]).T
        final_distances = np.delete(pruned_distances, outside_cloud, axis = 0)
        
        self.distances = np.hstack((x_values,y_values,final_distances)) # Combines the distance matrix (containing the l and b coordinate as well as the distance) with the x and y pixel coordinate
        
        self.n = len(self.distances[:,0]) # Calculates the final properties based on the distances that are actually inside the cloud
        self.mean = np.mean(self.distances[:,4])
        self.median = np.median(self.distances[:,4])
        self.var = np.var(self.distances[:,4])
        self.diff = np.max(self.distances[:,4], initial = -1)-np.min(self.distances[:,4], initial = 99)

    def property_list(self): # Prints out all the cloud properties
        print(f"l0 = {self.l0} - galactic longitude coordinate in the bottom left corner")
        print(f"b0 = {self.b0} - galactic latitude coordinate in the bottom left corner")
        print(f"l1 = {self.l1} - galactic longitude coordinate in the top right corner")
        print(f"b1 = {self.b1} - galactic latitude coordinate in the top right corner")
        print(f"dl = {self.dl} - galactic longitude change for each pixel moved")
        print(f"db = {self.db} - galactic latitude change for each pixel moved")
        print(f"lc = {self.lc} - central pixel galactic longitude value")
        print(f"bc = {self.bc} - central pixel galactic latitude value")
        print(f"xc = {self.xc} - central pixel x-coordinate index on the data array")
        print(f"yc = {self.yc} - central pixel y-coordinate index on the data array")
        print(f"delta_x = {self.delta_x} - pixel count along x-axis")
        print(f"delta_y = {self.delta_y} - pixel count along y-axis")
        print(f"delta_l = {self.delta_l} - width of the file in degrees")
        print(f"delta_b = {self.delta_b} - height of the file in degrees")

    def show_data(self): # Plots the data as well as the distance points within the cloud
        import matplotlib.pyplot as plt
        plt.imshow(self.data, cmap = "YlGn_r", origin = "lower")
        x = self.distances[:,0]
        y = self.distances[:,1]
        plt.scatter(x,y, c = self.distances[:,4], cmap = "cool")
        plt.colorbar()
        plt.show()

    def remove_row(self, i): # Removes a distance measurement from the cloud and then updates all relevant values
        import numpy as np
        self.distances = np.delete(self.distances, i, axis = 0)
        self.n = len(self.distances[:,0])
        self.mean = np.mean(self.distances[:,4])
        self.median = np.median(self.distances[:,4])
        self.var = np.var(self.distances[:,4])
        self.diff = np.max(self.distances[:,4], initial = -1)-np.min(self.distances[:,4], initial = 99)       


def saveclouds(clouds, filename): # Saves a cloud or a set of clouds
    import pickle
    with open(filename, "wb") as file:
        pickle.dump(clouds, file)

def loadclouds(filename): # Opens a cloud or a set of clouds
    import pickle
    with open(filename, "rb") as file:
        clouds = pickle.load(file)
    return clouds