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
    #   mean - mean value of the distances withing the mask
    #   var - variance of the distances within the mask
    #   n - number of distance measurements withing the mask
    #   diff - difference between furthest distance data point and closest within mask (-1 if there is not data)


    def __init__(self,hdul,position_data):
        import numpy as np
        self.data = hdul[0].data
        hdr = hdul[0].header
        self.lc = hdr["CRVAL1"]
        self.dl = hdr["CD1_1"]
        self.xc = hdr["CRPIX1"]
        self.delta_x = hdr["NAXIS1"]

        self.bc = hdr["CRVAL2"]
        self.db = hdr["CD2_2"]
        self.yc = hdr["CRPIX2"]
        self.delta_y = hdr["NAXIS2"]
        self.l0 = self.lc - self.dl*self.xc
        self.b0 = self.bc - self.db*self.yc

        self.l1 = self.l0 + self.delta_x*self.dl
        self.b1 = self.b0 + self.delta_y*self.db

        self.delta_l = np.abs(self.l0-self.l1)
        self.delta_b = np.abs(self.b0-self.b1)

        remove = np.where(np.any((position_data[:,0] < min(self.l0,self.l1), position_data[:,0] > max(self.l0,self.l1), position_data[:,1] < min(self.b0,self.b1), position_data[:,1] > max(self.b0,self.b1)), axis = 0))
        pruned_data = np.delete(position_data, remove, axis = 0)

        l = pruned_data[:,0]
        b = pruned_data[:,1]

        x = (self.xc + (l-self.lc)/self.dl).astype(int)
        y = (self.yc + (b-self.bc)/self.db).astype(int)
        remove2 = np.where(self.data[y, x] == 0)

        x_values = np.array([np.delete(x, remove2)]).T
        y_values = np.array([np.delete(y, remove2)]).T
        final_data = np.delete(pruned_data, remove2, axis = 0)
        self.distances = np.hstack((x_values,y_values,final_data))
        self.n = len(self.distances[:,0])
        self.mean = np.mean(self.distances[:,4])
        self.var = np.var(self.distances[:,4])
        self.diff = np.max(self.distances[:,4], initial = -1)-np.min(self.distances[:,4], initial = -1)
    

    def property_list(self):
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

    def show_data(self, plt):
        plt.imshow(self.data)
        x = self.distances[:,0]
        y = self.distances[:,1]
        plt.scatter(x,y, color = "green")
        plt.show()


def saveclouds(clouds, filename):
    import pickle
    with open(filename, "wb") as file:
        pickle.dump(clouds, file)

def loadclouds(filename):
    import pickle
    with open(filename, "rb") as file:
        clouds = pickle.load(file)
    return clouds