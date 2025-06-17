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
    # xc, yc - central pixel coordinate indices on the data array (Note: will often fall outside the array size)

    # length - number of pixels in the x-direction
    # height - number of pixels in the y-direction
    # data - 2D array representing the cloud and its values
    # distances - 2D array where each row represents a distance and the columns are 
    #   (galactic longitude, 
    #   galactic latitude, 
    #   pixel index in x-direction,
    #   pixel index in y-direction,
    #   distance)


    def __init__(self,hdul,position_data):
        import numpy as np
        self.data = hdul[0].data
        hdr = hdul[0].header
        self.lc = hdr["CRVAL1"]
        self.dl = hdr["CD1_1"]
        self.xc = hdr["CRPIX1"]
        self.length = hdr["NAXIS1"]

        self.bc = hdr["CRVAL2"]
        self.db = hdr["CD2_2"]
        self.yc = hdr["CRPIX2"]
        self.height = hdr["NAXIS2"]
        self.l0 = self.lc - self.dl*self.xc
        self.b0 = self.bc - self.db*self.yc

        self.l1 = self.l0 + self.length*self.dl
        self.b1 = self.b0 + self.height*self.db

        remove = np.where(np.any((position_data[:,2] < min(self.l0,self.l1), position_data[:,2] > max(self.l0,self.l1), position_data[:,3] < min(self.b0,self.b1), position_data[:,3] > max(self.b0,self.b1)), axis = 0))
        pruned_data = np.delete(position_data, remove, axis = 0)

        l = pruned_data[:,2]
        b = pruned_data[:,3]

        x = (self.xc + (l-self.lc)/self.dl).astype(int)
        y = (self.yc + (b-self.bc)/self.db).astype(int)

        remove2 = np.where(self.data[y, x] == 0)

        x_values = np.array([np.delete(x, remove2)]).T
        y_values = np.array([np.delete(y, remove2)]).T
        final_data = np.delete(pruned_data, remove2, axis = 0)
        print(x_values,y_values,final_data[:,2:5])
        self.distances = np.hstack((x_values,y_values,final_data[:,2:5]))
        print(self.distances)
    

    def property_list(self):
        print(self.dl)
