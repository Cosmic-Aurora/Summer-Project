# ID, l, b (of highest extinction value within the cloud), Area (physical area), distance, flag, mass (sum over pixel magnitudes*area of each pixel (use distance)) (assume square pixels, small enough angels that b and l dont matter), # of points

# Comparing topdown to other images (sent by Jouni)
import cloud

clouds = cloud.loadclouds("Data/clouds.pkl")

cl = clouds["cloud_1"]
cl.property_list()
cl.show_data()
