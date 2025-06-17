from astropy.io import fits
import numpy as np

hdul_original_masked = fits.open(r"Data/PROMISE/Masked_Data.fits") # Ladda in FITS
data_masked = hdul_original_masked[0].data # Data från original
hdr_masked = hdul_original_masked[0].header

print("Here")

sections_vertical = []
#data_trans = np.transpose(data_masked)

lastx = 0
for x in range(0, np.shape(data_masked)[1]-1):
    if sum(data_masked[:,x]) == 0:
        if lastx + 8000 < x:
            print(f"thing{x}")
            sections_vertical.append(x)
            lastx = x
print("Here2")
sections_vertical.append(120000)
lastx=0
for i, x in enumerate(sections_vertical):
    vertical_slice_masked = data_masked[:,lastx:x]
    hdr_masked = hdul_original_masked[0].header.copy()
    center = hdr_masked["CRPIX1"]
    hdr_masked["CRPIX1"] = center - lastx
    hdu_new_masked = fits.PrimaryHDU(vertical_slice_masked, hdr_masked)
    hdul_new_masked = fits.HDUList([hdu_new_masked])
    hdul_new_masked.writeto(fr"Data/Promise/Slices/Masked_Data_Slice_{i}.fits", overwrite = True) 
    lastx = x

hdul_original_masked.close()
print("Here3")
hdul_original_full = fits.open(r"Data/PROMISE/Full_Data.fits") # Ladda in FITS
data_full = hdul_original_full[0].data # Data från original
print("Here4")
lastx = 0
for i, x in enumerate(sections_vertical):
    vertical_slice_full = data_full[:,lastx:x]
    hdr_full = hdul_original_full[0].header.copy()
    center = hdr_full["CRPIX1"]
    hdr_full["CRPIX1"] = center - lastx
    hdu_new_full = fits.PrimaryHDU(vertical_slice_full, hdr_full)
    hdul_new_full = fits.HDUList([hdu_new_full])
    hdul_new_full.writeto(fr"Data/Promise/Slices/Full_Data_Slice_{i}.fits", overwrite = True)
    lastx = x


hdul_original_full.close()