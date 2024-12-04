from PIL import Image
import glob
from PIL import EpsImagePlugin

EpsImagePlugin.gs_windows_binary = r'C:\Program Files\Ghostscript'
folder = 'C:/Users/Werk/surfdrive/Proefschrift/CO2_on_CO2_physisorption/Frontiers/Images/'

ext_in = ".eps"
ext_out = ".jpg"

paths = glob.glob(folder+"*"+ext_in)

for image_path in paths:
    eps_image = Image.open(image_path)
    eps_image.load(scale=10)
    out_path = image_path[:-len(ext_in)]+ext_out
    print(out_path)
    # eps_image.save(out_path)