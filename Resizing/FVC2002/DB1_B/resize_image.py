
#!/usr/bin/env python

from PIL import Image
from resizeimage import resizeimage
import os, sys

def resizeImage(infile, output_dir="", size=(1024,768)):
     outfile = os.path.splitext(infile)[0]
     extension = os.path.splitext(infile)[1]

     if infile != outfile:
        try :
            im = Image.open(infile)
            cover = resizeimage.resize_cover(im, [256, 256])
            cover.save(output_dir+outfile+extension, im.format)
        except IOError:
            print "cannot reduce image for ", infile


if __name__=="__main__":
    output_dir = "resized"
    dir = os.getcwd()

    if not os.path.exists(os.path.join(dir,output_dir)):
        os.mkdir(output_dir)

    for file in os.listdir(dir):
        if os.path.splitext(file) != "resized":
            resizeImage(file,output_dir+"/")