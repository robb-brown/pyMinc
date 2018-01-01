import pyMinc
import pylab
pylab.ion()

im = pyMinc.VIOVolume('00100.mnc')

print('Image shape is: {}'.format(im.data.shape))
pylab.figure(); pylab.imshow(im.data[im.data.shape[0]//2])
