"""Use Maya's shipped OCIO Python binding, without initializing Maya or rendering."""
from pathlib import Path
import subprocess
import sys

import numpy as np
import PyOpenColorIO as ocio

folder = Path(sys.argv[1])
config = ocio.Config.CreateFromFile('C:/Program Files/Autodesk/Maya2027/resources/OCIO-configs/Maya2022-default/config.ocio')
transform = ocio.DisplayViewTransform(src='ACEScg', display='sRGB', view='ACES 1.0 SDR-video')
processor = config.getProcessor(transform).getDefaultCPUProcessor()
for name in ('beauty', 'beauty_raw', 'specular', 'transmission', 'specular_indirect'):
    with (folder / (name + '.pfm')).open('rb') as stream:
        assert stream.readline().strip() == b'PF'
        width, height = map(int, stream.readline().split())
        scale = float(stream.readline())
        pixels = np.frombuffer(stream.read(), dtype='<f4' if scale < 0 else '>f4').copy()
    pixels = np.ascontiguousarray(pixels.reshape(height, width, 3))
    processor.apply(ocio.PackedImageDesc(pixels, width, height, 3))
    target = folder / (name + '_display.pfm')
    with target.open('wb') as stream:
        stream.write(f'PF\n{width} {height}\n-1.0\n'.encode())
        stream.write(pixels.astype('<f4').tobytes())
    subprocess.run(['C:/Program Files/Autodesk/Arnold/Maya2027/bin/oiiotool.exe',
                    str(target), '--attrib', 'oiio:ColorSpace', 'sRGB', '-d', 'uint8',
                    '-o', str(folder / (name + '_view.png'))], check=True)
print('OCIO_VIEW_PASS', ocio.__version__, 'ACEScg -> ACES 1.0 SDR-video (sRGB)')
