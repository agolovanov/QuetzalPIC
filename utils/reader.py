import h5py as _h5py
import os as _os
import numpy as _np

class Simulation:
    def __init__(self, path):
        if not _os.path.exists(path):
            raise ValueError(f"Path [{path}] does not exist")

        if _os.path.isdir(path):
            self.dir = path
        else:
            raise ValueError(f"Path [{path}] is not a directory")

    def read_field(self, field_name, iteration=0, plane='xy', include_mesh=False):
        filename = _os.path.join(self.dir, f'Fields_{plane}_{iteration:03d}.h5')
        with _h5py.File(filename, 'r') as f:
            field = f[field_name]
            result = field[()].T

            if include_mesh:
                ny, nx = result.shape

                x0, y0 = field.attrs['origin']
                dx, dy = field.attrs['steps']

                x = _np.linspace(x0, x0 + (nx - 1) * dx, nx)
                y = _np.linspace(y0, y0 + (ny - 1) * dy, ny)
                
                xx, yy = _np.meshgrid(x, y)
                return xx, yy, result
            else:
                return result