#! /usr/bin/python


from mayavi import mlab
import numpy as np
import math

X, Y, Z = np.ogrid[-2:2:20j, -2:2:20j, -2:2:20j]
N =  np.sqrt(1/math.pi)
print(N)
s = N*np.exp(-np.sqrt(X**2 + Y**2 + Z**2))

mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(s),
                                    plane_orientation='x_axes',
                                                                slice_index=20,
                                                                                        )
mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(s),
                                    plane_orientation='y_axes',
                                                                slice_index=20,
                                                                                        )
mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(s),
                                    plane_orientation='z_axes',
                                                                slice_index=20,
                                                                                        )
mlab.outline()
mlab.show()


