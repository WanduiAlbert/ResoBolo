
import numpy as np
import matplotlib.pyplot as pl


a = 6.

y = np.linspace(-6,6,1000)

y0 = y - a/(1 + 4*y*y)

pl.plot(y0,y-y0)
pl.show()
