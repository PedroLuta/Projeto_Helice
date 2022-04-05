import numpy as nup
import Auxiliary

x = [1, 2, 4]
y = [1, 2, 4]

area = Auxiliary.area_under_curve(x, y)
area_nup = nup.trapz(y, x)

print(area)
print(area_nup)