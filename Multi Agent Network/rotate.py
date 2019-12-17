def rotate(pivot, p, rot):
    x1, y1, z1 = pivot
    x, y, z = p
    x = x - x1
    y = y - y1
    z = z - z1
    xdelta, ydelta, zdelta = rot
    if xdelta == -1 and ydelta == 1:
      return (x1 + y,y1 - x,z1 + z)
    if xdelta == 1 and ydelta == -1:
      return (x1 - y,y1 + x,z1 + z)
    if xdelta == -1 and zdelta == 1:
      return (x1 + z,y1 + y,z1 - x)
    if xdelta == 1 and zdelta == -1:
      return (x1 - z,y1 + y,z1 + x)
    if ydelta == -1 and zdelta == 1:
      return (x1 + x,y1 + z,z1 - y)
    if ydelta == 1 and zdelta == -1:
      return (x1 + x,y1 - z,z1 + y)
    return(-x, -y, z)