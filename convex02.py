# http://scipy.github.io/devdocs/generated/scipy.spatial.ConvexHull.html





import numpy as np
import matplotlib.pyplot as plt
from numpy import ones,zeros,sort


from scipy.spatial import ConvexHull,Delaunay
# points = np.random.rand(4, 2)   # 30 random points in 2-D
# points = np.array([[0,1.530],[1,0.936],[2, 0.649],[3,0.527],[4,0.406],[5,0.212],[6,0.082],[7,0.040],[8,0.029]])
points = np.array([[0,1.020],[1,0.533],[2,0.324],[3, 0.178],[4,0.103],[5,0.057],[6,0.033],[7,0.018],[8,0.007]])
print 'points.shape',points.shape
print 'points'
print points
hull = ConvexHull(points)
print 'hull.points'
print hull.points
print 'hull.simplices'
print hull.simplices
print 'hull.vertices'
print hull.vertices
print 'hull.neighbors'
print hull.neighbors
print 'hull.equations'
print hull.equations

y = np.sort(hull.vertices)


plt.plot(points[:,0], points[:,1], 'o')
for simplex in hull.simplices:
    print simplex
    print points[simplex, 0],points[simplex, 1]
    plt.plot(points[simplex, 0], points[simplex, 1])


plt.show()

