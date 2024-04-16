import numpy as np
from scipy.spatial import ConvexHull


def inhull(X: np.ndarray,
           hull_pts: np.ndarray,
           eps: float = np.finfo(float).eps
          )->np.ndarray:
  
  '''
  Check for all points in X whether they lie in the convex hull defined by hull_pts.
  Adapted from https://stackoverflow.com/questions/31404658/check-if-points-lies-inside-a-convex-hull

  Parameters
  -----------
  X: np.ndarray
    An array of observations. Must be of shape (n_samples, n_dimensions).
  
  hull_pts: np.ndarray
    An array of points from which the convex hull is determined. Must be of shape (n_points, n_dimensions).
  
  eps: np.float64
    The tolerance to be used when checking whether a given point is inside the hull.
    Choose > 0 to avoid numerical issues, default = np.finfo(float).eps

  Returns
  --------
  in_hull: np.ndarray
    A boolean array of shape (n_samples,) indicating which points in X are inside the hull
    
  '''

  if X.shape[1] != hull_pts.shape[1]:
    raise ValueError(f'Hull points and test points must have the same dimensions.')
  
  hull = ConvexHull(hull_pts)

  # A is shape (f, d) and b is shape (f, 1).
  A, b = hull.equations[:, :-1], hull.equations[:, -1]

  # The hull is defined as all points x for which Ax + b <= 0.
  in_hull = np.array([np.all(A @ x + b < eps) for x in X])

  return in_hull