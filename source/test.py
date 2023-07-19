import numpy as np
import matplotlib.pyplot as plt
from dge.fitter_utils import create_mask
from dge.helpers import center_pad
from matplotlib.colors import LogNorm
from scipy.ndimage import affine_transform, zoom
from scipy.optimize import minimize, basinhopping, LinearConstraint
from scipy.fftpack import fftn, ifftn
from scipy import ndimage
from scipy.special import binom
from skimage.draw import polygon
from PIL import Image, ImageDraw
from scipy.ndimage import rotate
from dge import Map, ProteinBlurrer
from scipy.spatial.transform import Rotation as R


center = (15,15)
shape = (31,31)
mask = create_mask("sphere", radius = 5, center = center, shape = shape)
mask2 = create_mask("sphere", radius = 10, center = center, shape = shape)
mask3 = create_mask("sphere", radius = 15, center = center, shape = shape)

# fig, ax = plt.subplots(nrows = 1, ncols = 3, sharex = True, sharey = True)

# ax[0].imshow(np.fft.fftshift(np.abs(np.fft.fftn(mask))), norm=LogNorm())
# ax[1].imshow(np.fft.fftshift(np.abs(np.fft.fftn(mask2))), norm=LogNorm())
# ax[2].imshow(np.fft.fftshift(np.abs(np.fft.fftn(mask3))), norm=LogNorm())
# plt.show()

def find_scaling(f, g):
    """
    Find the ideal scaling between two objects f and g
    using cross-correlation as the similarity measure.

    Parameters
    ----------
    f : array_like
        First input array.
    g : array_like
        Second input array. Should have the same number of dimensions as f.

    Returns
    -------
    ideal_scale : tuple
        The scaling factors that maximize the cross-correlation between f and g.
    """
    F = fftn(f)
    def objective(scale):
        if np.any(scale < 0):
            return 0
        g_scaled = zoom(g, scale)
        g_scaled = center_pad(g_scaled, shape = f.shape)
        G_scaled = fftn(g_scaled)
        correlation = np.real(np.dot(F.flatten(), G_scaled.flatten()))
        correlation /= (np.linalg.norm(F) * np.linalg.norm(G_scaled))
        return -correlation

    initial_guess = np.ones(f.ndim)
    bounds = tuple((0.1, 10) for _ in range(f.ndim))
    linear_constraint = LinearConstraint(
        np.eye(len(bounds)), np.min(bounds, axis=1), np.max(bounds, axis=1)
    )
    result = basinhopping(objective, initial_guess,
                minimizer_kwargs={
                    "method": "COBYLA",
                    "constraints": linear_constraint
                },
    )
    return result


def cost_func(G1_flat, F, D):
    G1 = np.fft.fftn(G1_flat.reshape(F.shape))
    # G1 = G1_flat.reshape(F.shape)
    D1 = F * G1
    error = np.linalg.norm(D - D1)
    error -= np.linalg.norm(mask - G1_flat.reshape(F.shape))
    print(error)
    return error

# Compute F, G and D
# F = fftn(mask)
# G = fftn(mask)
# D = F * G

# initial_guess = np.random.rand(mask.size)
# result = minimize(cost_func, initial_guess, args=(F, D), tol = 1)
# synthesized_object = result.x.reshape(F.shape)
# plt.imshow(synthesized_object);plt.show()


def bezier_curve(control_points : np.ndarray, resolution : int):
    number_points = control_points.shape[0]
    sampling = np.linspace(0, 1, resolution)

    coefficient = np.power(sampling, np.arange(number_points)[:, None])
    coefficient *= np.power(1 - sampling, np.arange(number_points)[::-1, None])

    factorials = np.array([np.math.factorial(x) for x in range(number_points)])
    scaling = np.math.factorial(number_points - 1) / (factorials * factorials[::-1])
    coefficient *= scaling[:, None]
    return control_points.T @ coefficient



def rasterize_curve(curve, shape):
    coordinates = np.flip(curve.T, axis = 1)
    coordinates -= coordinates.min(axis = 0)
    scaling = np.divide(coordinates.max(axis = 0),
        np.subtract(shape, 3)
    )
    scaling[scaling < 1] = 1
    coordinates /= scaling
    coordinates = coordinates.astype(int)
    image = Image.new('L', shape, 0)
    ImageDraw.Draw(image).polygon(
        [tuple(p) for p in coordinates], outline=1, fill=1
    )
    arr = np.array(image)

    coordinates = np.array(np.where(arr > 0))
    coordinates_mod = coordinates.copy()
    axis_max = coordinates_mod.max(axis=1)
    axis_min = coordinates_mod.min(axis=1)
    axis_difference = axis_max - axis_min
    half_box_size = np.array(arr.shape) // 2
    new_origin = half_box_size - axis_max + (axis_difference // 2)
    coordinates_mod += new_origin[:, None]

    out = np.zeros(shape, dtype=np.float32)

    in_box = np.logical_and(
        coordinates_mod < np.array(out.shape)[:, None], coordinates_mod >= 0
    ).min(axis=0)
    coordinates = coordinates[:, in_box]
    coordinates_mod = coordinates_mod[:, in_box]

    coordinates = tuple(coordinates)
    coordinates_mod = tuple(coordinates_mod)
    np.add.at(out, coordinates_mod, arr[coordinates])

    return out


def get_rotation_matrix(angle_in_degrees):
    rotation = R.from_rotvec(np.radians(angle_in_degrees))
    rotation_matrix = rotation.as_matrix()
    return rotation_matrix

def create_volume(curve, shape):
    # Rasterize the 2D curve
    arr = rasterize_curve(curve, shape)

    max_dim = np.max(arr.shape)
    volume = np.zeros((max_dim, *arr.shape), dtype = np.float32)
    volume[max_dim // 2, : ,:] = arr

    image_map = Map(volume, apix = 1, origin = (0,0,0))
    image_map.pad(new_shape = np.repeat(max_dim, 3))

    output_map = image_map.empty
    for angle in range(0, 361, 5):
        rotmat = get_rotation_matrix((0, angle, 0))
        output_map.fullMap = output_map + image_map.rotate(
            rotmat = rotmat,
            center_rotation = True,
            threshold = 0
        )

    return output_map, image_map

def energy_difference(points_flat, input_fft_abs, shape, resolution):
    print(points_flat)

    points = points_flat.reshape(-1, 2)
    points = np.append(points, [points[0]], axis=0)

    curve = bezier_curve(points, resolution)
    synth_volume, _ = create_volume(curve, shape)
    synth_volume.fullMap[synth_volume.fullMap > 0] = 1

    synth_fft = np.fft.fftn(synth_volume.fullMap)

    energy_diff = np.sum(
        np.abs(input_fft)**2) - np.sum(np.abs(synth_fft)**2)/input_fft.size

    return energy_diff


shape, nPoints, resolution = (31, 31), 3, 5000
xmin, xmax = 0, shape[0]

points = np.random.rand(nPoints, 2) * 100
points = np.random.uniform(xmin, xmax, size = (nPoints, 2))
points = np.append(points, [points[0]], axis=0)
input_volume = create_mask("sphere", radius = 15,
    center = (15, 15, 15), shape = (31, 31, 31)
)


points_flat = points[:-1].flatten()

input_fft = np.fft.fftn(input_volume)
input_fft_abs = np.abs(input_fft)

# result = minimize(energy_difference,
#     points_flat, args=(input_fft_abs, shape, resolution))

bounds = [(xmin, xmax)] * len(points_flat)
linear_constraint = LinearConstraint(
    np.eye(len(bounds)), np.min(bounds, axis=1), np.max(bounds, axis=1)
)
result = basinhopping(
    x0 = points_flat,
    func=lambda x: energy_difference(
        x,
        input_fft_abs,
        shape,
        resolution,
    ),
    minimizer_kwargs={"method": "COBYLA", "constraints": linear_constraint},
    niter=100
)


optimized_points = result.x.reshape(-1, 2)
optimized_curve = bezier_curve(
    np.append(optimized_points, [optimized_points[0]], axis=0),
    resolution
)
output_map, image_map = create_volume(optimized_curve, shape)
output_map.fullMap[output_map.fullMap > 0] = 255

# curve = bezier_curve(points, 1000)

blurrer = ProteinBlurrer()
# output_map, image_map = create_volume(curve, size=20, angle=5, axes=(1, 0))
output_map.fullMap = blurrer.gaussian_blur(
    template = output_map.fullMap, apix = 1, sigma = 1
)
output_map.write_map("vol.mrc")
image_map.write_map("curve.mrc")


