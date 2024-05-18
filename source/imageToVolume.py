from tme import Density, Preprocessor
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R


def get_rotation_matrix(angle_in_degrees):
    rotation = R.from_rotvec(np.radians(angle_in_degrees))
    rotation_matrix = rotation.as_matrix()
    return rotation_matrix


# 20, 32, 40, 48, 64, 96, 160 available
font_size = 160
char = chr(0x1F4A9)
char = "\u2665"
font = ImageFont.truetype("/System/Library/Fonts/Apple Color Emoji.ttc", size=font_size)

img = Image.new("L", (font_size, font_size), color="black")

draw = ImageDraw.Draw(img)
draw.text((0, 0), char, fill="white", font=font)


img = np.array(img)
fig, ax = plt.subplots(nrows=1, ncols=2)
ax[0].imshow(img)
# img[img < 40] = 0
ax[1].imshow(img)
plt.show()
output_volume = np.zeros((font_size, font_size, font_size))
output_volume[font_size // 2, :, :] = np.array(img)
output_volume = output_volume.astype(np.float32)
image_map = Density(data=output_volume, sampling_rate=1, origin=0)
image_map, _ = image_map.centered(0)
output_map = image_map.empty

# image_map.write_map("image.mrc")

for angle in range(0, 361, 1):
    rotmat = get_rotation_matrix((0, angle, 0))
    output_map.data = (
        output_map.data
        + image_map.rigid_transform(
            rotation_matrix=rotmat,
            use_geometric_center=False,
            order=1,
        ).data
    )
blurrer = Preprocessor()
output_map.data[output_map.data > 0] = 1
output_map.data = blurrer.gaussian_filter(template=output_map.data, sigma=1)
print(output_map.shape)
# output_map.to_file("../templates/volume.mrc")
output_map = output_map.resample(8)
output_map.sampling_rate = 13.481
output_map.pad(new_shape=(51, 51, 51))
# output_map.to_file("../templates/emoji.mrc")
