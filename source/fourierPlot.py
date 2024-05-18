import numpy as np
from tme import Density
import matplotlib.pyplot as plt

templates = (
    "/Users/vmaurer/src/ribosomeSpheres/templates/emoji_10/template.mrc",
    "/Users/vmaurer/src/ribosomeSpheres/templates/sphere_10/template.mrc",
    "/Users/vmaurer/src/ribosomeSpheres/templates/emd3228_10/template.mrc",
)
templates_class = ("emoji", "sphere", "emd3228")

with open("radialAverages.tsv", "w", encoding="utf-8") as ofile:
    ofile.write("class\tbin\tvalue\n")
    for template, template_class in zip(templates, templates_class):
        template = Density.from_file(template).data
        fourier_transform = np.fft.fftshift(np.fft.fftn(template))
        fourier_spectrum = np.abs(fourier_transform)

        half_shape = np.divide(template.shape, 2).astype(int)
        center = np.array(half_shape).reshape((-1,) + (1,) * template.ndim)
        distances = np.linalg.norm(
            np.indices(
                template.shape,
            )
            - center,
            axis=0,
        )

        bins = np.rint(distances).astype(int)
        for i in np.unique(bins):
            ofile.write(
                f"{template_class}\t{i}\t{fourier_spectrum[bins == i].mean()}\n"
            )

colors = dict(zip(templates_class, ("g", "r", "b")))
fig, ax = plt.subplots(nrows=1, ncols=3)
legends = []
all_phases = {}

for template, template_class in zip(templates, templates_class):
    template = Density.from_file(template).data
    fourier_transform = np.fft.fftshift(np.fft.fftn(template))

    mag = np.abs(fourier_transform)
    phase = np.angle(fourier_transform)

    half_shape = np.divide(template.shape, 2).astype(int)
    center = np.array(half_shape).reshape((-1,) + (1,) * template.ndim)
    distances = np.linalg.norm(
        np.indices(
            template.shape,
        )
        - center,
        axis=0,
    )

    bins = np.rint(distances).astype(int)
    unique_bins = np.sort(np.unique(bins))
    unique_bins = unique_bins[unique_bins < np.max(half_shape)]

    mags = np.array([np.mean(mag[bins == i]) for i in unique_bins])
    mags /= mags.max()

    phases = np.array([np.mean(np.abs(phase[bins == i])) for i in unique_bins])
    all_phases[template_class] = phases
    ax[0].plot(unique_bins, mags, c=colors[template_class])
    ax[1].plot(unique_bins, phases, c=colors[template_class])
    legends.append(f"{template_class}")

ax[0].legend(legends, loc="upper right")

# Set plot titles
ax[0].set_title("Mag")
ax[1].set_title("Phase")
ax[2].set_title("Phase Difference")

all_phases["sphere"] -= all_phases["emd3228"]
all_phases["emoji"] -= all_phases["emd3228"]
all_phases["emd3228"] -= all_phases["emd3228"]

for template_class, phases in all_phases.items():
    ax[2].plot(unique_bins, phases, c=colors[template_class])
    legends.append(f"{template_class}")

plt.show()
