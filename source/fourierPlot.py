import numpy as np
from dge import Map

templates = (
    "/Users/vmaurer/src/ribosomeSpheres/templates/emoji_10/template.mrc",
    "/Users/vmaurer/src/ribosomeSpheres/templates/sphere_10/template.mrc",
    "/Users/vmaurer/src/ribosomeSpheres/templates/7p6z_10/template.mrc",
    )
templates_class = ("emoji", "sphere", "ribosome")

with open("radialAverages.tsv", "w", encoding = "utf-8") as ofile:
    ofile.write("class\tbin\tvalue\n")
    for template, template_class in zip(templates, templates_class):
        template = Map.from_file(template).fullMap
        fourier_transform = np.fft.fftshift(np.fft.fftn(template))
        fourier_spectrum = np.abs(fourier_transform)

        half_shape = np.divide(template.shape, 2).astype(int)
        center = np.array(half_shape).reshape((-1,) + (1,) * template.ndim)
        distances = np.linalg.norm(np.indices(template.shape,) - center, axis=0)

        bins = np.rint(distances).astype(int)
        radial_averages = [
            fourier_spectrum[bins == i].mean() for i in np.unique(bins)
        ]
        for i in np.unique(bins):
            ofile.write(
                f"{template_class}\t{i}\t{fourier_spectrum[bins == i].mean()}\n"
            )

    # mask = np.ones((20,20,20))
    # template = Map(mask, apix = 1, origin = (0,0,0))
    # template.pad((51,51,51))
    # template = template.fullMap
    # fourier_transform = np.fft.fftshift(np.fft.fftn(template))
    # fourier_spectrum = np.abs(fourier_transform)

    # half_shape = np.divide(template.shape, 2).astype(int)
    # center = np.array(half_shape).reshape((-1,) + (1,) * template.ndim)
    # distances = np.linalg.norm(np.indices(template.shape,) - center, axis=0)

    # bins = np.rint(distances).astype(int)
    # radial_averages = [
    #     fourier_spectrum[bins == i].mean() for i in np.unique(bins)
    # ]
    # for i in np.unique(bins):
    #     ofile.write(
    #         f"rectangle\t{i}\t{fourier_spectrum[bins == i].mean()}\n"
    #     )