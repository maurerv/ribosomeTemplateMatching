import textwrap
from os import makedirs
from os.path import join
from copy import deepcopy

import numpy as np
from tme import Density
from tme.matching_utils import create_mask

BASEDIR = "/g/kosinski/vmaurer/ribosomePaper/templates_pytme_2k"
TARGET_PATH = "/g/kosinski/vmaurer/ribosomePaper/targets/TS_037_inverted.mrc"
sampling_rate = 15
BASEDIR = "/Users/vmaurer/src/ribosomeSpheres/templates_pytme"
TARGET_PATH = "/Users/vmaurer/Downloads/tomograms/TS_037.mrc"


# def make_sbatch(target_path : str,
#     template_path : str, template_mask_path : str, output_path : str,
#     angular_sampling : float) -> str:
#     return textwrap.dedent(f"""\
#         #!/bin/bash
#         #SBATCH --mem 60G
#         #SBATCH -p gpu-el8
#         #SBATCH -N 1
#         #SBATCH -t 12:00:00
#         #SBATCH -C gpu=3090
#         #SBATCH --gres=gpu:1
#         #SBATCH --export=NONE
#         #$BATCH --qos=normal
#         conda activate dge
#         match_template.py \
#             -m {target_path} \
#             -i {template_path} \
#             --template_mask {template_mask_path} \
#             --cutoff_template 0 \
#             -s FLC \
#             -n 1 \
#             --use_gpu \
#             -a {angular_sampling} \
#             --memory_scaling 0.7 \
#             --no_fourier_padding \
#             -o {output_path}
#     """)


def make_sbatch(
    target_path: str,
    template_path: str,
    template_mask_path: str,
    output_path: str,
    angular_sampling: float,
) -> str:
    return textwrap.dedent(
        f"""\
        #!/bin/bash
        #SBATCH --mem 150G
        #SBATCH -p htc-el8
        #SBATCH -N 1
        #SBATCH --cpus-per-task=10
        #SBATCH --qos=low
        #SBATCH --ntasks=1
        #SBATCH -t 48:00:00

        module load Anaconda3

        source activate dge
        match_template.py \
            -m {target_path} \
            -i {template_path} \
            --template_mask {template_mask_path} \
            --cutoff_template 0 \
            -s FLCSphericalMask \
            -n 10 \
            -a {angular_sampling} \
            --no_fourier_padding \
            -o {output_path}
    """
    )


if __name__ == "__main__":
    original = Density.from_structure(
        "../templates/7p6z.cif",
        sampling_rate=13.481,
        filter_by_residues=None,
    )
    original.pad(new_shape=(51, 51, 51))

    radius_lower, radius_upper, run_files = 1, 20, []
    # for radius in range(radius_lower, radius_upper):
    #     template = create_mask(
    #         mask_type = "ellipse",
    #         shape = (51, 51, 51),
    #         center = (25, 25, 25),
    #         radius = radius
    #     )
    #     center = Density.center_of_mass(template, 0).astype(int)
    #     # Alternativel this should be mask_type = "box" for cube masks.
    #     mask = create_mask(
    #         mask_type = "ellipse",
    #         shape = (51, 51, 51),
    #         center = center,
    #         radius = radius + 2
    #     )

    #     template = template.astype(np.float32)
    #     mask = mask.astype(np.float32)

    #     outdir = join(BASEDIR, f"sphere_{radius}")
    #     makedirs(outdir, exist_ok = True)

    #     temp = Density(template, sampling_rate = 13.481, origin = (0,0,0))
    #     temp.origin = deepcopy(original.origin)
    #     temp.to_file(
    #         join(outdir, "template.mrc")
    #     )
    #     Density(
    #         mask, sampling_rate = 13.481, origin = deepcopy(original.origin)
    #     ).to_file(join(outdir, "templateMask.mrc"))

    #     sbatch_file = join(outdir, "submit_job.sbatch")
    #     with open(sbatch_file, "w", encoding = "utf-8") as ofile:
    #         ofile.write(make_sbatch(
    #             target_path = TARGET_PATH,
    #             template_path = join(outdir, "template.mrc"),
    #             template_mask_path = join(outdir, "templateMask.mrc"),
    #             output_path = f"{outdir}/output.pickle",
    #             angular_sampling = sampling_rate,
    #         ))
    #     run_files.append(f"sbatch {sbatch_file}\n")

    # for radius in range(radius_lower, radius_upper):
    #     template = Density.from_file("../templates/volume.mrc")
    #     resampling_rate = np.min(np.divide(template.shape, 2*radius))
    #     template = template.resample(resampling_rate)
    #     template.pad(new_shape = (51, 51, 51))
    #     template.origin = deepcopy(original.origin)
    #     template.sampling_rate = (13.481, 13.481, 13.481)

    #     mask = create_mask(
    #         mask_type = "ellipse",
    #         shape = (51, 51, 51),
    #         center = (25,25,25),
    #         radius = radius + 2
    #     )
    #     mask = mask.astype(np.float32)

    #     outdir = join(BASEDIR, f"emoji_{radius}")
    #     makedirs(outdir, exist_ok = True)

    #     Density(mask, sampling_rate = 13.481,
    #         origin = template.origin).to_file(
    #         join(outdir, "templateMask.mrc")
    #     )
    #     template.to_file(join(outdir, "template.mrc"))
    #     template = template.data

    #     sbatch_file = join(outdir, "submit_job.sbatch")
    #     with open(sbatch_file, "w", encoding = "utf-8") as ofile:
    #         ofile.write(make_sbatch(
    #             target_path = TARGET_PATH,
    #             template_path = join(outdir, "template.mrc"),
    #             template_mask_path = join(outdir, "templateMask.mrc"),
    #             output_path = f"{outdir}/output.pickle",
    #             angular_sampling = sampling_rate,
    #         ))
    #     run_files.append(f"sbatch {sbatch_file}\n")

    # original = Density.from_file("../templates/emd_3228.map.gz")
    # for radius in range(radius_lower, radius_upper):
    #     resampling_rate = 20 * 13.481 / (2 * radius)
    #     template = original.resample(resampling_rate)
    #     template.pad(new_shape = (51, 51, 51))

    #     template.origin = deepcopy(original.origin)
    #     template.sampling_rate = (13.481, 13.481, 13.481)
    #     outdir = join(BASEDIR, f"emd3228_{radius}")
    #     makedirs(outdir, exist_ok = True)

    #     center = Density.center_of_mass(template.data).astype(int)
    #     mask = create_mask(
    #         mask_type = "ellipse",
    #         shape = (51, 51, 51),
    #         center = center,
    #         radius = radius + 2
    #     )
    #     mask = mask.astype(np.float32)

    #     template.to_file(join(outdir, "template.mrc"))
    #     template = template.data
    #     Density(mask, sampling_rate = 13.481, origin = (0,0,0)).to_file(
    #         join(outdir, "templateMask.mrc")
    #     )

    #     sbatch_file = join(outdir, "submit_job.sbatch")
    #     with open(sbatch_file, "w", encoding = "utf-8") as ofile:
    #         ofile.write(make_sbatch(
    #             target_path = TARGET_PATH,
    #             template_path = join(outdir, "template.mrc"),
    #             template_mask_path = join(outdir, "templateMask.mrc"),
    #             output_path = f"{outdir}/output.pickle",
    #             angular_sampling = sampling_rate,
    #         ))
    #     run_files.append(f"sbatch {sbatch_file}\n")

    original = Density.from_structure(
        "../templates/HA_ranked_0.pdb", sampling_rate=13.481, filter_by_residues=None
    )
    original.pad(new_shape=(51, 51, 51))
    original.origin = (0, 0, 0)

    for radius in range(radius_lower, radius_upper):
        resampling_rate = 20 * 13.481 / (2 * radius)
        template = original.resample(resampling_rate)
        template.pad(new_shape=(51, 51, 51))

        template.origin = deepcopy(original.origin)
        template.sampling_rate = (13.481, 13.481, 13.481)
        outdir = join(BASEDIR, f"ha_{radius}")
        makedirs(outdir, exist_ok=True)

        center = Density.center_of_mass(template.data, 0).astype(int)
        mask = create_mask(
            mask_type="ellipse",
            shape=(51, 51, 51),
            center=(25, 25, 25),
            radius=radius + 2,
        )
        mask = mask.astype(np.float32)

        template.to_file(join(outdir, "template.mrc"))
        template = template.data
        Density(mask, sampling_rate=13.481, origin=(0, 0, 0)).to_file(
            join(outdir, "templateMask.mrc")
        )

        sbatch_file = join(outdir, "submit_job.sbatch")
        with open(sbatch_file, "w", encoding="utf-8") as ofile:
            ofile.write(
                make_sbatch(
                    target_path=TARGET_PATH,
                    template_path=join(outdir, "template.mrc"),
                    template_mask_path=join(outdir, "templateMask.mrc"),
                    output_path=f"{outdir}/output.pickle",
                    angular_sampling=sampling_rate,
                )
            )
        run_files.append(f"sbatch {sbatch_file}\n")

    run_script = join(BASEDIR, "../source/runJobs.sh")
    with open(run_script, "w", encoding="utf-8") as ofile:
        ofile.writelines(run_files)
