import textwrap
from os import makedirs
from os.path import join
from copy import deepcopy

import numpy as np
from dge import Map
from dge.fitter_utils import create_mask

BASEDIR = "/g/kosinski/vmaurer/ribosomePaper/templates"
TARGET_PATH = "/g/kosinski/vmaurer/ribosomePaper/targets/TS_037.mrc"
BASEDIR = "/Users/vmaurer/src/ribosomeSpheres/templates"

def make_jobxml(target : str, template : str, template_mask : str,
    outdir : str, angles : str = "angles_90_2.em") -> str:
    return textwrap.dedent(f"""\
        <JobDescription Destination="{outdir}" ID="0" Members="1">
                <Volume Sampling="[0, 0, 0]" Subregion="[0, 0, 0, 0, 0, 0]" Binning="[0, 0, 0]" Filename="{target}"/>
          <Reference PreWedge="" File="{template}" Weighting="">
          </Reference>
          <Mask Filename="{template_mask}" Binning="1" isSphere="True"/>
          <SingleTiltWedge Smooth="0.0" Angle1="30" CutoffRadius="0.0" Angle2="30">
            <TiltAxisRotation Z1="0.0" Z2="0.0" X="0.0"/>
          </SingleTiltWedge>
          <Angles Type="FromEMFile" File="{angles}"/>
          <Score Type="FLCFScore" Value="-10000000000.0">
            <PeakPrior Smooth="-1.0" Radius="0.0" Filename=""/>
          </Score>
          <BandPassFilter LowestFrequency="3.0" Smooth="0.0" HighestFrequency="15.0"/>
        </JobDescription>
    """)

def make_sbatch(outpath : str, xml_path : str) -> str:
    return textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --mem 200G
        #SBATCH -p htc-el8
        #SBATCH -N 1
        #SBATCH -n 6
        #SBATCH -t 48:00:00
        #SBATCH -o {outpath}/slurm.%N.%j.out
        #SBATCH -e {outpath}/slurm.%N.%j.err
        module purge
        module load PyTom/1.0b-foss-2021a-CUDA-11.3.1
        time mpirun -n 6 ` which localization.py ` -j {xml_path} -x 3 -y 3 -z 3
    """)

def write_em(path_to_emfile, value):
    """Function that writes .em file (to the tom format).
    :param path_to_emfile: file path to the .em file
    :param header: dictionary containing header information
    :param value: numpy array containing data
    :return: None
    """
    with open(path_to_emfile, "wb") as f:
        np.array([6], dtype = np.int8).tofile(f)
        np.array([0], dtype = np.int8).tofile(f)
        np.array([0], dtype = np.int8).tofile(f)
        np.array([5], dtype = np.int8).tofile(f)
        np.array(value.shape[::-1], dtype=np.int32).tofile(f)
        np.zeros(496, dtype = np.int8).tofile(f)
        value.tofile(f)

if __name__ == "__main__":

    original = Map.from_structure(
        "../templates/7p6z.cif",
        apix = 13.481,
        keep_residues = None,
    )
    original.pad(new_shape = (51, 51, 51))

    radius_lower, radius_upper, run_files = 1, 20, []
    for radius in range(radius_lower, radius_upper):
        template = create_mask(
            mask_type = "sphere",
            shape = (51, 51, 51),
            center = (25, 25, 25),
            radius = radius
        )
        mask = np.ones_like(template)
        outdir = join(BASEDIR, f"sphere_{radius}")
        makedirs(outdir, exist_ok = True)
        temp = Map(template, apix = 13.481, origin = (0,0,0))
        temp.origin = deepcopy(original.origin)
        temp.write_map(
            join(outdir, "template.mrc")
        )

        template_path = join(outdir, "template.em")
        templateMask_path = join(outdir, "templateMask.em")
        write_em(template_path, template)
        write_em(templateMask_path, mask)

        xml_path = join(outdir, "TM_job.xml")
        with open(xml_path, "w", encoding = "utf-8") as ofile:
            ofile.write(make_jobxml(
                target = TARGET_PATH,
                template = template_path,
                template_mask = templateMask_path,
                outdir = outdir
            ))
        sbatch_file = join(outdir, "TM_submit_CPU.sbatch")
        with open(sbatch_file, "w", encoding = "utf-8") as ofile:
            ofile.write(make_sbatch(
                outpath = outdir,
                xml_path = xml_path
            ))
        run_files.append(f"sbatch {sbatch_file}\n")

    for radius in range(radius_lower, radius_upper):
        template = Map.from_file("../templates/volume.mrc")
        resampling_rate = np.min(np.divide(template.box_size, 2*radius))
        template = template.resample(resampling_rate)
        template.pad(new_shape = (51, 51, 51))
        template.origin = deepcopy(original.origin)
        template.apix = 13.481
        template.write_map(join(outdir, "template.mrc"))
        template = template.fullMap

        mask = np.ones_like(template)
        outdir = join(BASEDIR, f"emoji_{radius}")
        makedirs(outdir, exist_ok = True)

        template_path = join(outdir, "template.em")
        templateMask_path = join(outdir, "templateMask.em")
        write_em(template_path, template)
        write_em(templateMask_path, mask)

        xml_path = join(outdir, "TM_job.xml")
        with open(xml_path, "w", encoding = "utf-8") as ofile:
            ofile.write(make_jobxml(
                target = TARGET_PATH,
                template = template_path,
                template_mask = templateMask_path,
                outdir = outdir,
                angles = "angles_25.25_980.em"
            ))
        sbatch_file = join(outdir, "TM_submit_CPU.sbatch")
        with open(sbatch_file, "w", encoding = "utf-8") as ofile:
            ofile.write(make_sbatch(
                outpath = outdir,
                xml_path = xml_path
            ))
        run_files.append(f"sbatch {sbatch_file}\n")

    for radius in range(radius_lower, radius_upper):
        resampling_rate = 20 * 13.481 / (2 * radius)
        template = Map.from_structure(
            "../templates/7p6z.cif",
            apix = resampling_rate,
            keep_residues = None,
        )
        template.pad(new_shape = (51, 51, 51))

        template.origin = deepcopy(original.origin)
        template.apix = 13.481
        outdir = join(BASEDIR, f"7p6z_{radius}")
        makedirs(outdir, exist_ok = True)

        template.write_map(join(outdir, "template.mrc"))
        template = template.fullMap

        mask = np.ones_like(template)

        template_path = join(outdir, "template.em")
        templateMask_path = join(outdir, "templateMask.em")
        write_em(template_path, template)
        write_em(templateMask_path, mask)

        xml_path = join(outdir, "TM_job.xml")
        with open(xml_path, "w", encoding = "utf-8") as ofile:
            ofile.write(make_jobxml(
                target = TARGET_PATH,
                template = template_path,
                template_mask = templateMask_path,
                outdir = outdir,
                angles = "angles_25.25_980.em"
            ))
        sbatch_file = join(outdir, "TM_submit_CPU.sbatch")
        with open(sbatch_file, "w", encoding = "utf-8") as ofile:
            ofile.write(make_sbatch(
                outpath = outdir,
                xml_path = xml_path
            ))
        run_files.append(f"sbatch {sbatch_file}\n")

    run_script = join(BASEDIR, "../source/runJobs.sh")
    with open(run_script, "w", encoding = "utf-8") as ofile:
        ofile.writelines(run_files)
