#!/usr/bin/env python3

# =============================================================================
# CIRCexplorer2 Input List Preparation Pipeline
# =============================================================================
# Purpose:
#     Prepares input list files for the CIRCexplorer2 detection step from STAR
#     mapping outputs. For each sample directory, the script locates the joint
#     chimeric junction file, mate1 chimeric junction file, mate2 chimeric
#     junction file, and coordinate-sorted BAM file, converts each path to the
#     Docker-mounted path format, and writes the corresponding list files.
#
# Inputs:
#     mapping_dir          :            Directory containing sample-level STAR
#                                       mapping output folders.
#     out_dir              :            Directory where CIRCexplorer2 input list
#                                       files will be written.
#
# Outputs:
#     out_dir/samplesheet  :            Docker paths to joint Chimeric.out.junction
#                                       files, one sample per line.
#     out_dir/mate1        :            Docker paths to mate1 Chimeric.out.junction
#                                       files, one sample per line.
#     out_dir/mate2        :            Docker paths to mate2 Chimeric.out.junction
#                                       files, one sample per line.
#     out_dir/bamlist      :            Docker paths to coordinate-sorted STAR
#                                       BAM files, one sample per line.
# =============================================================================

import sys
from pathlib import Path

DOCKER_PREFIX = "/host_os" # needed if circtools run by their docker

# =============================================================================
# Function: pick_junction
# =============================================================================
# Selects the STAR chimeric junction file from a sample or mate-specific output
# folder. The function accepts either the exact STAR filename or one matching
# prefixed variant, and fails loudly if the folder contains none or more than
# one possible match.
#
# Parameters:
#     folder               :            Directory expected to contain a STAR
#                                       Chimeric.out.junction file.
#
#
# Returns:
#     Path                 :            Path to the selected chimeric junction file.
# =============================================================================

def pick_junction(folder: Path):
    """
    Pick Chimeric.out.junction file.
    Accepts either exact name or single prefixed variant.
    """
    exact = folder / "Chimeric.out.junction"
    if exact.exists():
        return exact

    matches = [f for f in folder.iterdir()
               if f.name.endswith("Chimeric.out.junction")]

    if len(matches) == 1:
        return matches[0]

    if len(matches) == 0:
        raise FileNotFoundError(f"[ERROR] No junction file in {folder}")

    raise RuntimeError(f"[ERROR] Multiple junction files in {folder}: {matches}")


# =============================================================================
# Function: pick_bam
# =============================================================================
# Selects the coordinate-sorted BAM produced by the joint paired-end STAR
# alignment for a sample. The function requires exactly one matching BAM file.
#
# Parameters:
#     folder               :            Sample-level STAR output directory.
#
#
# Returns:
#     Path                 :            Path to the selected coordinate-sorted BAM.
# =============================================================================

def pick_bam(folder: Path):
    """
    Pick STAR coordinate-sorted BAM from joint mapping.
    """
    matches = list(folder.glob("*Aligned.sortedByCoord.out.bam"))

    if len(matches) == 1:
        return matches[0]

    if len(matches) == 0:
        raise FileNotFoundError(f"[ERROR] No BAM found in {folder}")

    raise RuntimeError(f"[ERROR] Multiple BAMs found in {folder}: {matches}")


# =============================================================================
# Function: to_docker_path
# =============================================================================
# Converts a host filesystem path into the path expected inside the Docker
# container by prefixing the resolved path with the configured host mount point.
#
# Parameters:
#     p                    :            Host filesystem path.
#
#
# Returns:
#     str                  :            Docker-visible path string.
# =============================================================================

def to_docker_path(p: Path) -> str:
    return DOCKER_PREFIX + str(p.resolve())


# =============================================================================
# Function: main
# =============================================================================
# Walks all sample directories in the STAR mapping output folder, collects the
# required joint, mate-specific, and BAM files, converts them to Docker paths,
# and writes the four CIRCexplorer2 input list files.
#
# Parameters:
#     mapping_dir          :            Directory containing sample-level STAR
#                                       mapping output folders.
#     out_dir              :            Directory where list files will be written.
#
#
# Returns:
#     None                 :            Writes samplesheet, mate1, mate2, and
#                                       bamlist files to out_dir.
# =============================================================================

def main(mapping_dir, out_dir):
    mapping = Path(mapping_dir).resolve()
    out = Path(out_dir).resolve()

    if not mapping.is_dir():
        sys.exit(f"[FATAL] Mapping dir not found: {mapping}")

    out.mkdir(parents=True, exist_ok=True)

    samples = sorted([d for d in mapping.iterdir() if d.is_dir()])

    if not samples:
        sys.exit("[FATAL] No sample directories found")

    samplesheet = []
    mate1_list = []
    mate2_list = []
    bam_list = []

    for s in samples:
        print(f"[INFO] {s.name}")

        # Collect the joint and mate-specific chimeric junction files plus the
        # coordinate-sorted BAM generated by the previous STAR mapping step.
        joint = pick_junction(s)
        m1 = pick_junction(s / "mate1")
        m2 = pick_junction(s / "mate2")
        bam = pick_bam(s)

        samplesheet.append(to_docker_path(joint))
        mate1_list.append(to_docker_path(m1))
        mate2_list.append(to_docker_path(m2))
        bam_list.append(to_docker_path(bam))

    # Write one Docker-visible input path per sample for each CIRCexplorer2
    # input category.
    (out / "samplesheet").write_text("\n".join(samplesheet) + "\n")
    (out / "mate1").write_text("\n".join(mate1_list) + "\n")
    (out / "mate2").write_text("\n".join(mate2_list) + "\n")
    (out / "bamlist").write_text("\n".join(bam_list) + "\n")

    print("[INFO] Files written:")
    print(f"  - {out / 'samplesheet'}")
    print(f"  - {out / 'mate1'}")
    print(f"  - {out / 'mate2'}")
    print(f"  - {out / 'bamlist'}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: script <01Mapping_dir> <02Detect/input>")

    main(sys.argv[1], sys.argv[2])
