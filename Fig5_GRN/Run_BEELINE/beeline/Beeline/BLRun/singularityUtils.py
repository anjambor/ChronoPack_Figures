"""
Utility functions to convert Docker commands to Singularity commands for BEELINE.
"""
from pathlib import Path
import os

# Path to Singularity images
SINGULARITY_IMAGE_DIR = Path("/new-stg/home/alex/projects/grn_inference/BEELINE/singularity_images")

# Mapping from Docker Hub image names to local Singularity image files
DOCKER_TO_SINGULARITY_MAP = {
    "grnbeeline/arboreto:base": "arboreto.sif",
    "grnbeeline/pidc:base": "pidc.sif", 
    "grnbeeline/genie3:base": "arboreto.sif",  # GENIE3 uses same image as arboreto
    "grnbeeline/ppcor:base": "ppcor.sif",
    "grnbeeline/grnvbem:base": "grnvbem.sif",
    "grnbeeline/scode:base": "scode.sif",
    "grnbeeline/scns:base": "scns.sif",
    "grnbeeline/sincerities:base": "sincerities.sif",
    "grnbeeline/leap:base": "leap.sif",
    "grnbeeline/grisli:base": "grisli.sif",  # Note: may need to be pulled separately
    "grnbeeline/singe:base": "singe.sif",
    "grnbeeline/scribe:base": "scribe.sif",
    "grnbeeline/scsgl:base": "scsgl.sif",
    "scinge:base": "scinge.sif",
    "jump3:base": "jump3.sif"
}

def docker_to_singularity_cmd(docker_cmd):
    """
    Convert a Docker command to equivalent Singularity command.
    
    Args:
        docker_cmd (str): Docker command string
        
    Returns:
        str: Equivalent Singularity command
    """
    # Parse the docker command
    parts = docker_cmd.split()
    
    if "docker run" not in docker_cmd:
        return docker_cmd
    
    # Find the docker run part
    run_idx = parts.index("run")
    
    # Extract components
    docker_image = None
    bind_mounts = []
    command_parts = []
    
    i = run_idx + 1
    while i < len(parts):
        if parts[i] == "--rm":
            i += 1
            continue
        elif parts[i] == "-v":
            # Volume mount
            bind_mounts.append(parts[i+1])
            i += 2
            continue
        elif parts[i].startswith("--expose"):
            # Skip expose option (not needed for Singularity)
            i += 1
            continue
        elif parts[i].startswith("--entrypoint"):
            # Skip entrypoint (handle differently in Singularity)
            i += 1
            continue
        elif parts[i].startswith("grnbeeline/") or parts[i].endswith(":base") or parts[i] in DOCKER_TO_SINGULARITY_MAP:
            # This is the docker image
            docker_image = parts[i]
            # Everything after this is the command
            command_parts = parts[i+1:]
            break
        else:
            i += 1
    
    # Get corresponding Singularity image
    if docker_image in DOCKER_TO_SINGULARITY_MAP:
        sif_file = DOCKER_TO_SINGULARITY_MAP[docker_image]
        sif_path = SINGULARITY_IMAGE_DIR / sif_file
    else:
        # Fallback: try to guess from image name
        if "arboreto" in docker_image.lower():
            sif_path = SINGULARITY_IMAGE_DIR / "arboreto.sif"
        elif "pidc" in docker_image.lower():
            sif_path = SINGULARITY_IMAGE_DIR / "pidc.sif"
        else:
            raise ValueError(f"Unknown Docker image: {docker_image}")
    
    # Build Singularity command
    sing_cmd_parts = ["singularity", "exec"]
    
    # Add bind mounts
    for mount in bind_mounts:
        if ":" in mount:
            # Format: host_path:container_path
            host_path, container_path = mount.split(":", 1)
            sing_cmd_parts.extend(["--bind", f"{host_path}:{container_path}"])
        else:
            # Simple mount
            sing_cmd_parts.extend(["--bind", mount])
    
    # Add image path
    sing_cmd_parts.append(str(sif_path))
    
    # Add command
    sing_cmd_parts.extend(command_parts)
    
    return " ".join(sing_cmd_parts)

def get_singularity_image_path(algorithm_name):
    """
    Get the path to the Singularity image for a given algorithm.
    
    Args:
        algorithm_name (str): Name of the algorithm or Docker image name
        
    Returns:
        Path: Path to the Singularity image file
    """
    # Handle Docker image names by mapping them to algorithm names
    docker_to_algorithm = {
        "grnbeeline/leap:base": "LEAP",
        "grnbeeline/grnvbem:base": "GRNVBEM", 
        "scsgl:base": "SCSGL",
        "grnbeeline/scribe:base": "SCRIBE",
        "grnbeeline/singe:0.4.1": "SINGE",
        "grnbeeline/grisli:base": "GRISLI",
        "grnbeeline/scode:base": "SCODE",
        "grnbeeline/scns:base": "SCNS",
        "jump3:base": "JUMP3",
        "grnbeeline/sincerities:base": "SINCERITIES",
        "grnbeeline/ppcor:base": "PPCOR",
        "grnbeeline/genie3:base": "GENIE3"
    }
    
    # Convert Docker image name to algorithm name if needed
    if algorithm_name in docker_to_algorithm:
        algorithm_name = docker_to_algorithm[algorithm_name]
    
    algorithm_to_sif = {
        "GRNBOOST2": "arboreto.sif",
        "GENIE3": "arboreto.sif", 
        "PIDC": "pidc.sif",
        "PPCOR": "ppcor.sif",
        "GRNVBEM": "grnvbem.sif",
        "SCODE": "scode.sif",
        "SCNS": "scns.sif",
        "SINCERITIES": "sincerities.sif",
        "LEAP": "leap.sif",
        "GRISLI": "grisli.sif",  # May need to be pulled separately if missing
        "SINGE": "singe.sif",
        "SCRIBE": "scribe.sif",
        "SCSGL": "scsgl.sif",
        "SCINGE": "scinge.sif",
        "JUMP3": "jump3.sif"
    }
    
    if algorithm_name.upper() in algorithm_to_sif:
        sif_file = algorithm_to_sif[algorithm_name.upper()]
        sif_path = SINGULARITY_IMAGE_DIR / sif_file
        
        # Check if the file exists
        if not sif_path.exists():
            print(f"WARNING: Singularity image {sif_path} not found!")
            print(f"You may need to pull it manually:")
            print(f"  singularity pull {sif_path} docker://grnbeeline/{algorithm_name.lower()}:base")
            
            # For GRISLI, which might not have a pre-built image, suggest alternatives
            if algorithm_name.upper() == "GRISLI":
                print("Note: GRISLI might not have a pre-built image. You may need to:")
                print("1. Build it from the Dockerfile in Beeline/Algorithms/GRISLI/")
                print("2. Or try pulling: singularity pull grisli.sif docker://grnbeeline/grisli:base")
        
        return sif_path
    else:
        raise ValueError(f"Unknown algorithm: {algorithm_name}")

def check_all_images():
    """
    Check which Singularity images are available and which are missing.
    """
    print("Checking Singularity image availability:")
    print("=" * 50)
    
    algorithm_to_sif = {
        "GRNBOOST2": "arboreto.sif",
        "GENIE3": "arboreto.sif", 
        "PIDC": "pidc.sif",
        "PPCOR": "ppcor.sif",
        "GRNVBEM": "grnvbem.sif",
        "SCODE": "scode.sif",
        "SCNS": "scns.sif",
        "SINCERITIES": "sincerities.sif",
        "LEAP": "leap.sif",
        "GRISLI": "grisli.sif",
        "SINGE": "singe.sif",
        "SCRIBE": "scribe.sif",
        "SCSGL": "scsgl.sif",
        "SCINGE": "scinge.sif"
    }
    
    available = []
    missing = []
    
    for algo, sif_file in algorithm_to_sif.items():
        sif_path = SINGULARITY_IMAGE_DIR / sif_file
        if sif_path.exists():
            size = sif_path.stat().st_size / (1024*1024*1024)  # Size in GB
            available.append(f"{algo:12} -> {sif_file:15} ({size:.1f}GB)")
        else:
            missing.append(f"{algo:12} -> {sif_file:15} (MISSING)")
    
    print("Available images:")
    for img in available:
        print(f"  ✓ {img}")
    
    if missing:
        print(f"\nMissing images ({len(missing)}):")
        for img in missing:
            print(f"  ✗ {img}")
        
        print(f"\nTo pull missing images:")
        for algo, sif_file in algorithm_to_sif.items():
            sif_path = SINGULARITY_IMAGE_DIR / sif_file
            if not sif_path.exists():
                print(f"  singularity pull {sif_path} docker://grnbeeline/{algo.lower()}:base")
    else:
        print("\n✓ All images are available!")
    
    return len(available), len(missing) 