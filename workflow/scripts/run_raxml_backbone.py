import subprocess
import sys

def run_raxml_backbone(alignment, tree, model, log_file):
    """
    Run the raxml-ng command with a given alignment, tree, and model, and handle specific error cases.
    Parameters:
    alignment (str): Path to the alignment file.
    tree (str): Path to the tree file to be used as a constraint.
    model (str): Substitution model to be used.
    log_file (str): Path to the log file where the output will be written.
    Returns:
    int: 0 if the initial command was successful, 1 if branch length optimization was run due to a specific error.
    """
    # Run the initial raxml-ng command
    cmd = [
        "raxml-ng",
        "--redo",
        "--msa", alignment,
        "--model", model,
        "--tree-constraint", tree,
        "--search"
    ]
    with open(log_file, "w") as log:
        result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    
    # Check if the log file contains the specific error message
    with open(log_file, "r") as log:
        log_content = log.read()
        if "ERROR: You provided a comprehensive, fully-resolved tree as a topological constraint." in log_content:
            print("The provided tree is fully resolved. Running branch length optimization instead.")
            # Run the branch length optimization
            cmd = [
                "raxml-ng",
                "--evaluate",
                "--msa", alignment,
                "--model", model,
                "--tree", tree,
                "--brlen", "scaled"
            ]
            with open(log_file, "a") as log:
                log.write("\nRunning branch length optimization...\n")
                result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
            return 1  # Indicate that the branch length optimization was run
    return 0  # Indicate that the initial command was successful

if __name__ == "__main__":
    alignment = sys.argv[1]
    tree = sys.argv[2]
    model = sys.argv[3]
    log_file = sys.argv[4]
    run_raxml_backbone(alignment, tree, model, log_file)