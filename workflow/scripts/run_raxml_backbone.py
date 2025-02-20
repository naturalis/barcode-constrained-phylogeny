import subprocess
import sys
import os

def run_raxml_backbone(alignment, tree, model, log_file):
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
        log.write("Running initial raxml-ng command:\n")
        log.write(" ".join(cmd) + "\n")
        try:
            result = subprocess.run(cmd, stdout=log, stderr=log, check=True)
        except subprocess.CalledProcessError as e:
            log.write(f"\nCommand failed with return code {e.returncode}\n")
            log.write(f"Command output: {e.output}\n")
            return 2  # Indicate that an error occurred

    # Check if the log file contains the specific error message
    with open(log_file, "r") as log:
        log_content = log.read()
        if "ERROR: You provided a comprehensive, fully-resolved tree as a topological constraint." in log_content:
            # Run the branch length optimization
            cmd = [
                "raxml-ng",
                "--evaluate",
                "--redo",
                "--msa", alignment,
                "--model", model,
                "--tree", tree,
                "--brlen", "scaled"
            ]
            with open(log_file, "a") as log:
                log.write("\nRunning branch length optimization:\n")
                log.write(" ".join(cmd) + "\n")
                try:
                    result = subprocess.run(cmd, stdout=log, stderr=log, check=True)
                except subprocess.CalledProcessError as e:
                    log.write(f"\nCommand failed with return code {e.returncode}\n")
                    log.write(f"Command output: {e.output}\n")
                    return 2  # Indicate that an error occurred
            return 1  # Indicate that the branch length optimization was run

    # Check if the expected output file is created
    output_file = alignment + ".raxml.bestTree"
    if not os.path.exists(output_file):
        with open(log_file, "a") as log:
            log.write(f"\nExpected output file {output_file} not found.\n")
        return 3  # Indicate that the output file is missing

    return 0  # Indicate that the initial command was successful

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python run_raxml_backbone.py <alignment> <tree> <model> <log_file>")
        sys.exit(1)
    
    alignment = sys.argv[1]
    tree = sys.argv[2]
    model = sys.argv[3]
    log_file = sys.argv[4]
    exit_code = run_raxml_backbone(alignment, tree, model, log_file)
    sys.exit(exit_code)