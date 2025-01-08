import argparse
import subprocess
import os

def run_raxml(alignment, tree, output, model, num_outgroups, log_file):
    # Count the number of taxa in the alignment
    with open(alignment, 'r') as f:
        taxon_count = sum(1 for line in f if line.startswith('>'))
    
    # Extract the outgroup names
    with open(alignment, 'r') as f:
        lines = [line.strip() for line in f if line.startswith('>')]
        outgroups = ",".join(line[1:] for line in lines[-num_outgroups:])
    
    # Check the constraint tree's properties
    constraint_tree_exists = os.path.exists(tree) and os.path.getsize(tree) > 0
    constraint_tree_valid = False
    if constraint_tree_exists:
        with open(tree, 'r') as f:
            content = f.read()
            constraint_tree_valid = ('(' in content and ')' in content and
                                     content.count('(') > 2 and
                                     content.count(',') < (taxon_count - 1))

    # Prepare the raxml-ng command
    if constraint_tree_valid:
        cmd = [
            "raxml-ng", "--redo", "--msa", alignment, "--outgroup", outgroups,
            "--model", model, "--tree-constraint", tree, "--search"
        ]
        log_message = "Running RAxML-NG with constraint tree"
    else:
        cmd = [
            "raxml-ng", "--redo", "--msa", alignment, "--model", model,
            "--search"
        ]
        log_message = "Constraint tree fully-resolved, running RAxML-NG without it"

    # Log the process and execute the command
    with open(log_file, 'w') as log:
        log.write(log_message + "\n")
        try:
            subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            log.write(f"RAxML-NG failed: {e}\n")
            raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run RAxML-NG with or without a constraint tree.")
    parser.add_argument("--alignment", required=True, help="Path to the alignment file.")
    parser.add_argument("--tree", required=True, help="Path to the constraint tree file.")
    parser.add_argument("--output", required=True, help="Path to the output tree file.")
    parser.add_argument("--model", required=True, help="Model to use for RAxML-NG.")
    parser.add_argument("--num_outgroups", type=int, required=True, help="Number of outgroups.")
    parser.add_argument("--log_file", required=True, help="Path to the log file.")
    args = parser.parse_args()

    run_raxml(
        alignment=args.alignment,
        tree=args.tree,
        output=args.output,
        model=args.model,
        num_outgroups=args.num_outgroups,
        log_file=args.log_file
    )
