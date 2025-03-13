#!/usr/bin/env python3

import sys
import os
import time
import concurrent.futures
import multiprocessing
from Bio import SeqIO
import dendropy
from tqdm import tqdm

def process_alignment_chunk(chunk_file):
    """Process a chunk of the alignment file and return the set of taxa IDs."""
    taxa = set()
    try:
        with open(chunk_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                taxa.add(record.id)
    except Exception as e:
        print(f"Error processing chunk {chunk_file}: {str(e)}")
    return taxa

def process_tree_nodes(nodes):
    """Process a subset of tree nodes and return the set of taxa labels."""
    return {node.taxon.label for node in nodes if node.taxon is not None}

def split_fasta_file(input_file, num_chunks):
    """Split a FASTA file into approximately equal chunks."""
    # Count sequences and determine chunk size
    seq_count = 0
    with open(input_file, "r") as f:
        for line in f:
            if line.startswith('>'):
                seq_count += 1
    
    if seq_count == 0:
        return []
    
    # Adjust number of chunks if there are fewer sequences than requested chunks
    num_chunks = min(num_chunks, seq_count)
    seqs_per_chunk = max(1, seq_count // num_chunks)
    
    # Create chunks
    chunk_files = []
    current_chunk = []
    current_count = 0
    current_chunk_idx = 0
    
    with open(input_file, "r") as f:
        current_seq = []
        for line in f:
            if line.startswith('>'):
                if current_seq:
                    current_chunk.extend(current_seq)
                    current_seq = []
                    current_count += 1
                    
                    # Start a new chunk if needed
                    if current_count >= seqs_per_chunk and current_chunk_idx < num_chunks - 1:
                        chunk_file = f"{input_file}.chunk_{current_chunk_idx}"
                        with open(chunk_file, "w") as chunk_out:
                            chunk_out.writelines(current_chunk)
                        chunk_files.append(chunk_file)
                        current_chunk = []
                        current_count = 0
                        current_chunk_idx += 1
                
            current_seq.append(line)
        
        # Add the last sequence if any
        if current_seq:
            current_chunk.extend(current_seq)
    
    # Write the final chunk
    if current_chunk:
        chunk_file = f"{input_file}.chunk_{current_chunk_idx}"
        with open(chunk_file, "w") as chunk_out:
            chunk_out.writelines(current_chunk)
        chunk_files.append(chunk_file)
    
    return chunk_files

def check_tree_alignment_compatibility(alignment_file, tree_file, log_file=None, threads=None):
    """
    Check if the taxa in the constraint tree match those in the alignment file using multithreading.
    Returns True if they match, False otherwise.
    """
    start_time = time.time()
    
    # Determine number of threads
    if threads is None:
        threads = multiprocessing.cpu_count()
    threads = max(1, min(threads, multiprocessing.cpu_count()))
    
    # Prepare logging
    if log_file:
        log = open(log_file, "w")
        write_log = lambda msg: log.write(f"{msg}\n")
    else:
        write_log = lambda msg: print(msg)
    
    write_log(f"Checking compatibility between tree and alignment...")
    write_log(f"Alignment file: {alignment_file}")
    write_log(f"Tree file: {tree_file}")
    write_log(f"Using {threads} threads")
    
    # Process the alignment file in parallel
    try:
        write_log(f"Splitting alignment file into {threads} chunks for parallel processing...")
        chunk_files = split_fasta_file(alignment_file, threads)
        write_log(f"Created {len(chunk_files)} chunk files")
        
        # Process chunks in parallel
        alignment_taxa = set()
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            future_to_chunk = {executor.submit(process_alignment_chunk, chunk): chunk for chunk in chunk_files}
            for future in tqdm(concurrent.futures.as_completed(future_to_chunk), total=len(chunk_files), desc="Processing alignment chunks"):
                chunk = future_to_chunk[future]
                try:
                    chunk_taxa = future.result()
                    alignment_taxa.update(chunk_taxa)
                except Exception as e:
                    write_log(f"Error processing chunk {chunk}: {str(e)}")
        
        # Clean up chunk files
        for chunk_file in chunk_files:
            try:
                os.remove(chunk_file)
            except:
                pass
                
        write_log(f"Found {len(alignment_taxa)} unique sequence IDs in alignment file.")
    except Exception as e:
        write_log(f"Error reading alignment file: {str(e)}")
        if log_file:
            log.close()
        return False
    
    # Read taxa from the tree file
    try:
        write_log("Parsing tree file...")
        tree = dendropy.Tree.get(path=tree_file, schema="newick")
        
        # Extract taxa with parallel processing
        leaf_nodes = list(tree.leaf_nodes())
        write_log(f"Found {len(leaf_nodes)} leaf nodes in tree")
        
        # Divide leaf nodes into chunks
        chunks = []
        chunk_size = max(1, len(leaf_nodes) // threads)
        for i in range(0, len(leaf_nodes), chunk_size):
            chunks.append(leaf_nodes[i:i+chunk_size])
        
        # Process chunks in parallel
        tree_taxa = set()
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            future_to_chunk = {executor.submit(process_tree_nodes, chunk): chunk for chunk in chunks}
            for future in tqdm(concurrent.futures.as_completed(future_to_chunk), total=len(chunks), desc="Processing tree nodes"):
                try:
                    chunk_taxa = future.result()
                    tree_taxa.update(chunk_taxa)
                except Exception as e:
                    write_log(f"Error processing tree nodes: {str(e)}")
            
        write_log(f"Found {len(tree_taxa)} unique taxa in tree file.")
    except Exception as e:
        write_log(f"Error reading tree file: {str(e)}")
        if log_file:
            log.close()
        return False
    
    # Check for mismatches
    write_log("Comparing taxa sets...")
    alignment_only = alignment_taxa - tree_taxa
    tree_only = tree_taxa - alignment_taxa
    common_taxa = alignment_taxa.intersection(tree_taxa)
    
    write_log(f"Summary:")
    write_log(f"  - Total taxa in alignment: {len(alignment_taxa)}")
    write_log(f"  - Total taxa in tree: {len(tree_taxa)}")
    write_log(f"  - Taxa in both: {len(common_taxa)} ({(len(common_taxa)/max(1, len(alignment_taxa.union(tree_taxa))))*100:.1f}%)")
    
    if alignment_only:
        write_log(f"  - Taxa in alignment but not in tree: {len(alignment_only)} ({(len(alignment_only)/len(alignment_taxa))*100:.1f}% of alignment)")
        sample = list(alignment_only)[:10]
        write_log(f"    Sample: {', '.join(sample)}{'...' if len(alignment_only) > 10 else ''}")
        
        # Save full list of alignment-only taxa to a file if log file is specified
        if log_file:
            alignment_only_file = log_file + ".alignment_only_taxa.txt"
            with open(alignment_only_file, "w") as f:
                for taxon in sorted(alignment_only):
                    f.write(f"{taxon}\n")
            write_log(f"    Full list saved to: {alignment_only_file}")
    else:
        write_log(f"  - All alignment taxa are present in the tree")
    
    if tree_only:
        write_log(f"  - Taxa in tree but not in alignment: {len(tree_only)} ({(len(tree_only)/len(tree_taxa))*100:.1f}% of tree)")
        sample = list(tree_only)[:10]
        write_log(f"    Sample: {', '.join(sample)}{'...' if len(tree_only) > 10 else ''}")
        
        # Save full list of tree-only taxa to a file if log file is specified
        if log_file:
            tree_only_file = log_file + ".tree_only_taxa.txt"
            with open(tree_only_file, "w") as f:
                for taxon in sorted(tree_only):
                    f.write(f"{taxon}\n")
            write_log(f"    Full list saved to: {tree_only_file}")
    else:
        write_log(f"  - All tree taxa are present in the alignment")
    
    # Determine compatibility status
    if len(common_taxa) == 0:
        compatibility_status = "NOT compatible at all"
        is_compatible = False
    elif len(alignment_only) == 0 and len(tree_only) == 0:
        compatibility_status = "FULLY compatible"
        is_compatible = True
    else:
        compatibility_status = "NOT FULLY compatible"
        is_compatible = False
    
    write_log(f"\nTree and alignment are {compatibility_status}.")
    
    # Print a summary of what needs to be done
    if not is_compatible:
        write_log("\nRecommendation:")
        if compatibility_status == "NOT compatible at all":
            write_log("  - The tree and alignment have no taxa in common!")
            write_log("  - Check if the taxa naming schemes are different between the files")
            write_log("  - Ensure you're using the correct input files")
        else:  # Partially compatible
            if tree_only:
                write_log("  - Prune the tree to remove taxa not in the alignment")
                write_log("    (Use a tree manipulation tool like 'nw_prune' from Newick Utilities)")
                write_log(f"    Command example: nw_prune {tree_file} `cat {log_file}.tree_only_taxa.txt` > pruned_tree.tre")
            if alignment_only:
                write_log("  - Either add missing taxa to the tree or filter the alignment")
                write_log("    (Use a sequence filtering tool to keep only sequences with IDs in the tree)")
    else:
        write_log("\nNo action needed - tree and alignment have identical taxa.")

    elapsed_time = time.time() - start_time
    write_log(f"Compatibility check completed in {elapsed_time:.2f} seconds.")
    
    if log_file:
        log.close()
    
    return is_compatible

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python check_tree_alignment.py <alignment_file> <tree_file> [log_file] [threads]")
        sys.exit(1)
    
    alignment_file = sys.argv[1]
    tree_file = sys.argv[2]
    log_file = sys.argv[3] if len(sys.argv) > 3 else None
    threads = int(sys.argv[4]) if len(sys.argv) > 4 else None
    
    compatible = check_tree_alignment_compatibility(alignment_file, tree_file, log_file, threads)
    
    # Exit with appropriate code (0 for compatible, 1 for incompatible)
    sys.exit(0 if compatible else 1)