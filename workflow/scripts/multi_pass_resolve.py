#!/usr/bin/env python3
"""
Multi-Pass Polytomy Resolution with Multithreading

Resolves polytomies in multiple passes with optimized parameters using parallel processing.
Includes rate limiting to prevent API throttling.
"""

import subprocess
import argparse
import logging
import os
import re
import tempfile
import sys
import multiprocessing
import time
import signal
import random
from functools import partial

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger('multi_pass_resolve')

def count_polytomies(tree_file):
    """Count polytomies in the tree"""
    # Create a temporary file to store polytomy info
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp:
        temp_file = temp.name
    
    try:
        # Get script path - handles running from any directory
        script_dir = os.path.dirname(os.path.abspath(__file__))
        find_polytomies_script = os.path.join(script_dir, "find_polytomies.py")
        
        # Run find_polytomies.py and redirect output to the temp file
        subprocess.run([
            "python", find_polytomies_script,
            "-t", tree_file
        ], stdout=open(temp_file, 'w'), check=True)
        
        # Parse the output to get the count
        with open(temp_file, 'r') as f:
            first_line = f.readline().strip()
            match = re.search(r'Found (\d+) polytomies', first_line)
            count = int(match.group(1)) if match else 0
            
        return count
    except Exception as e:
        logger.error(f"Error in count_polytomies: {e}")
        return 0
    finally:
        # Clean up
        if os.path.exists(temp_file):
            os.remove(temp_file)

def process_chunk(chunk_data):
    """Process a chunk of polytomies in a separate process"""
    input_tree, output_tree, min_depth, min_size, max_size, chunk_id, chunk_start, chunk_end, delay = chunk_data
    
    try:
        # Get script paths
        script_dir = os.path.dirname(os.path.abspath(__file__))
        find_polytomies_script = os.path.join(script_dir, "find_polytomies.py")
        resolve_polytomies_script = os.path.join(script_dir, "resolve_polytomies.py")
        
        # Create a temporary file for the polytomies
        with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp:
            polytomies_file = temp.name
        
        # Add initial delay with some jitter to avoid all processes starting at once
        if delay > 0:
            jitter = random.uniform(0, delay * 0.5)  # Add up to 50% jitter
            logger.info(f"Thread {chunk_id}: Waiting {delay + jitter:.2f}s before starting")
            time.sleep(delay + jitter)
        
        # Find all polytomies first
        subprocess.run([
            "python", find_polytomies_script,
            "-t", input_tree,
            "--min-depth", str(min_depth),
            "--min-size", str(min_size), 
            "--max-size", str(max_size)
        ], stdout=open(polytomies_file, 'w'), check=True)
        
        # Process only this chunk
        logger.info(f"Thread {chunk_id}: Processing polytomies {chunk_start} to {chunk_end}")
        
        # Add --delay parameter if delay is specified
        cmd = [
            "python", resolve_polytomies_script,
            "-t", input_tree,
            "-p", polytomies_file,
            "-o", output_tree,
            "--min-depth", str(min_depth),
            "--min-size", str(min_size),
            "--max-size", str(max_size),
            "--chunk-start", str(chunk_start),
            "--chunk-end", str(chunk_end)
        ]
        
        # Add delay parameter if needed
        if delay > 0:
            cmd.extend(["--delay", str(delay)])
        
        subprocess.run(cmd, check=True)
        
        # Clean up
        os.remove(polytomies_file)
        
        logger.info(f"Thread {chunk_id}: Completed processing")
        return True
        
    except Exception as e:
        logger.error(f"Thread {chunk_id} error: {e}")
        return False

def run_pass_parallel(input_tree, output_tree, pass_number, min_depth, min_size, max_size, threads, delay=0):
    """Run a single pass of polytomy resolution using parallel processing"""
    # Get script paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    find_polytomies_script = os.path.join(script_dir, "find_polytomies.py")
    
    # Create a temporary file for the polytomies
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp:
        polytomies_file = temp.name
    
    try:
        # Find polytomies for this pass
        subprocess.run([
            "python", find_polytomies_script,
            "-t", input_tree,
            "--min-depth", str(min_depth),
            "--min-size", str(min_size), 
            "--max-size", str(max_size)
        ], stdout=open(polytomies_file, 'w'), check=True)
        
        # Count initial polytomies for this pass
        polytomy_count = 0
        with open(polytomies_file, 'r') as f:
            first_line = f.readline().strip()
            match = re.search(r'Found (\d+) polytomies', first_line)
            polytomy_count = int(match.group(1)) if match else 0
        
        if polytomy_count == 0:
            logger.info(f"Pass {pass_number}: No polytomies found with specified parameters")
            # Copy input tree to output
            subprocess.run(["cp", input_tree, output_tree])
            return 0, 0, 0
        
        # If delay is enabled, apply a more reasonable thread cap
        if delay > 0:
            # Allow more threads even with longer delays
            max_threads = max(1, min(threads, 4))  # Use up to 4 threads with delay
            threads_to_use = min(max_threads, polytomy_count)
            logger.info(f"Using {threads_to_use} threads with {delay}s delay between requests")
        else:
            threads_to_use = min(threads, polytomy_count, multiprocessing.cpu_count())
        
        if threads_to_use <= 1 or polytomy_count <= 5:
            # For small numbers of polytomies, just use single process
            logger.info(f"Processing {polytomy_count} polytomies in a single thread")
            resolve_script = os.path.join(script_dir, "resolve_polytomies.py")
            
            cmd = [
                "python", resolve_script,
                "-t", input_tree,
                "-p", polytomies_file,
                "-o", output_tree,
                "--min-depth", str(min_depth),
                "--min-size", str(min_size),
                "--max-size", str(max_size)
            ]
            
            # Add delay parameter if needed
            if delay > 0:
                cmd.extend(["--delay", str(delay)])
                
            subprocess.run(cmd, check=True)
        else:
            # Use multithreading - divide polytomies into chunks
            logger.info(f"Processing {polytomy_count} polytomies with {threads_to_use} threads")
            
            # Calculate chunk size
            chunk_size = max(1, polytomy_count // threads_to_use)
            
            # Create temp output files for each thread
            temp_outputs = []
            for i in range(threads_to_use):
                with tempfile.NamedTemporaryFile(delete=False, suffix=f'.thread{i}.tre') as temp_out:
                    temp_outputs.append(temp_out.name)
            
            # Prepare arguments for each thread
            chunk_args = []
            for i in range(threads_to_use):
                chunk_start = i * chunk_size
                chunk_end = min((i + 1) * chunk_size, polytomy_count)
                
                # Skip empty chunks
                if chunk_start >= chunk_end:
                    continue
                    
                chunk_args.append((
                    input_tree,
                    temp_outputs[i],
                    min_depth,
                    min_size,
                    max_size,
                    i + 1,  # chunk id for logging
                    chunk_start,
                    chunk_end,
                    delay   # Add delay parameter
                ))
            
            # Process chunks in parallel
            start_time = time.time()
            with multiprocessing.Pool(processes=threads_to_use) as pool:
                results = pool.map(process_chunk, chunk_args)
            
            logger.info(f"Parallel processing completed in {time.time() - start_time:.2f} seconds")
            
            # Merge the results by finding a successful output
            merged = False
            for i, success in enumerate(results):
                if success and os.path.exists(temp_outputs[i]):
                    # Copy this file to the output
                    subprocess.run(["cp", temp_outputs[i], output_tree])
                    merged = True
                    break
            
            # If no successful output, fall back to the input tree
            if not merged:
                logger.warning("No successful thread outputs, using input tree")
                subprocess.run(["cp", input_tree, output_tree])
            
            # Clean up temp files
            for tmp_file in temp_outputs:
                if os.path.exists(tmp_file):
                    os.remove(tmp_file)
        
        # Count remaining polytomies
        final_count = count_polytomies(output_tree)
        
        # Calculate how many were resolved
        resolved = polytomy_count - final_count
        
        return polytomy_count, resolved, final_count
    
    except Exception as e:
        logger.error(f"Error in run_pass_parallel: {e}")
        # Ensure we have an output file
        subprocess.run(["cp", input_tree, output_tree])
        return 0, 0, 0
    finally:
        # Clean up
        if os.path.exists(polytomies_file):
            os.remove(polytomies_file)

def multi_pass_resolve(input_tree, output_tree, threads=1, save_intermediate=False, keep_temp=False, delay=0):
    """Run multiple passes of polytomy resolution"""
    # Define passes: most reliable to least reliable
    passes = [
        {"name": "Deep genus/family groups", "min_depth": 5, "min_size": 3, "max_size": 40},
        {"name": "Mid-level groups", "min_depth": 3, "min_size": 3, "max_size": 60},
        {"name": "Shallow order/class groups", "min_depth": 1, "min_size": 3, "max_size": 80},
        {"name": "All remaining groups", "min_depth": 0, "min_size": 2, "max_size": 150}
    ]
    
    # Count initial polytomies
    total_initial = count_polytomies(input_tree)
    logger.info(f"Starting with {total_initial} polytomies in the tree")
    logger.info(f"Using {threads} threads for parallel processing")
    if delay > 0:
        logger.info(f"Rate limiting: {delay} seconds delay between API requests")
    
    # Create a temporary file for intermediate results
    current_tree = input_tree
    intermediate_trees = []
    temp_files = []  # Track all temporary files
    
    # Run each pass
    total_resolved = 0
    
    for i, pass_config in enumerate(passes):
        pass_num = i + 1
        pass_name = pass_config["name"]
        min_depth = pass_config["min_depth"]
        min_size = pass_config["min_size"]
        max_size = pass_config["max_size"]
        
        logger.info(f"\n--- PASS {pass_num}: {pass_name} ---")
        logger.info(f"Parameters: min_depth={min_depth}, min_size={min_size}, max_size={max_size}")
        
        # Create output file for this pass
        if save_intermediate:
            # Create a named intermediate file
            output_base = os.path.splitext(output_tree)[0]
            temp_output = f"{output_base}_pass{pass_num}.tre"
            intermediate_trees.append(temp_output)
        else:
            # Use a temporary file
            with tempfile.NamedTemporaryFile(delete=False, suffix=f'.pass{pass_num}.tre') as temp:
                temp_output = temp.name
                temp_files.append(temp_output)
                logger.info(f"Pass {pass_num} output tree: {temp_output}")
        
        # Run this pass with parallel processing
        initial, resolved, remaining = run_pass_parallel(
            current_tree, temp_output, pass_num, min_depth, min_size, max_size, threads, delay
        )
        
        # Update total resolved
        total_resolved += resolved
        
        logger.info(f"Pass {pass_num} results:")
        logger.info(f"  - Targeted polytomies: {initial}")
        logger.info(f"  - Resolved in this pass: {resolved}")
        logger.info(f"  - Resolution rate for pass: {(resolved/initial*100) if initial > 0 else 0:.1f}%")
        logger.info(f"  - Total resolved so far: {total_resolved}")
        logger.info(f"  - Overall resolution: {(total_resolved/total_initial*100) if total_initial > 0 else 0:.1f}%")
        
        if save_intermediate:
            logger.info(f"  - Saved intermediate tree to: {temp_output}")
        
        # Update current tree for next pass
        current_tree = temp_output
    
    # Final pass - ensure tree is fully bifurcating with format_constraint_tree.py
    logger.info("\n--- FINAL PASS: Force Bifurcation ---")
    
    # Get script path
    script_dir = os.path.dirname(os.path.abspath(__file__))
    format_script = os.path.join(script_dir, "format_constraint_tree.py")
    
    subprocess.run([
        "python", format_script,
        "-i", current_tree,
        "-o", output_tree,
        "-m", "bifurcate"
    ], check=True)
    
    # Count final polytomies to verify
    final_polytomies = count_polytomies(output_tree)
    logger.info(f"\n===== FINAL RESULTS =====")
    logger.info(f"Initial polytomies: {total_initial}")
    logger.info(f"Resolved by OpenToL: {total_resolved}")
    logger.info(f"Remaining after forced bifurcation: {final_polytomies}")
    logger.info(f"Resolution rate: {(1 - final_polytomies/total_initial)*100 if total_initial > 0 else 100:.1f}%")
    
    # List all temporary files if requested
    if keep_temp and temp_files:
        logger.info("\n===== TEMPORARY FILES =====")
        for i, file in enumerate(temp_files):
            if os.path.exists(file):
                logger.info(f"{i+1}. {file}")
        logger.info("==========================\n")
    
    # Clean up temp files if not keeping them
    if not keep_temp:
        for file in temp_files:
            if os.path.exists(file) and file != input_tree:
                os.remove(file)
    else:
        logger.info("Temporary files have been kept")

def main():
    parser = argparse.ArgumentParser(description="Multi-pass polytomy resolution")
    parser.add_argument("-t", "--tree", required=True, help="Input tree file")
    parser.add_argument("-o", "--output", required=True, help="Output tree file")
    parser.add_argument("-j", "--threads", type=int, default=multiprocessing.cpu_count(),
                        help=f"Number of threads to use (default: {multiprocessing.cpu_count()})")
    parser.add_argument("-i", "--intermediate", action="store_true", 
                        help="Save intermediate trees")
    parser.add_argument("--keep-temp", action="store_true", 
                        help="Keep temporary files (don't delete them)")
    parser.add_argument("--slow", type=float, default=0,
                        help="Add delay between API requests in seconds (default: 0, no delay)")
    parser.add_argument("-v", "--verbose", action="store_true", 
                        help="Enable verbose output")
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Create log file
    log_file = f"multi_pass_resolve_{time.strftime('%Y%m%d_%H%M%S')}.log"
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(file_handler)
    
    logger.info(f"Starting multi-pass resolution with {args.threads} threads")
    if args.slow > 0:
        logger.info(f"Rate limiting enabled: {args.slow}s delay between API requests")
    logger.info(f"Log file: {log_file}")
    
    start_time = time.time()
    
    try:
        multi_pass_resolve(args.tree, args.output, args.threads, args.intermediate, args.keep_temp, args.slow)
        elapsed = time.time() - start_time
        logger.info(f"Total execution time: {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")
    except Exception as e:
        logger.error(f"Critical error: {e}")
        import traceback
        logger.error(traceback.format_exc())

if __name__ == "__main__":
    main()