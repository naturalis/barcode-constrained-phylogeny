#!/usr/bin/env python3

import sys
import re
import os
import gzip
import time
import multiprocessing
from tqdm import tqdm
from functools import partial
from concurrent.futures import ProcessPoolExecutor, as_completed

def read_process_ids_from_phylip(phylip_file):
    """Extract process IDs from a Phylip file."""
    process_ids = set()
    with open(phylip_file, 'r') as f:
        # Skip the first line (contains dimensions)
        next(f)
        for line in f:
            # In Phylip format, the sequence name comes before the first space
            match = re.match(r'^(\S+)', line)
            if match:
                process_id = match.group(1)
                process_ids.add(process_id)
    return process_ids

def identify_process_id_column(header):
    """Identify the process ID column in the BCDM file header."""
    columns = header.strip().split('\t')
    
    # Check for common process ID column names
    for i, col in enumerate(columns):
        col_lower = col.lower()
        if 'processid' in col_lower or 'process_id' in col_lower:
            return i
    
    # If not found, check for other possible ID columns
    for i, col in enumerate(columns):
        col_lower = col.lower()
        if col_lower == 'id' or 'sequenceid' in col_lower:
            return i
    
    # Return None if no suitable column found
    return None

def get_file_line_count(file_path):
    """Count the number of lines in a file efficiently."""
    is_gzipped = file_path.endswith('.gz')
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    print(f"Counting lines in {file_path}...")
    line_count = 0
    chunk_size = 1024 * 1024  # 1MB chunks
    
    with open_func(file_path, mode) as file:
        # Skip header
        next(file)
        line_count = 1  # Start with 1 for the header
        
        with tqdm(unit='MB', desc="Counting lines") as pbar:
            while True:
                chunk = file.read(chunk_size)
                if not chunk:
                    break
                line_count += chunk.count('\n')
                pbar.update(1)
    
    return line_count

def split_file_into_chunks(file_path, num_chunks):
    """Return file offsets for chunks."""
    is_gzipped = file_path.endswith('.gz')
    if is_gzipped:
        print("Warning: Multithreaded processing is not as efficient with gzipped files")
        # For gzipped files, we'll return just one chunk for simplicity
        return [(0, -1)]
    
    total_size = os.path.getsize(file_path)
    chunk_size = total_size // num_chunks
    
    offsets = []
    with open(file_path, 'rb') as f:
        # Read header
        header = f.readline()
        header_offset = len(header)
        
        # First chunk starts after header
        start_offset = header_offset
        
        for i in range(num_chunks - 1):
            # Jump to approximate chunk boundary
            f.seek(start_offset + chunk_size)
            
            # Find the next newline
            while f.read(1) != b'\n':
                pass
            
            # Record end offset for current chunk
            end_offset = f.tell()
            offsets.append((start_offset, end_offset))
            
            # Set start of next chunk
            start_offset = end_offset
        
        # Last chunk goes to end of file
        offsets.append((start_offset, -1))
        
    return offsets

def process_chunk(file_path, offset_range, process_ids, process_id_col, results_queue=None):
    """Process a chunk of the BCDM file."""
    start_offset, end_offset = offset_range
    is_gzipped = file_path.endswith('.gz')
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    matching_records = []
    count = 0
    
    with open_func(file_path, mode) as f:
        # If not the first chunk, get header first
        if start_offset > 0:
            header = f.readline()  # Read and discard header
            f.seek(start_offset)
        
        # Read until end_offset or EOF
        while True:
            if end_offset > 0 and f.tell() >= end_offset:
                break
                
            line = f.readline()
            if not line:
                break
                
            count += 1
            
            fields = line.strip().split('\t')
            if len(fields) > process_id_col:
                pid = fields[process_id_col]
                if pid in process_ids:
                    matching_records.append(line)
    
    if results_queue:
        results_queue.put((matching_records, count))
    return matching_records, count

def filter_bcdm_file(bcdm_file, process_ids, output_file, threads=None):
    """Filter the BCDM file to only include rows with matching process IDs."""
    if threads is None:
        threads = max(1, multiprocessing.cpu_count() - 1)
    
    # Determine if the file is gzipped
    is_gzipped = bcdm_file.endswith('.gz')
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    start_time = time.time()
    print(f"Starting to filter {bcdm_file} using {threads} threads")
    
    # Read the header to identify the process ID column
    with open_func(bcdm_file, mode) as f:
        header = f.readline()
    
    process_id_col = identify_process_id_column(header)
    if process_id_col is None:
        print("Error: Could not identify process ID column in BCDM file")
        return False
    
    print(f"Using column {process_id_col} for process IDs")
    
    # Get approximate number of records for progress tracking
    total_lines = get_file_line_count(bcdm_file)
    print(f"File has approximately {total_lines:,} records (including header)")
    
    # Split file into chunks
    chunks = split_file_into_chunks(bcdm_file, threads)
    print(f"Split file into {len(chunks)} processing chunks")
    
    # Process chunks in parallel
    matching_records = []
    total_processed = 0
    
    # Setup progress bar
    pbar = tqdm(total=total_lines-1, desc="Filtering", unit="records")
    
    # Create a manager and queue for progress updates
    manager = multiprocessing.Manager()
    results_queue = manager.Queue()
    
    # Process chunks in parallel
    chunk_processors = []
    with ProcessPoolExecutor(max_workers=threads) as executor:
        for chunk in chunks:
            processor = executor.submit(
                process_chunk, bcdm_file, chunk, process_ids, process_id_col, results_queue
            )
            chunk_processors.append(processor)
        
        # Monitor the queue and update progress
        all_matches = []
        all_counts = 0
        completed = 0
        
        while completed < len(chunks):
            if not results_queue.empty():
                matches, count = results_queue.get()
                all_matches.extend(matches)
                all_counts += count
                pbar.update(count)
                completed += 1
            else:
                time.sleep(0.1)
    
    pbar.close()
    
    # Write results
    with open(output_file, 'w') as outfile:
        outfile.write(header)
        for record in all_matches:
            outfile.write(record)
    
    elapsed = time.time() - start_time
    print(f"Completed filtering in {elapsed:.2f} seconds")
    print(f"Processed {all_counts:,} records, found {len(all_matches):,} matches")
    
    return True

def main():
    if len(sys.argv) < 4:
        print("Usage: python filter_bcdm.py <phylip_file> <bcdm_file> <output_file> [threads]")
        sys.exit(1)
    
    phylip_file = sys.argv[1]
    bcdm_file = sys.argv[2]
    output_file = sys.argv[3]
    threads = int(sys.argv[4]) if len(sys.argv) > 4 else None
    
    print(f"Reading process IDs from {phylip_file}...")
    process_ids = read_process_ids_from_phylip(phylip_file)
    print(f"Found {len(process_ids):,} unique process IDs")
    
    success = filter_bcdm_file(bcdm_file, process_ids, output_file, threads)
    
    if success:
        print(f"Filtered BCDM saved to {output_file}")
    else:
        print("Failed to filter BCDM file")
        sys.exit(1)

if __name__ == "__main__":
    main()