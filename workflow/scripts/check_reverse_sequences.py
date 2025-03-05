#!/usr/bin/env python3

import argparse
import tempfile
import os
import logging
import sys
import concurrent.futures
from functools import partial
from copy import deepcopy
from Bio import SeqIO
from Bio.AlignIO import read as read_alignment
from subprocess import run
from tqdm import tqdm  # Add this for progress bars


def setup_logger(name, level_str):
    """Set up and return a logger with the specified name and level"""
    level = getattr(logging, level_str.upper())
    logger = logging.getLogger(name)
    logger.setLevel(level)
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger

def align_score(record, hmmfile, logger):
    """
    Uses a Hidden Markov Model to align a sequence using hmmalign.
    Returns the score based on posterior probabilities.
    """
    try:
        # Make a copy of the record to avoid modifying the original
        clean_record = deepcopy(record)
        
        # Remove dashes from the sequence
        clean_record.seq = clean_record.seq.replace('-', '')
        
        # Use unique prefix for temp files to ensure thread safety
        prefix = f"seq_{record.id.replace('|', '_')}_"
        
        # Open temporary files for the sequence and alignment
        with tempfile.NamedTemporaryFile(mode='w+', prefix=prefix) as temp_fasta, \
             tempfile.NamedTemporaryFile(mode='w+', prefix=prefix) as temp_stockholm:

            # Save the cleaned sequence to the temporary file
            SeqIO.write(clean_record, temp_fasta.name, 'fasta')

            # Run hmm align, read the aligned sequence
            run(['hmmalign', '--trim', '-o', temp_stockholm.name, hmmfile, temp_fasta.name], check=True)
            alignment = read_alignment(temp_stockholm.name, "stockholm")

        # Get posterior probability string
        quality_string = alignment.column_annotations.get('posterior_probability', '')

        count = 0
        # Count . and * characters
        dot_count = quality_string.count('.')
        star_count = quality_string.count('*')

        # Give value 0 to . and value 10 to *
        count += star_count * 10

        # Add all numbers
        digit_sum = sum(int(char) for char in quality_string if char.isdigit())
        count += digit_sum

        # Calculate average count
        average_count = count / len(quality_string) if quality_string else 0
        return average_count, alignment[0]
    
    except Exception as e:
        logger.error(f"Error in align_score for {record.id}: {e}")
        return 0, None

def process_sequence(record, hmmfile, logger):
    """
    Process a single sequence, testing both orientations.
    Returns the record in the correct orientation and whether it was reversed.
    """
    # Score original orientation
    count_fwd, alignment_fwd = align_score(record, hmmfile, logger)
    logger.debug(f'Forward alignment score {count_fwd} for {record.id}')
    
    # Make a copy for reverse complement
    rev_record = deepcopy(record)
    rev_record.seq = rev_record.seq.reverse_complement()
    
    # Score reverse orientation
    count_rev, alignment_rev = align_score(rev_record, hmmfile, logger)
    logger.debug(f'Reverse alignment score {count_rev} for {record.id}')
    
    # Keep the orientation with higher score
    if count_fwd >= count_rev:
        logger.debug(f'Keeping forward orientation for {record.id}')
        return record, False
    else:
        logger.debug(f'Keeping reverse orientation for {record.id}')
        return rev_record, True

def correct_revcom(hmmfile, sequences, logger, threads=4):
    """
    Check each sequence with hmmalign in both orientations
    and keep the orientation with the higher score.
    Multi-threaded version with progress reporting.
    """
    corrected_seqs = []
    reversed_count = 0
    
    # Create a partial function for processing sequences
    process_func = partial(process_sequence, hmmfile=hmmfile, logger=logger)
    
    # Process sequences in parallel with progress bar
    logger.info(f"Starting parallel processing with {threads} threads")
    
    # Using ThreadPoolExecutor with tqdm for progress tracking
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        # Submit all tasks and get futures
        futures = [executor.submit(process_func, seq) for seq in sequences]
        
        # Process results as they complete with progress bar
        for i, future in enumerate(tqdm(concurrent.futures.as_completed(futures), 
                                        total=len(futures), 
                                        desc="Processing sequences")):
            record, is_reversed = future.result()
            corrected_seqs.append(record)
            if is_reversed:
                reversed_count += 1
            
            # Log progress periodically (e.g., every 1000 sequences)
            if (i + 1) % 1000 == 0 or i == 0:
                logger.info(f"Processed {i + 1}/{len(sequences)} sequences, {reversed_count} reversed so far")
    
    logger.info(f'Corrected {reversed_count} reverse complemented sequences out of {len(sequences)}')
    return corrected_seqs

def main():
    parser = argparse.ArgumentParser(description='Check and correct reverse complemented sequences using HMM')
    parser.add_argument('fasta', help='Input FASTA file')
    parser.add_argument('hmm', help='HMM model file')
    parser.add_argument('-o', '--output', help='Output FASTA file (default: corrected_output.fa)', 
                        default="corrected_output.fa")
    parser.add_argument('-v', '--verbosity', default='INFO', 
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], 
                        help='Log level (default: INFO)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='Number of threads to use (default: 4)')
    parser.add_argument('-l', '--log', help='Log file (optional, if not specified logs to stdout)')
    
    args = parser.parse_args()
    
    # Set up logging to file if specified
    logger = setup_logger('check_reverse', args.verbosity)
    if args.log:
        file_handler = logging.FileHandler(args.log)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
    
    # Read sequences
    sequences = list(SeqIO.parse(args.fasta, "fasta"))
    logger.info(f"Read {len(sequences)} sequences from {args.fasta}")
    
    # Correct sequences with multiple threads
    corrected_sequences = correct_revcom(args.hmm, sequences, logger, args.threads)
    
    # Write corrected sequences
    SeqIO.write(corrected_sequences, args.output, "fasta")
    logger.info(f"Wrote {len(corrected_sequences)} corrected sequences to {args.output}")

if __name__ == "__main__":
    main()