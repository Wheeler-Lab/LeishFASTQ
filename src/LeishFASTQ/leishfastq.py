import gzip
from collections import defaultdict
import re
import regex
import itertools
import concurrent.futures

DEFAULT_BARCODE_LENGTH = 17
DEFAULT_LEFT_FLANKING_SEQUENCE = "GTGTATCGGATGTCAGTTGC"
DEFAULT_RIGHT_FLANKING_SEQUENCE = "GTATAATGCAGACCTGCTGC"

_complement_map = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G',
}
def _reverse_complement(sequence):
    return ''.join([_complement_map[c] for c in sequence][::-1])

def find_flanked_sequences_in_fastq_file(filename,
                                         flank_sequences=(DEFAULT_LEFT_FLANKING_SEQUENCE, DEFAULT_RIGHT_FLANKING_SEQUENCE),
                                         sequence_length=DEFAULT_BARCODE_LENGTH,
                                         allow_substitutions=0):
    """
    find_flanked_sequences_in_fastq_file(
        filename,
        flank_sequences=('GTGTATCGGATGTCAGTTGC', 'GTATAATGCAGACCTGCTGC'),
        sequence_length=17,
        allow_single_substitution=False)
    
    Find all sequences that are bracketed by the given left and right flank sequences
    (defaults to the LeishGEdit Barcode flanks) in the gzipped fastq file in `filename`.
    The sequences have to have to be `sequence_length` nucleotides long (default 17).
    Optionally allow single subsitutions in the flanking sequences.
    
    Returns a `dict` with the found sequences as keys and their read counts as values.
    """
    flank_left, flank_right = flank_sequences

    # We need to reverse complement the flanking sequences for the reverse reads.
    reverse_read = '_R2_' in filename.stem
    if reverse_read:
        flank_left, flank_right = _reverse_complement(flank_right), _reverse_complement(flank_left)

    # Construct regular expressions
    if allow_substitutions > 0:
        # Use the slightly slower `regex` package if substitutions are allowed.
        # This packages allows fuzzy matching.
        sequence_regex = regex.compile(f'[ATGC]*({flank_left}){{s<={allow_substitutions}}}([ATGC]{{{sequence_length}}})({flank_right}){{s<={allow_substitutions}}}[ATGC]*')
    else:
        # If not, use the builtin regex module.
        sequence_regex = re.compile(f'[ATGC]*({flank_left})([ATGC]{{{sequence_length}}})({flank_right})[ATGC]*')
    
    # Iterate over all 4-line blocks in the file (no checking for corrupt files is done!)
    # Match each sequence in the file with the regex, and record any found flanked sequence.
    found_sequence_count = defaultdict(lambda: 0)
    with gzip.open(filename, 'rt') as f:
        sequences = itertools.islice(f, 1, None, 4)
        for sequence in sequences:
            # Strip the newline off the sequence.
            m = sequence_regex.match(sequence.rstrip('\n'))
            if m is not None:
                # The found sequence is the second group.
                found_sequence_count[m.group(2)] += 1
            else:
                # If no match is found, count this in the "garbage" group
                found_sequence_count[''] += 1
    
    # The barcode sequences need to be reverse complemented for the reverse reads.
    if reverse_read:
        found_sequence_count = {_reverse_complement(seq): count for seq, count in found_sequence_count.items()}

    # Return a regular dict.
    return dict(found_sequence_count)

def match_sequences(known_sequences, candidate_sequence_counts, allow_substitutions):
    """
    match_sequence(known_sequences, candidate_sequence_counts, allow_single_substitution=False)
    
    Match the candidate sequences in `candidate_sequence_counts` against he `known_sequences` list.
    Sequences that are not in `known_sequences` are discarded. If single substitutions are allowed,
    a fuzzy match is performed, and the counts for any matches found during this are added to any
    regular matches.
    
    Returns a `dict` with the matched sequences as keys and their read counts as values.
    """
    mapped = defaultdict(lambda: {'count': 0, 'mismatched': 0})
    for sequence, count in candidate_sequence_counts.items():
        if sequence not in known_sequences and sequence != '':
            sequence_regex = regex.compile(f'({sequence}){{s<={allow_substitutions}}}')
            fuzzy_matches = [candidate for candidate in known_sequences if sequence_regex.match(candidate) is not None]
            if len(fuzzy_matches) == 1:
                mapped[fuzzy_matches[0]]['count'] += count
                mapped[fuzzy_matches[0]]['mismatched'] += count
            elif len(fuzzy_matches) > 1:
                raise ValueError('Multiple matches found - this should not be possible.')
            else:
                mapped[sequence]['count'] += count
        else:
            mapped[sequence]['count'] += count        
    
    return {sequence: mapped[sequence] for sequence in sorted(mapped, key=lambda k: -mapped[k]['count'])}

def process_fastq(filename, known_sequences, **kwargs):
    """
    process_fastq(filename, known_sequences, allow_single_substitution=False)
    
    Find flanked sequences in the gzipped FASTQ file given in `filename` and match them against
    the list `known_sequences`. Optionally allow a single nucleotide substitution in both the flank
    sequences and the matching sequences.
    
    Returns a tuple containing the filename and the dict with the matched sequence read counts.
    """
    found_sequences = find_flanked_sequences_in_fastq_file(filename, **kwargs)
    if (known_sequences is not None) and (kwargs["allow_substitutions"] > 0):
        if (kwargs['allow_substitutions'] > 1):
            print("WARNING: Allowing more than one substitution in the barcode is not allowed because it leads to inconsistencies. Automatically reduced this to one.")
        found_sequences = match_sequences(known_sequences, found_sequences, 1)
    else:
        found_sequences = {seq: {'count': count, 'mismatched': 0} for seq, count in found_sequences.items()}
    return filename.stem, found_sequences

def process_fastq_files(filenames, known_sequences=None, progress_wrap=None, nr_workers=None, **kwargs):
    """
    process_fastq_files(filenames, known_sequences, allow_single_substitution=False, progress_wrap=None)
    
    Batch process the gzipped FASTQ files in `filenames` to find flanked sequences and match them
    against the list `known_sequences`. Optionally allow a single nucleotide substitution in both the flank
    sequences and the matching sequences.
    
    Uses `concurrent.futures.ProcessPoolExecutor` to allow parallel processing. It is possible to specify a
    progress callback such as `tqdm` using `progress_wrap` (default `None`).
    
    Returns a `dict` with the filenames as keys and their sequence read counts as values.
    """
    if progress_wrap is None:
        progress_wrap = lambda x, *a, **kw: x
    with concurrent.futures.ProcessPoolExecutor(max_workers=nr_workers) as executor:
        futures = [executor.submit(process_fastq, f, known_sequences, **kwargs)
                   for f in filenames]
        results = [r.result() for r in progress_wrap(concurrent.futures.as_completed(futures), total=len(futures))]
    return {r[0]: r[1] for r in results}

def command_line():
    import os
    import pathlib
    import argparse

    import pandas as pd
    from tqdm.auto import tqdm

    parser = argparse.ArgumentParser(description="Find barcodes in gzipped fastq files.")
    parser.add_argument('fastq-gz-file', type=pathlib.Path, nargs='+', help="List of FASTQ files to process")
    parser.add_argument('--allow-mismatches', dest="mismatches", type=int, default=0, help="Allow up to the given number of mismatches in each flank and the barcode, if a list of known barcodes is given (default %(default)s)")
    parser.add_argument('--known-barcodes', dest='barcode_list', type=pathlib.Path, default=None, help="File containing a list of known barcodes to match against. This helps if mismatches are allowed.")
    parser.add_argument('--left-flanking-sequence', dest='left_flanking_sequence', type=str, default=DEFAULT_LEFT_FLANKING_SEQUENCE, help="Left barcode flanking sequence (default: %(default)s)")
    parser.add_argument('--right-flanking-sequence', dest='right_flanking_sequence', type=str, default=DEFAULT_RIGHT_FLANKING_SEQUENCE, help="Right barcode flanking sequence (default: %(default)s)")
    parser.add_argument('--barcode-length', dest='barcode_length', type=int, default=17, help="Barcode length in nucleotides (default: %(default)d)")
    parser.add_argument('--nr-workers', dest='nr_workers', type=int, default=None, help="Number of parallel processes to use (default: same as the number of cores)")
    args = parser.parse_args()

    if args.barcode_list is not None:
        with open(args.barcode_list, 'r') as f:
            barcode_list = [line.split(',')[0].strip() for line in f]
    else:
        barcode_list = None

    results = process_fastq_files(
        vars(args)['fastq-gz-file'],
        allow_substitutions=args.mismatches,
        known_sequences=barcode_list,
        nr_workers=args.nr_workers,
        flank_sequences=(args.left_flanking_sequence, args.right_flanking_sequence),
        sequence_length=args.barcode_length,
        progress_wrap=lambda *a, **kw: tqdm(*a, desc="Finding barcodes", **kw)
    )

    for filename, counts in tqdm(results.items(), total=len(results), desc="Saving data"):
        if len(counts) > 0:
            pd.DataFrame([{'barcode_sequence': seq, 'count': v['count'], 'mismatched': v['mismatched']} for seq, v in counts.items()]).sort_values('count', ascending=False).to_csv(filename.replace('.fastq.gz', '') + '.csv', index=False)