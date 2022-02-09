import gzip
from collections import defaultdict
import re
import regex
import itertools
import concurrent.futures

def find_flanked_sequences_in_fastq_file(filename,
                                         flank_sequences=('GTGTATCGGATGTCAGTTGC', 'GTATAATGCAGACCTGCTGC'),
                                         sequence_length=17,
                                         allow_single_substitution=False):
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
    # Construct regular expressions
    if allow_single_substitution:
        # Use the slightly slower `regex` package if substitutions are allowed.
        # This packages allows fuzzy matching.
        sequence_regex = regex.compile(f'[ATGC]*({flank_left}){{s<=1}}([ATGC]{{{sequence_length}}})({flank_right}){{s<=1}}[ATGC]*')
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

    # Return a regular dict.
    return dict(found_sequence_count)

def match_sequences(known_sequences, candidate_sequence_counts, allow_single_substitution=False):
    """
    match_sequence(known_sequences, candidate_sequence_counts, allow_single_substitution=False)
    
    Match the candidate sequences in `candidate_sequence_counts` against he `known_sequences` list.
    Sequences that are not in `known_sequences` are discarded. If single substitutions are allowed,
    a fuzzy match is performed, and the counts for any matches found during this are added to any
    regular matches.
    
    Returns a `dict` with the matched sequences as keys and their read counts as values.
    """
    mapped = defaultdict(lambda: 0)
    for sequence, count in candidate_sequence_counts.items():
        if sequence in known_sequences:
            mapped[sequence] += count
        else:
            if allow_single_substitution:
                sequence_regex = regex.compile(f'({sequence}){{s<=1}}')
                fuzzy_matches = [candidate for candidate in known_sequences if sequence_regex.match(candidate) is not None]
                if len(fuzzy_matches) == 1:
                    mapped[fuzzy_matches[0]] += count
                elif len(fuzzy_matches) > 1:
                    raise ValueError('Multiple matches found - this should not be possible.')
    
    return {sequence: mapped[sequence] for sequence in sorted(mapped, key=lambda k: -mapped[k])}

def process_fastq(filename, known_sequences, allow_single_substitution=False):
    """
    process_fastq(filename, known_sequences, allow_single_substitution=False)
    
    Find flanked sequences in the gzipped FASTQ file given in `filename` and match them against
    the list `known_sequences`. Optionally allow a single nucleotide substitution in both the flank
    sequences and the matching sequences.
    
    Returns a tuple containing the filename and the dict with the matched sequence read counts.
    """
    found_sequences = find_flanked_sequences_in_fastq_file(filename, allow_single_substitution=allow_single_substitution)
    return filename.stem, match_sequences(known_sequences, found_sequences, allow_single_substitution=allow_single_substitution)

def process_fastq_files(filenames, known_sequences, allow_single_substitution=False, progress_wrap=None):
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
        progress_wrap = lambda x: x
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_fastq, f, known_sequences, allow_single_substitution=allow_single_substitution)
                   for f in filenames]
        results = [r.result() for r in progress_wrap(concurrent.futures.as_completed(futures), total=len(futures))]
    return {r[0]: r[1] for r in results}