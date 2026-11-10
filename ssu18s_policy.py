"""User-defined inclusive input length limits, applied before every analysis."""
MIN_18S_BP=800
MAX_18S_BP=4000


def in_length_range(row):
    return MIN_18S_BP<=len(row['sequence'])<=MAX_18S_BP


def filter_records(records):
    return [r for r in records if in_length_range(r)]


def gc_percent(sequence):
    if not sequence or set(sequence)-set('ACGT'):
        raise ValueError('GC calculation requires nonempty unambiguous DNA')
    return 100*(sequence.count('G')+sequence.count('C'))/len(sequence)
