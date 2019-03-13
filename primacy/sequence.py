from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from primacy.primer import ambiguous_count

class Sequence(object):

    alphabet = IUPACAmbiguousDNA()
    def __init__(
            self, seq_id, seq,
            target_start, target_end,
            primer_len_min,
            primer_len_max, **qwargs):
        """
        Takes sequence info and returns all upstream
        and downstream primers in a range of sizes defined
        in the primer_len_range

        seq_id: (str) unique identifier of the sequence
        seq: (str) sequence of target region with upstream and downstream
                flanking regions. Will be converted into Seq object with
                ambiguous bases allowed, X represents
                regions to avoid when getting primers.
        target_start: (int) the index of the start of the target region
        target_end: (int) the index of the last base in the target region
        primer_len_min: (int) the min size of primers that should be generated.
        primer_len_max: (int) the max size of primers that should be generated. 
        The primer_len_min and primer_len_max and target_start and target_end
        are converted into a ranges (primer_len_range and target_range,
        respectively) with the max incremented by one so that it is included. 
        """
        self.seq_id = seq_id
        self.seq = seq
        self.target_range = (target_start, target_end)
        self.primer_len_range = (primer_len_min, primer_len_max)
    
    @property
    def seq_id(self):
        return self.__seq_id

    @seq_id.setter
    def seq_id(self, seq_id):
        self.__seq_id = seq_id

    @property
    def seq(self):
        return self.__seq

    @seq.setter
    def seq(self, seq):
        self.__seq = Seq(seq, Sequence.alphabet).upper()
    
    @property
    def target_range(self):
        return self.__target_range
    
    
    @target_range.setter
    def target_range(self, range_):
        """Tuple representing the target range within the sequence"""
        start, stop = range_
        try:
            assert 0 <= start <= stop < len(self.seq)
            self.__target_range = (start, stop + 1)
        except AssertionError:
            raise("Target sequence out of range")
    
    
    @property
    def primer_len_range(self):
        return self.__primer_len_range

    @primer_len_range.setter
    def primer_len_range(self, range_):
        """Tuple representing the primer length range"""
        min_, max_ = range_
        try:
            assert min_ <= max_
            self.__primer_len_range = (min_, max_ + 1)
        except AssertionError:
            raise("Primer range is not valid")

    @property
    def left_flank(self):
        return self.seq[:self.target_range[0]]

    @property
    def right_flank(self):
        return self.seq.reverse_complement()[
            :len(self.seq) - self.target_range[1]]

    @property
    def primers(self):
        return self._get_primers_by_flank("F") + \
        self._get_primers_by_flank("R")


    def _get_primers_by_flank(self, flank):
        """
        Returns a list of all possible primers in a sequence over range of 
        primer sizes [ SEQ_ID, PRIMER_ID, FLANK, SEQ]. 
        The list contains each primer id representing
        the seq_id, the position in the sequence, the length of the primer and 
        the flank (SEQID_POS_LEN_FLANK). The position is defined such that 
        POS_L + POS_R + len(target) = total amplicon length.
        """
        seq = self.right_flank if flank == "R" else self.left_flank
        primers = []
        for primer_size in range(
            self.primer_len_range[0],
            self.primer_len_range[1]):
            for i in range(0, len(seq) - primer_size + 1):
                primer_seq = seq[i: i + primer_size]
                if "X" not in primer_seq and ambiguous_count(
                    primer_seq) < 1024:
                    primers.append(
                        [self.seq_id,
                        len(seq) - i,
                        primer_size,
                        flank,
                        primer_seq])
        return primers

        

        
 
