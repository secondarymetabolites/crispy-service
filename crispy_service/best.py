"""Implementation of CRISPR-BEST related features."""

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation

from . import utils


class EditMode(object):
    """CRISPR-BEST edit mode base class."""
    same_strand: dict[str, str] = {}
    opposite_strand: dict[str, str] = {}

    @staticmethod
    def can_edit(seq):
        raise NotImplementedError


class CtoT(EditMode):
    same_strand = {'C': 'T'}
    opposite_strand = {'G': 'A'}

    @staticmethod
    def can_edit(seq):
        return "C" in seq

    def __str__(self):
        return self.name()

    @staticmethod
    def name():
        return "CtoT"


class AtoG(EditMode):
    same_strand = {'A': 'G'}
    opposite_strand = {'T': 'C'}

    @staticmethod
    def can_edit(seq):
        return "A" in seq

    def __str__(self):
        return self.name()

    @staticmethod
    def name():
        return "AtoG"


class FastFeatureLocation(object):
    __slots__ = (
        'start',
        'end',
        'strand'
    )

    def __init__(self, start: int, end: int, strand: int) -> None:
        self.start = start
        self.end = end
        self.strand = strand

    def __iter__(self):
        yield from range(self.start, self.end)

    def __eq__(self, other):
        if not isinstance(other, (FeatureLocation, FastFeatureLocation)):
            return False

        return self.start == other.start and self.end == other.end and self.strand == other.strand


class BestEditWindow(object):
    """CRISPR-BEST edit window."""

    def __init__(self, record, cds, location, mode=CtoT, size=7, offset=13):
        """Set up a CRISPR-BEST edit window.

        Parameters:
            record: The SeqRecord object the CDS and gRNA are located on
            cds: The SeqFeature of the CDS potentially edited by the gRNA
            location: the FeatureLocation of the gRNA responsible for the editing
            mode: The edit mode to use (e.g. CtoT)
            size: The size of the edit window
            offset: Distance from PAM that the edit window starts at
        """
        self.record = record
        self.cds = cds
        self.location = location
        self.mode = mode
        self.size = size
        self.offset = offset

    def get_affected_codons(self):
        """Get all codons of the CDS potentially affected by the gRNA editing.

        Returns:
            A list of CodonChange objects containing all codons covered by the gRNA edit window
        """
        window_location = self.edit_window(self.location, self.size, self.offset)
        for before_codon in self.extract_codons():
            after_codon = before_codon.mutate(window_location, self.mode)
            if not after_codon:
                continue
            yield CodonChange(before_codon, after_codon)

    def get_mutations(self):
        """Get alll codons of the CDS mutated by the gRNA editing.

        Returns:
            A list of CodonChange objects containg mutations due to the editing
        """
        return filter(lambda x: not x.conservative, self.get_affected_codons())

    def extract_codons(self):
        """Extract codons from the CDS."""

        seq = self.cds.extract(self.record.seq)
        strand = self.cds.location.strand
        cds_start = self.cds.location.start
        cds_end = self.cds.location.end

        for idx, i in enumerate(range(0, len(seq), 3)):
            if strand == -1:
                start = cds_end - (i + 3)
                end = cds_end - i
            else:
                start = cds_start + i
                end = cds_start + i + 3
            loc = FastFeatureLocation(start, end, strand)
            yield Codon(seq[i:i + 3], loc, idx)

    @staticmethod
    def can_edit(record, location, mode=CtoT, size=7, offset=13):
        """Check if the edit window can edit, i.e. it contains a C base."""
        window = BestEditWindow.edit_window(location, size, offset)
        seq = window.extract(record.seq)

        return mode.can_edit(seq)

    @staticmethod
    def overlaps(cds, location, size=7, offset=13):
        """Check if the gRNA edit window overlaps with the CDS."""
        window = BestEditWindow.edit_window(location, size, offset)
        return (cds.location.start <= window.start <= cds.location.end) \
            or (cds.location.start <= window.end <= cds.location.end)

    @staticmethod
    def edit_window(location, size=7, offset=13):
        """Get the gRNA edit window coordinates, based on gRNA strand.

            Parameters:
                start: the location of the gRNA to get the edit window for

            Returns:
                A FeatureLocation representing the edit window

        """
        if location.strand == -1:
            return FeatureLocation(location.start + offset, location.start + (offset + size), -1)
        else:
            return FeatureLocation(location.end - (offset + size), location.end - offset, 1)


class Codon(object):
    def __init__(self, seq, location, position):
        self.seq = seq
        self.location = location
        self.position = position

    def translate(self):
        return self.seq.translate()

    def mutate(self, location, mode=CtoT):
        """Mutate this codon in a given location.

            Parameters:
                location: A FeatureLocation object representing the edit window

            Returns:
                A new Codon object representing the mutated codon, or None if location and Codon don't overlap
        """
        if self.location.strand == location.strand:
            base_changes = mode.same_strand
        else:
            base_changes = mode.opposite_strand

        if not utils.locations_overlap(self.location, location):
            return None

        m_start = max(0, location.start - self.location.start)
        m_end = min(3, location.end - self.location.start)

        new_seq_str = ""
        for i, base in enumerate(self.seq):
            if m_start <= i < m_end:
                new_seq_str += base_changes.get(base, base)
            else:
                new_seq_str += base

        return Codon(
            Seq(new_seq_str),
            FastFeatureLocation(int(self.location.start), int(self.location.end), self.location.strand),
            self.position,)

    def __eq__(self, other):
        return self.seq == other.seq and \
            self.location == other.location and \
            self.position == other.position

    def __str__(self):
        return "{s.seq}{{{s.location}}}({aa}{s.position})".format(s=self, aa=self.translate())

    def __repr__(self):
        return str(self)


class CodonChange(object):
    def __init__(self, before, after):
        """Record codon changes caused by CRISPR-BEST editing

            Parameters:
                before: Codon object representing the codon before the mutation
                after: Codon object representing the codon after the mutation

        """

        if before.position != after.position:
            raise ValueError("before.position ({}) != after.position ({})".format(before, after))

        self.before = before
        self.after = after

    @property
    def before_aa(self):
        return self.before.translate()

    @property
    def after_aa(self):
        return self.after.translate()

    @property
    def aa_position(self):
        """Human-readable (1-based) position of the changed amino acid."""
        return self.before.position + 1

    @property
    def conservative(self):
        return self.before_aa == self.after_aa

    def __str__(self):
        return "{s.before_aa}{s.aa_position}{s.after_aa}".format(s=self)

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return self.before == other.before and self.after == other.after
