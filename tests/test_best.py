from __future__ import print_function, division

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
import pytest

from crispy_service.best import (
    AtoG,
    BestEditWindow,
    Codon,
    CodonChange,
)
from crispy_service.crispy import crispy_scan


@pytest.fixture
def rec():
    seq = Seq("ATTTA"*10 +
              "ATGGCGCGATACATGGCCGGCATCTGCCACGCCACCTGA" +
              "ATTTA"*10 +
              "TCAGAACAGGGCGTCGTTGGCTCCGTCGTTCGCGTAGCGTCGCTCCAT" +
              "ATTTA"*10)
    mary = SeqFeature(FeatureLocation(50, 89, 1), type="CDS", id="mary")
    merry = SeqFeature(FeatureLocation(139, 187, -1), type="CDS", id="merry")
    rec = SeqRecord(seq, id="FAKE", features=[mary, merry])

    crispy_scan([rec], rec)
    rec.features.sort(key=lambda x: x.location.start)
    rec.features.sort(key=lambda x: int(x.qualifiers.get('1bpmm', ['0'])[0]))

    return rec


@pytest.fixture
def cdses(rec):
    return [f for f in rec.features if f.type == 'CDS']


@pytest.fixture
def grnas(rec):
    return [FeatureLocation(i["start"], i["end"], i["strand"]) for i in crispy_scan([rec], rec)]


def make_codon(seq, start, end, strand, position):
    return Codon(Seq(seq), FeatureLocation(start, end, strand, position), position)


def make_codon_change(before, after, position, strand=1):
    start = position
    end = position + 3

    before_codon = make_codon(before, start, end, strand, position // 3)
    after_codon = make_codon(after, start, end, strand, position // 3)
    return CodonChange(before_codon, after_codon)


def test_extract_codons_fw(rec, cdses, grnas):
    cds = cdses[0]
    grna = grnas[0]  # any will do, we're not really using these

    expected = []
    for idx, triplet in enumerate(["ATG", "GCG", "CGA", "TAC", "ATG", "GCC", "GGC", "ATC", "TGC",
                                   "CAC", "GCC", "ACC", "TGA"]):
        codon = Codon(Seq(triplet),
                      FeatureLocation(cds.location.start + (idx * 3),
                                      cds.location.start + (idx * 3) + 3, 1),
                      idx)
        expected.append(codon)

    window = BestEditWindow(rec, cds, grna)
    ret = list(window.extract_codons())
    assert ret == expected


def test_extract_codons_rv(rec, cdses, grnas):
    cds = cdses[1]
    grna = grnas[0]  # any will do, we're not really using these

    expected = []
    for idx, triplet in enumerate(["ATG", "GAG", "CGA", "CGC", "TAC", "GCG", "AAC", "GAC", "GGA", "GCC",
                                   "AAC", "GAC", "GCC", "CTG", "TTC", "TGA"]):
        codon = Codon(Seq(triplet),
                      FeatureLocation(cds.location.end - (idx * 3) - 3,
                                      cds.location.end - (idx * 3), -1),
                      idx)
        expected.append(codon)

    window = BestEditWindow(rec, cds, grna)
    ret = list(window.extract_codons())
    assert ret == expected
    atg = ret[0]
    assert str(rec.seq[atg.location.start:atg.location.end]) == "CAT"


def test_edit_window(rec, grnas):
    expected = FeatureLocation(34, 41, 1)
    ret = BestEditWindow.edit_window(grnas[0])
    assert ret == expected

    assert grnas[3].strand == -1
    expected = FeatureLocation(79, 86, -1)
    ret = BestEditWindow.edit_window(grnas[3])
    assert ret == expected


def test_overlaps(rec, cdses, grnas):
    cds = cdses[0]
    expected = [False for _ in range(len(grnas))]
    expected[1] = True
    expected[2] = True
    expected[3] = True
    expected[4] = True

    ret = list(map(lambda x: BestEditWindow.overlaps(cds, x), grnas))
    assert ret == expected


def test_can_edit(rec, grnas):
    expected = [False for _ in range(len(grnas))]
    expected[2] = True
    expected[3] = True
    expected[9] = True
    expected[10] = True


def test_get_affected_codons_fw_fw(rec, cdses, grnas):
    # CDS on +1, gRNA on +1
    edit_window = BestEditWindow(rec, cdses[0], grnas[2])
    expected = ["M1M", "A2V", "R3*"]

    ret = list(map(str, edit_window.get_affected_codons()))
    assert ret == expected


def test_get_affected_codons_fw_rv(rec, cdses, grnas):
    # CDS on +1, gRNA on -1
    edit_window = BestEditWindow(rec, cdses[0], grnas[3])
    expected = ["H10H", "A11T", "T12T"]

    ret = list(map(str, edit_window.get_affected_codons()))
    assert ret == expected


def test_get_affected_codons_rv_fw(rec, cdses, grnas):
    # CDS on -1, gRNA on +1
    edit_window = BestEditWindow(rec, cdses[1], grnas[9])
    expected = ["L14L", "F15F", "*16*"]

    ret = list(map(str, edit_window.get_affected_codons()))
    assert ret == expected


def test_get_affected_codons_rv_rv(rec, cdses, grnas):
    # CDS on -1, gRNA on -1
    edit_window = BestEditWindow(rec, cdses[1], grnas[10])
    expected = ["R3*", "R4C", "Y5Y"]

    ret = list(map(str, edit_window.get_affected_codons()))
    assert ret == expected


def test_mutations_codons_fw_fw(rec, cdses, grnas):
    # CDS on +1, gRNA on +1
    edit_window = BestEditWindow(rec, cdses[0], grnas[2])
    expected = ["A2V", "R3*"]

    ret = list(map(str, edit_window.get_mutations()))
    assert ret == expected


def test_get_mutations_fw_rv(rec, cdses, grnas):
    # CDS on +1, gRNA on -1
    edit_window = BestEditWindow(rec, cdses[0], grnas[3])
    expected = ["A11T"]

    ret = list(map(str, edit_window.get_mutations()))
    assert ret == expected


def test_get_mutations_rv_fw(rec, cdses, grnas):
    # CDS on -1, gRNA on +1
    edit_window = BestEditWindow(rec, cdses[1], grnas[9])
    expected = []

    ret = list(map(str, edit_window.get_mutations()))
    assert ret == expected


def test_get_mutations_rv_rv(rec, cdses, grnas):
    # CDS on -1, gRNA on -1
    edit_window = BestEditWindow(rec, cdses[1], grnas[10])
    expected = ["R3*", "R4C"]

    ret = list(map(str, edit_window.get_mutations()))
    assert ret == expected


def test_codon():
    atg = Codon(Seq("ATG"), FeatureLocation(0, 3, 1), 1)
    assert str(atg) == "ATG{[0:3](+)}(M1)"
    assert repr(atg) == str(atg)


def test_codon_mutate():
    aaa = Codon(Seq("AAA"), FeatureLocation(10, 13, 1), 1)
    ccc = Codon(Seq("CCC"), FeatureLocation(10, 13, 1), 1)

    nothing = ccc.mutate(FeatureLocation(6, 9, 1))
    assert nothing is None

    nothing = ccc.mutate(FeatureLocation(16, 19, 1))
    assert nothing is None

    full = ccc.mutate(FeatureLocation(10, 13, 1))
    assert str(full.seq) == "TTT"

    first_two = ccc.mutate(FeatureLocation(10, 12, 1))
    assert str(first_two.seq) == "TTC"

    last_two = ccc.mutate(FeatureLocation(11, 13, 1))
    assert str(last_two.seq) == "CTT"

    first = ccc.mutate(FeatureLocation(10, 11, 1))
    assert str(first.seq) == "TCC"

    last = ccc.mutate(FeatureLocation(12, 13, 1))
    assert str(last.seq) == "CCT"

    nothing = aaa.mutate(FeatureLocation(6, 9, 1), mode=AtoG)
    assert nothing is None

    nothing = aaa.mutate(FeatureLocation(16, 19, 1), mode=AtoG)
    assert nothing is None

    full = aaa.mutate(FeatureLocation(10, 13, 1), mode=AtoG)
    assert str(full.seq) == "GGG"

    first_two = aaa.mutate(FeatureLocation(10, 12, 1), mode=AtoG)
    assert str(first_two.seq) == "GGA"

    last_two = aaa.mutate(FeatureLocation(11, 13, 1), mode=AtoG)
    assert str(last_two.seq) == "AGG"

    first = aaa.mutate(FeatureLocation(10, 11, 1), mode=AtoG)
    assert str(first.seq) == "GAA"

    last = aaa.mutate(FeatureLocation(12, 13, 1), mode=AtoG)
    assert str(last.seq) == "AAG"


def test_codon_change():
    conservative = make_codon_change("CGC", "CGT", 18)
    assert conservative.conservative
    assert str(conservative) == "R7R"
    assert repr(conservative) == str(conservative)

    non_conservative = make_codon_change("CGA", "TGA", 0)
    assert not non_conservative.conservative
    assert str(non_conservative) == "R1*"

    assert conservative != non_conservative
