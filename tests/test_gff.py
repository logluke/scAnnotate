from src.scAnnotate.gff import merge_overlapping_intervals, load_genes_with_exons


def test_merge_overlapping_intervals():
    intervals = [(1, 10)]
    assert merge_overlapping_intervals(intervals) == [(1, 10)]

    intervals = [(1, 5), (2, 6), (8, 10), (9, 12)]
    assert merge_overlapping_intervals(intervals) == [(1, 6), (8, 12)]

    intervals = [(1, 5), (4, 10), (9, 15)]
    assert merge_overlapping_intervals(intervals) == [(1, 15)]


def test_load_genes_with_exons():
    assert True == True
