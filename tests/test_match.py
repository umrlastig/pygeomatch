import pytest
from pygeomatch.match import merge_groups
from pygeomatch.util import MatchingLink

def test_merge_groups():
    l1 = MatchingLink(0,0,0,{})
    l2 = MatchingLink(0,1,0,{})
    l3 = MatchingLink(1,1,1,{})
    l4 = MatchingLink(1,2,1,{})
    # 2 groups that should be merged
    links = merge_groups([l1, l2, l3, l4])
    print(links)
    for l in links:
        assert l.group == 0

def test_tuple():
    ref, comp, group, measures = MatchingLink(1,2,3,{"measure": 4.}).as_tuple()
    assert ref == 1
    assert comp == 2
    assert group == 3
    assert measures == {"measure": 4.}