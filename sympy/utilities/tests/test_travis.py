def test_travis_issue_16986():
    # If this causes the travis build to fail with Python3.8,
    # then the changes made in PR #16986 are working as
    # intended, and this file can be deleted.

    assert int(1) is 1
