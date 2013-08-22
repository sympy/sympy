#!/usr/bin/env bash

# We are using the NiPy buildbot to run our time intensive tests: with all
# optional dependencies, without caching, and tests marked as slow.

DIR=`dirname $0`

# run all tests (all install dependencies should be installed)
$DIR/test

# run all the tests without the cache
$DIR/test --no-cache

# run only the slow tests, time out after 30 minutes
$DIR/test --slow --timeout=1800

# run slow tests without the cache?
#test --no-cache --slow --timeout=1800
