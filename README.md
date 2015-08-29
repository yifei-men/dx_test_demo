[![Build
Status](https://travis-ci.org/mckinsel/dx_test_demo.svg?branch=master)](https://travis-ci.org/mckinsel/dx_test_demo)

# Testing Demo

## Use cases

### Testing a single applet

Each applet has a `test` directory that contains functional tests for the
applet. The tests can be run using [pytest](www.pytest.org).

For example,

    py.test applet1 --project_id project-BgKPX680Y028gX2xxVYzvy7b

will build applet1 in the given project. Then, it will look for tests in the
applet1 directory using the usual Python test discovery rules. For applet1, it
will find `applet1/test/test_applet1.py`. This contains a function,
`test_alignment_count`, that runs the newly-built applet with some sample
inputs, downloads the resulting BAM file to a temporary directory, and checks
the number of alignments in the BAM file against the expected value.

Tests can also be run against an applet that's already been built:

    py.test applet1 --applet_id applet-BgQx2100Y02JBJKF5540k7BX

Rather than building a new applet, this command will run the tests using the
given applet ID. The tests will be run in the project that contains the applet.

The previous commands all used the project's root folder, but a different
folder can be specified:


    py.test applet1 --project_id project-BgKPX680Y028gX2xxVYzvy7b --folder /a/test/folder

This command will build a new applet in the `/a/test/folder` folder of the
project and run the tests there. If the folder does not exist, it will create
it.

### Testing multiple applets

The script `run_tests.py` looks for applets in the repo, builds them, and runs
their tests.

    run_tests.py --project_id project-BgKPX680Y028gX2xxVYzvy7b

will walk through the current directory, finding applet directories. For each
applet directory, it will build the applet and run its tests. The applet and
test outputs are placed in a folder named like this:

    /applets/{applet_dir}/{date}-{time}-{git_hash}

where `applet_dir` is the relative path to the local directory that contains
the applet, `date` and `time` is when the test is run, and `git_hash` is the
hash of the latest commit of the repo that contains the applet.

Since some repos contain many applets, you may want to only run tests for the
applets that have changed recently. This can be done with the `commit_range`
parameter:

    run_tests.py --project_id project-BgKPX680Y028gX2xxVYzvy7b --commit_range 962dbf8..f3cf772

This will find the changed files in the commit range and run tests for applet
directories that contain at least one of those files.

### Continuous integration

The `.travis.yml` file of this repo is set up to build and run tests on each
push. When changes are pushed, it will find the applets that are different, run
their tests, and report the results using whatever Travis notification scheme
has been set up.

```
git clone git@github.com:mckinsel/dx_test_demo.git
cd dx_test_demo
sed -i 's/sleep 1/sleep 2/' applet1/src/applet1.sh
git commit applet1/src/applet1.sh
git push
```

This would cause Travis to build and test applet1 in the project specified in
the `.travis.yml` file. The project is currently a public project called
test-demo-CI.

### Static analysis

Nearly all applet code is written in python or bash. It is easy to write code
that has non-obvious problems in either language, but this is very, very true
for bash. To help find these, `run_tests.py` runs a linting tool on the applet
source file for applets that it tests. For python it runs
[pylint](http://www.pylint.org/), for bash it runs
[shellcheck](http://www.shellcheck.net/about.html). It reports the results of
these tools along with the test results, so you get something like this:

```
applet1/src/applet1.sh:10:17: note: Use $(..) instead of legacy `..`. [SC2006]
applet1/src/applet1.sh:13:16: warning: Quote this to prevent word splitting.
[SC2046]
applet1/src/applet1.sh:13:16: note: Use $(..) instead of legacy `..`. [SC2006]
applet1/src/applet1.sh:13:102: warning: Quote this to prevent word splitting.
[SC2046]
applet1/src/applet1.sh:13:102: note: Use $(..) instead of legacy `..`. [SC2006]
```

## Writing tests

### conftest<span></span>.py

Each test directory has a file called `conftest.py`. That file contains code
for parsing pytest command line arguments and creating arguments for test
functions. In its present form, it provides tests three arguments:

* `applet_id` - the ID of the applet that is going to be tested. Tests should
  usually call some sort of `dxpy.DXApplet(applet_id).run()`.
* `project_id` - the ID of the project where tests should be run. This should
  usually be passed as the `project` argument of `dxpy.DXApplet.run`.
* `folder` - the platform folder where test results should be written.

### Pytest

In pytest, most tests can just use `assert` statements. Pytest also lets
functions use a `tmpdir` argument which is a temporary directory that pytest
manages.

Pytest syntax is pretty liberal. You just write some functions or class methods
with asserts.

```python
def test_something(applet_id, project_id, folder):
    applet = dxpy.DXApplet(applet_id)
    job = applet.run({"input_file": input_file_link},
                     project=project_id,
                     folder=folder)
    job.wait_on_done()

    assert job.describe()['output']['output_value'] > 10 # or whatever
```

## Known issues

* Only a subset of applets are tested for each push, and Travis builds are
  independent. So, if a breaking change to applet1 is pushed, Travis will mark
  the build as broken. But if a non-breaking change to applet2 is then pushed,
  Travis will mark the build as fixed because only tests for applet2 were run,
  and they passed. But, applet1 is still broken.

* Even if the tests are very simple with very small inputs, they will still
  take a long time to run because of platform overhead. This makes this scheme
  less useful for a developer who wants to run tests quickly and frequently. 

* Everything is run in serial right now. This is easy enough to fix once an
  appropriate parallelization scheme is decided on.

* Some of the tests will look idle to Travis because all the work is happening
  on the platform. Travis will kill these tests after a relatively short period
  of time.

* Every test name has to be unique. If they are not, you get weird import
  errors from pytest.

* Tests can have lots of dependencies, and this scheme does not have a way to
  specify them.
