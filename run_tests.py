#!/usr/bin/env python
"""Run pytest tests on applets found in the pwd."""
import argparse
import json
import os
import re
import subprocess
import sys
import time
import dxpy

import pytest

def get_applet_dirs(rootdir):
    """Find all applet directories within rootdir.

    An applet directory is just a directory with a file called
    dxapp.json in it.
    """

    applet_dirs = []
    for root, dirs_ignored, files in os.walk(rootdir):
        if "dxapp.json" in files:
            applet_dirs.append(re.sub("^\.\/", '', root))

    return applet_dirs

def get_git_changes(commit_range):
    """Get a list of files changed in the given git commit range."""

    git_cmd = ['git', 'diff', '--name-only', commit_range]

    changes = filter(None, subprocess.check_output(git_cmd).split('\n'))

    return changes

def remove_unchanged_dirs(applet_dirs, git_changes):
    """Get the list of applet dirs that contain files changed
    in the relevant git range.
    """

    changed_dirs = []

    for dir_ in applet_dirs:
        if any(g.startswith(dir_) for g in git_changes):
            changed_dirs.append(dir_)

    return changed_dirs

def get_revision_tag():
    """Return a string that uniquely identifies the applet being tested.

    The applet is identified by the current date and time as well as the git
    commit hash of HEAD.
    """

    # Use a consistent time zone when naming folders
    os.environ["TZ"] = "America/Los_Angeles"
    time.tzset()

    date_string = time.strftime("%Y-%m-%d-%H.%M.%S")
    try:
        git_revision = subprocess.check_output(["git", "rev-parse", "HEAD"]).strip()
    except subprocess.CalledProcessError:
        git_revision = "none"

    return date_string + "-" + git_revision

def run_test(applet_dir, project_id):
    """Build the applet in applet_dir and run its test in the
    project specified by project_id.
    """

    test_folder = os.path.join("/applets", applet_dir,
                               get_revision_tag())
    test_args = [applet_dir, "--project_id", project_id, '--folder',
                 test_folder]

    print "test args:", ' '.join(test_args)
    returncode = pytest.main(test_args)

    return returncode

def run_static_analyzer(applet_dir):
    """Run a static analyzer, pylint or shellcheck, on the source file of
    the given applet.
    """
    dxapp_json = json.load(open(os.path.join(applet_dir, 'dxapp.json')))
    interpreter = dxapp_json['runSpec']['interpreter']
    source_file = os.path.join(applet_dir, dxapp_json['runSpec']['file'])

    if interpreter.startswith("python"): # Handle python2.7, etc.
        run_pylint(source_file)
    elif interpreter == 'bash':
        input_names = [str(e['name']) for e in dxapp_json['inputSpec']]
        run_shellcheck(source_file, input_names)

def run_pylint(source_file):
    """Run pylint on the source file."""
    pylint_cmd = ['pylint', source_file]
    # Disable heresies
    pylint_cmd.append("--disable=W0141,W0142")

    proc = subprocess.Popen(pylint_cmd, stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    print stdout

def run_shellcheck(source_file, supplied_inputs=None):
    """Run shellcheck on the source file.

    If supplied_inputs is specified, "referenced but not assigned errors" for those
    variables will be suppressed.
    """
    shellcheck_cmd = ["shellcheck", '--format', 'gcc', source_file]

    proc = subprocess.Popen(shellcheck_cmd, stdout=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    shellcheck_results = filter(None, stdout.split("\n"))

    # Shellcheck will complain about variables that are referenced but not
    # assigned. But some of those are inputs from dxapp.json. So we will filter
    # those out to reduce spurious errors.
    cleaned_results = []
    for result in shellcheck_results:
        if "[SC2154]" in result: # Is a referenced but not assigned error
            missing_name = re.search(
                "(?<=warning: )\w+(?= is referenced but not assigned)", result).group()
            if supplied_inputs and missing_name not in supplied_inputs:
                cleaned_results.append(result)
        else:
            cleaned_results.append(result)

    for result in cleaned_results:
        print result.strip()

def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--project_id",
                        required=True,
                        help=("The project in which to build the applets and "
                              "run their tests."))
    parser.add_argument("--commit_range",
                        required=False,
                        default=None,
                        help=("Limit tested applets to the those that were changed "
                              "in the given git commit range."))
    args = parser.parse_args()

    applet_dirs = get_applet_dirs(os.path.curdir)
    print "Found applet dirs:", applet_dirs

    if args.commit_range is not None:
        files_changed_in_commit_range = get_git_changes(args.commit_range)
        applet_dirs = remove_unchanged_dirs(
            applet_dirs, files_changed_in_commit_range)
        print "After filtering by commit range:", applet_dirs

    returncodes = []
    for applet_dir in applet_dirs:
        run_static_analyzer(applet_dir)
        returncode = run_test(applet_dir, args.project_id)
        returncodes.append(returncode)

    sys.exit(any(r != 0 for r in returncodes))

if __name__ == '__main__':
    main()
