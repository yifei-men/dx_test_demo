import json
import os
import subprocess

import dxpy

import pytest

class TestArgumentError(Exception):
    pass

def dx_mkdir(project_id_, folder_):
    """Make a folder in the given project."""

    project_qualified_folder = ':'.join([project_id_, folder_])
    cmd = ['dx', 'mkdir', '-p', project_qualified_folder]
    print ' '.join(cmd)
    subprocess.check_call(cmd)

def dx_build(local_applet_dir, project_id_, folder_=None):
    """Build the applet located in the local_applet_dir to the
    given project_id. If folder is supplied, build that applet to
    that folder in the project. Otherwise, build to the root
    folder of the project.
    """

    cmd = ['dx', 'build', local_applet_dir]
    if folder_:
        destination = ':'.join([project_id_, folder_])
        with open(os.path.join(local_applet_dir, 'dxapp.json')) as dxapp:
            applet_name = json.load(dxapp)['name']
        destination = os.path.join(destination, applet_name)
    else:
        destination = project_id_

    cmd.extend(['--destination', destination])
    print ' '.join(cmd)
    id_ = json.loads(subprocess.check_output(cmd))['id']

    return id_

def get_applet_dir():
    """An applet should have a directory structure something like
    this:

    applet_root/dxapp.json
    applet_root/src/code.py
    applet_root/test/conftest.py
    applet_root/test/test_applet.py

    Well, we're here in conftest.py (Hey!). So the applet dir is just
    the parent dir of the dir that contains this file.
    """

    this_file_path = os.path.abspath(__file__)

    return os.path.dirname(os.path.dirname(this_file_path))

def get_project_id(project_id_, applet_id_):
    """Return the project_id to use for the tests.

    If project_id is specified, just use that. If not, use the
    the project that contains the applet.
    """
    if project_id_:
        return project_id_

    applet = dxpy.DXApplet(applet_id_)
    return applet.describe()['project']

def pytest_addoption(parser):
    """Declare the command line options for running the
    tests for this applet.
    """
    parser.addoption("--applet_id", action="store", default=None,
                     help="dxid of applet to be tested")
    parser.addoption("--project_id", action="store", default=None)
    parser.addoption("--folder", action="store", default=None)

@pytest.fixture(scope="session")
def applet_id(request):
    """applet_id is the applet that will be run. If it's not
    specified, then the applet will be built in the project
    specified by project_id.
    """
    id_ = request.config.getoption("--applet_id")

    # If we didn't get an applet_id from the command line, then we're going to
    # have to build the applet and get the id.
    if id_ is None:
        project_id_ = request.config.getoption("--project_id")

        if project_id_ is None:
            raise TestArgumentError("Either --applet_id or --project_id must be provided.")

        folder_ = request.config.getoption("--folder")
        if folder_:
            dx_mkdir(project_id_, folder_)
        id_ = dx_build(get_applet_dir(), project_id_, folder_)

    return id_

@pytest.fixture
def project_id(request):
    """This is the project in which to run the tests and perhaps build
    the applet.
    """
    id_ = request.config.getoption("--project_id")

    if id_ is None:
        applet_id_ = request.config.getoption("--applet_id")

        if applet_id_ is None:
            raise TestArgumentError("Either --applet_id or --project_id must be provided.")
        id_ = get_project_id(id_, applet_id_)

    return id_

@pytest.fixture(scope="session")
def folder(request):
    """This is the folder in which the tests should be run and applet
    built.
    """
    folder_ = request.config.getoption("--folder")

    if folder_:
        project_id_ = get_project_id(request.config.getoption("--project_id"),
                                     request.config.getoption("--applet_id"))
        dx_mkdir(project_id_, folder_)
    else:
        folder_ = "/"
    return folder_
