import os
import subprocess
import time

import dxpy

import keen

# These files are stored in the test project.
SAMPLE_FASTQ = "file-BQbXKk80fPFj4Jbfpxb6Ffv2"
HS37D5_BWA_INDEX = "file-B6ZY4942J35xX095VZyQBk0v"

def test_alignment_count(applet_id, project_id, folder, tmpdir):
    """Run BWA on a FASTQ file and verify that the number of
    alignments produced is correct.
    """

    # Recall that applet_id is set in the associated conftest.py, which either
    # gets it from the command line or builds the applet and retrieves its id.

    # And tmpdir is some pytest magic. It's type is py.path.local.LocalPath.
    # It's strpath property just returns a string.

    applet = dxpy.DXApplet(applet_id)
    input_dict = {"fastq": dxpy.dxlink(SAMPLE_FASTQ),
                  "genomeindex_targz": dxpy.dxlink(HS37D5_BWA_INDEX)}

    start_time = time.time()
    job = applet.run(input_dict, instance_type="mem1_ssd1_x16",
                     folder=folder, project=project_id)
    job.wait_on_done()
    end_time = time.time()
    elapsed_time = end_time - start_time

    output_bam_dxfile = dxpy.DXFile(job.describe()["output"]["bam"])
    local_filename = os.path.join(tmpdir.strpath, "test.bam")
    dxpy.download_dxfile(output_bam_dxfile.get_id(),
                         local_filename)
    count_alignments_cmd = "samtools view {bam} | wc -l".format(
        bam=local_filename)
    num_alignments = int(subprocess.check_output(count_alignments_cmd,
                                                 shell=True))
    keen.add_event("applet1_tests", {"num_alignments": num_alignments,
                                     "running_time": elapsed_time})
    assert job.describe()['state'] == 'done'
