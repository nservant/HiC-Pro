import os
import tempfile
from iced.datasets import clear_data_home, get_data_home
from iced.datasets import load_sample_yeast
from nose.tools import assert_true, assert_false, assert_equal

DATA_HOME = tempfile.mkdtemp(prefix="hiclib_data_home_test_")
LOAD_FILES_ROOT = tempfile.mkdtemp(prefix="hiclib_load_files_test_")
TEST_CATEGORY_DIR1 = ""
TEST_CATEGORY_DIR2 = ""


def test_data_home():
    # get_data_home will point to a pre-existing folder
    data_home = get_data_home(data_home=DATA_HOME)
    assert_equal(data_home, DATA_HOME)
    assert_true(os.path.exists(data_home))

    # clear_data_home will delete both the content and the folder it-self
    clear_data_home(data_home=data_home)
    assert_false(os.path.exists(data_home))

    # if the folder is missing it will be created again
    data_home = get_data_home(data_home=DATA_HOME)
    assert_true(os.path.exists(data_home))

def test_data_sub_yeast():
    counts, lengths = load_sample_yeast()
