import datetime
from pathlib import Path
from unittest.mock import patch

import pytest

from rtl2 import regtest


@pytest.fixture
def runners():
    (args, test_list, suite) = regtest.setup(['MFIX-tests.ini'])
    suite._buildDir = Path("/top/build")
    with patch("rtl2.suite.datetime") as dt:
        dt.date.today.return_value = datetime.date(2020, 2, 2)
        assert suite.full_test_dir == Path("/top/build/2020-02-02-001")
    return regtest.get_runners(args, test_list, suite, None)
