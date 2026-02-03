# MFiX Regression Testing (Level 2)

regtest.py is the main regression testing script for MFiX-Exa Level 2 Regression Tests.
Typing './regtest.py -h' will give a verbose
description of usage and setup.  The regtest set of tools is
fully documented in the AMReX Users Guide, but here's a quick
overview of usage.

## Quickstart

```
python regtest.py MFIX-tests.ini
```

This will build MFiX under `rt-MFIX-Exa/build` and run the [RTL2] tests in the
config file. When the tests finish, an HTML report of the results will be
generated at `rt-MFIX-Exa/web/index.html`

## Details

1. Each regression test suite is defined in a separate directory
2. MFIX-tests.ini is the default, which builds the rtl2 tests under mfix/benchmarks/rtl2.

 a. In the [source] config blocks of the .ini file, the location of the MFiX
    tests are set as "dir".

 b. Parameters that define each regression test appear in labeled blocks below
    the repository blocks. The number of options allowed for specifying these
    tests will certainly grow over time; the -h option to the tester should list
    everything currently available.

 c. Each labeled test will be built, run and compared to "benchmark" results.
    The results of these tasks will be formatted into a web page in a folder
    defined as "webTopDir" in the .ini file. Make sure that variable points to
    something you like and can access.

 d. Tests are run by typing

    ./regtest.py MFIX-tests.ini

