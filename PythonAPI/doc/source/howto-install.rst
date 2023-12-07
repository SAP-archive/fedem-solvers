Installing fedempy
==================

To install the latest available ``fedempy`` package on Windows,
just give the following command in the console window::

    pip install fedempy

For this to work, you first need to copy the file
`pip.conf <https://github.com/SAP/fedem-solvers/blob/main/PythonAPI/pip.conf>`_
into `%APPDATA%\\pip\\pip.ini`, where `%APPDATA%` is the environment variable connected
to your user profile that points to where applications will have data and settings stored.
Notice the extension of the copied file needs to be `.ini` (as opposed to `.conf` on Linux).
If you you are running python in a virtual environment, the file can also be placed in the
root folder of the virtual environment file system instead.

Before using the ``fedempy`` modules, the following environment variables need to be set:

| **FEDEM_REDUCER** = path to the Fedem reducer shared object library (`fedem_reducer_core.dll`)
| **FEDEM_SOLVER** = path to the Fedem solver shared object library (`fedem_solver_core.dll`)
| **FEDEM_MDB** = path to the Fedem mechanism database shared object library (`FedemDB.dll`)

The current version of ``fedempy`` is #FEDEMPY_VERSION# and is compatible with **Fedem #FEDEM_VERSION#**.
