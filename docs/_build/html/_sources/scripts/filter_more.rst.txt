******************
:mod:`filter_more`
******************

Filter a csv file with a given trajectory and get the ``n`` most important derivatives
for a given out_param.
Usage:

- ``-i``, ``--input``: Path to csv file.
- ``-n``: The number of derivatives on output. If n=-1 no filtering is being done.
- ``-o``, ``--output``: Path and name where to store the filtered file.
- ``-p``, ``--param``: The out parameter to filter for.
- ``t``, ``--trajectory``: The trajectory to filter for.

.. automodule:: scripts.filter_more
   :members: