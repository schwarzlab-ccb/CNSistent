.. CNSistent documentation master file, created by
   sphinx-quickstart on Tue Oct 29 18:11:36 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _introduction:

.. image:: files/Logo.png
   :alt: CNSistent Logo
   :width: 800px
   :align: center

|

Introduction
============

CNSistent is a Python tool for processing and analyzing copy number data. It is designed to work with data from a variety of sources. 
The tool is designed to be easy to use, and to provide a comprehensive set of analyses and visualizations.
CNSistent can be used as a Python package, or downloaded together with the respective data (PCAWG, TRACERx, TCGA, genomic locations):

Installation links
------------------
|PyPI version| |Documentation Status|

#. `Full GitHub repository with ~1GB of data. <https://github.com/ICCB-Cologne/CNSistent>`_
#. `PIP package only. <https://pypi.org/project/cnsistent/>`_


.. |PyPI version| image:: https://badge.fury.io/py/CNSistent.svg
   :target: https://badge.fury.io/py/CNSistent

.. |Documentation Status| image:: https://readthedocs.org/projects/cnsistent/badge/?version=latest
   :target: https://cnsistent.readthedocs.io/en/latest/?badge=latest

Example of the API
------------------

Files used below were adapted from `TRACERx Zenodo archive <https://zenodo.org/records/7649257>`_.

1. Load CNS Data and Display Heatmap
````````````````````````````````````

Load CNS data from a CSV file and visualize the first 5 rows using a heatmap.

.. code-block:: python

    import cns
    import cns.data_utils as cdu
    samples_df, raw_df = cdu.main_load("raw", "TRACERx")
    cns.fig_heatmap(cns.cns_head(raw_df, 5), max_cn=6)


.. image:: files/intro_1.png
    :alt: Raw Data Heatmap
    :width: 800px

2. Impute Missing Segments
``````````````````````````

Fill in missing segments in the data, impute using the extension method, and display a heatmap for the first 5 rows.

.. code-block:: python

    imp_df = cns.main_impute(raw_df, print_info=True)
    cns.fig_heatmap(cns.cns_head(imp_df, 5), max_cn=6)

.. image:: files/intro_2.png
    :alt: Imputed Data Heatmap
    :width: 800px

3. Create 3 mb Segments and convert to a feature array
``````````````````````````````````````````````````````

Aggregate the imputed CNS data into 3 MB segments and convert it into a feature array.

.. code-block:: python

    segs = cns.main_segment(split_size=3_000_000)
    seg_df = cns.main_aggregate(imp_df, segs, print_info=True)
    features, rows, columns = cns.bins_to_features(seg_df)
    print("Samples: {0}, Alleles: {1}, Bins: {2}.".format(*features.shape))

Printed output:

``Alleles: 2, samples: 403, bins: 960.``

4. Group Segments by Cancer Type
````````````````````````````````

Group the CNS data by cancer type, calculate the total CN, and visualize mean linear profiles.

.. code-block:: python

    type_groups = {c: cns.select_cns_by_type(seg_df, samples_df, c, "type") for c in ["LUAD", "LUSC"]}
    groups_df = cns.stack_groups([cns.group_samples(v, group_name=k) for k, v in type_groups.items()])
    cns.fig_lines(cns.add_total_cn(groups_df), cn_columns="total_cn")

.. image:: files/intro_3.png
    :alt: Grouped Data Heatmap
    :width: 800px

The example code is also in `example_API.py <https://github.com/ICCB-Cologne/CNSistent/blob/main/example_API.py>`_.

Example in terminal
-------------------

CNSistent reads SCNA profiles as ``.tsv`` files. Have an example file ``data.tsv``:

.. code-block:: python

    sample_id    chrom   start   end     total_cn
    sample1      chr1    100     200     1       
    ...

.. note::
    Column naming is fully describe in the :ref:`input_format` section.

To preprocess the segments:

.. code-block:: bash

    cns impute data.tsv --out imputed.csv

To create statistics:

.. code-block:: bash

    cns coverage data.tsv --out samples.tsv
    cns ploidy imputed.tsv --samples samples.tsv --out samples.tsv
    cns signatures imputed.tsv --samples samples.tsv --out samples.tsv

To calculate the mean ploidy per chromosome arm:

.. code-block:: bash

    segment arms --out arms.bed
    cns aggregate imputed.tsv --segments arms.bed --out a_bins.tsv

To conduct breakpoint clustering with 1 mb distance:

.. code-block:: bash

    segment imputed.tsv --merge 1000000 --out clust.bed
    cns aggregate imputed.tsv  --segments clust.bed --out c_bins.tsv

To conduct segmentation using 5 mb bins:

.. code-block:: bash

    segment whole --step 5000000 --out clust.bed
    cns aggregate data.tsv  --segments clust.bed --out c_bins.tsv

Extension of the example is in `example_CLI.sh <https://github.com/ICCB-Cologne/CNSistent/blob/main/example_CLI.sh>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   self
   quickstart
   CLI
   API
   reference


`LICENSE <https://github.com/ICCB-Cologne/CNSistent/blob/main/LICENSE.txt>`_
----------------------------------------------------------------------------