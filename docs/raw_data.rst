Results of Benchmak with Raw Data
=================================

**SuMMIT** is launched directly on the data over the 14 views described in `./dataset`.

For each experiment described here, we ran each algorithm on 5 train/test splits to avoid a lucky split.

Moreover, for each task, we present several figures :

* An accuracy bar plot that shows the train and test accuracy for each classifier
* The same graph for f1-score, to take unbalancing into account,
* An error graph that shows on which example each classifier failed.

.. note::
    All the following figures are interactive, so you can zoom in and hover to get the information you need.

.. _emf_vs_mono:
First task : 'EMF' versus `mono-clustered`
------------------------------------------

This task is the most unbalanced one (5306 `mono_clustered` vs 181 `EMF`), so the following accuracy graph has not much meaning as a classifier that classify every protein as `mono_clustered` has an accuracy score of 0.97.


.. raw:: html
    :file: ./images/raw_results/mono_vs_emf_acc.html


Therefor, the f1-score is a much more relevant metric on an unbalanced dataset,
to be able to see which classifier can classify the minority of EMF. The high f1-score of the monoview classifiers on the `PPINetwork_topology` view are explained by the fact that biologists classify `mono_clustered` protein based on the PPI network, so the information of this view is sufficient to classify `mono_clustered` proteins vs `EMF`


.. raw:: html
    :file: ./images/raw_results/mono_vs_emf_f1.html

The following figure is the representation of a matrix that shows on it rows the examples of the dataset and on its columns the classifiers of the benchmark. For each example and classifier, a black rectangle means that the classifier was always wrong when it had to classify the example, and a white one means it was always right, and between them are shades of grey.

Here, we can clearly see that the vast majority of classifiers have are frequently wrong on the `EMF` (the first rows), but that some of them classify well all the examples (the multiview ones, and the `imbalance_bagging` mainly.)

.. raw:: html
    :file: ./images/raw_results/mono_vs_emf_err.html

As this task is highly unbalanced and not really biologically relevant , we study a more evenly-balanced case : `mono_clustered` vs `multi_clustered`.

.. _multi_vs_mono:
Second task : `multi_clustered` vs `mono_clustered`
---------------------------------------------------

This task is more evenly balanced than the previous with 5306 `mono_clustered` and 2205 `multi_clustered`, so the accuracy is a bit more relevant, but analyzing both metrics is still mandatory.

.. raw:: html
    :file: ./images/raw_results/mono_vs_multi_acc.html


Similarly to the previous task, as the biological classification between mono and multi clustered proteins is done using the PPI network, the task is easily solved with only the `PPINetwork_topology` view.


.. raw:: html
    :file: ./images/raw_results/mono_vs_multi_f1.html

The following matrix shows similar results to the previous task even if the dataset is less un-balanced, some classifiers are unable to classify any multi_clustered (nearly all decision trees). Where the random forest on PPINetwork_topology is nearly perfect, similarly to the multiview algorithms.

.. raw:: html
    :file: ./images/raw_results/mono_vs_multi_err.html

The last task is the most difficult one as it is the most relevant biologically, however it is still highly unbalanced (less thant the first).

.. _emf_vs_multi:
Last task : EMF versus multi-clustered
--------------------------------------

This task is difficult as, while the EMF labelled proteins are surely EMFs, the multi_clustered class may contain some EMFs that were not biologically discovered.

The following figure shows the accuracy of every classifier, however, this metric is barely relevant as the task is highly unbalanced.

.. raw:: html
    :file: ./images/raw_results/multi_vs_emf_acc.html


The relevant metric in this case is the f1-score.
We can see below that it is here particularly low, meaning that barely any classifier is able to extract information about the asked task form the data.

In order to tackle this issue, re-balancing the dataset could remove a part of the problem's difficulty and improve scores.

.. raw:: html
    :file: ./images/raw_results/multi_vs_emf_f1.html

It is clear on the following figure that the un-balancing of the dataset causes algorithms that does not take it into account to mis-classify nearly all the EMFs and causes the ones that take it into account to mis-classify a lot of multi-clustered to be able to build a model including information about the EMFs.

.. raw:: html
    :file: ./images/raw_results/multi_vs_emf_err.html

To try to facilitate the task, we propose two ways of re-balancing the dataset : `sub-sampling <sub_sampling>`_ and `over-sampling <oversampling>`_.