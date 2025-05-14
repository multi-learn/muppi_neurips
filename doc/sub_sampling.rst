Sub-sampling
============

To re-balance the dataset, and try to obtain more information with usual classifiers that are not built to deal with unbalanced problems, we used sub-sampling.

This method consists on selecting a subset of the majorrty class (here `multi_clustered`) to generate a new balanced task. To do so, we rebalanced with one multi_clustered for each EMF and we generated 10 randomly sub-sampled datasets on which we learned the classifiers.

The results shown here are the mean and standard deviation over these 10 sub-samplings, to avoid any lucky pick.



Accuracy
--------

The classifiers are sorted by test accuracy below, the best classifier is `MuCombo <http://dev.pages.lis-lab.fr/scikit-multimodallearn/reference/api.html#module-multimodal.boosting.cumbo>`_, but more importantly, the multiview classifiers are slightly better than all the monoview ones except for MVML. Indeed, as it is a kernel-based algorithm, hyper-parameter optimization plays a major role in performance, however, they were not optimized for each sub sampled dataset, to reduce computation time.

.. raw:: html
    :file: ./images/sub_sampling_results/multi_vs_emf_sub_acc.html


.. note::
    As the 10 sub-sampled datasets are re-balanced to 1:1, the accuracy is now relevant.

F1-score
--------

F1-score is still relevant as it is based on the the positive class examples (EMFs) that are well classified, and does not take the well-classified multi-clustered into account.

.. raw:: html
    :file: ./images/sub_sampling_results/multi_vs_emf_sub_f1.html

The order of the classifiers is still nearly the same, MVML has a far higher rank, and the top five algorithms are still the same, even if they are ordered differently.



Comparison with un-balanced data
--------------------------------

Performance
<<<<<<<<<<<

The gap in performance between the un-balanced and the sub-sampled datasets is obvious. The best f1-score on the raw data is 0.46, as in the sub-sampling setting, the fifth lowest is already higher, and the best f1-score is 0.77.

Sub-sampling seems to be a good mean to extract meaningful information from this dataset.

Duration
<<<<<<<<

However, it is far longer to run each algorithm on the sub-sampled datasets. Indeed, running **SuMMIT** on one took 37 mins, so if not parallelization is available, running on the 10 sub-sampled datasets takes 370 mins (6 hours and 10 mins), while running **SuMMIT** on the raw dataset took 54 mins.

.. note::
    The durations are for reference only as they depend on the machine and the algorithms used for the benchmark.

Features Importance ??
----------------------

This experiment was designed to be able to extract the features on which the classifiers based their models to be able to interpret the decisions. However, ...



Samples analysis ??
-------------------

Don't remember what that was about ?!