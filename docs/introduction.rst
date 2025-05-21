Introduction Benchmark of MuPPI Dataset
========================================
MuPPI dataset (in NoNA release) aggregates 14 description on 8749 proteins.

Here, we present several benchmarks on this dataset, performed with **SUpervised MultiModal Integration Tool**
(**SuMMIT** platform). Its documentation is available on
`<https://gitlab.lis-lab.fr/baptiste.bauvin/multiview-machine-learning-omis>`_
while the sources are available `here <http://dev.pages.lis-lab.fr/multiview-machine-learning-omis/>`_

This platform aims at running multiple state-of-the-art classifiers to get a baseline on common algorithms for any classification task.

Algorithms
<<<<<<<<<<<<<<<<<<<<<<<<<<

To set a baseline, the selected algorithms are:

+ Monoview (on each view)
    - from `**scikit-learn** <https://scikit-learn.org/stable/index.html>`_
        * `DecisionTreeClassifier <https://scikit-learn.org/stable/modules/generated/sklearn.tree.DecisionTreeClassifier.html>`_
        * `RandomForestClassifier  <https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html#sklearn.ensemble.RandomForestClassifier>`_
    - from `imbalanced-learn <https://imbalanced-learn.readthedocs.io/en/stable/>`_
        * `Bagging <https://imbalanced-learn.readthedocs.io/en/stable/generated/imblearn.ensemble.BalancedBaggingClassifier.html#imblearn.ensemble.BalancedBaggingClassifier>`_ (only for un-balanced tasks)

+ Multiview algorithms
    - late fusion, that learns a monoview algorithm for each view and combine their decisions with a majority vote,
    - early fusion, that learns a monoview algorithm on the concatenation of all the views,
    - from  `scikit-multimodallearn <http://dev.pages.lis-lab.fr/scikit-multimodallearn/>`_
        * Mumbo,
        * MuCumbo.

.. note::
    Bagging, Mumbo, and MuCombo are not available in the master branch of **SuMMIT**, but are upon request.

Learning Tasks
<<<<<<<<<<<<<<<

The proteins are labeled as `mono_clustered` (6139 proteins), `multi_clustered` (2413 proteins)
and `EMF` (197 proteins)

As showed in the following results, we analyzed three tasks :

+ :ref:`EMF versus mono_clustered <emf_vs_mono>`, highly un-balanced, but fairly easy with only one view,

+ :ref:`multi_clustered versus mono-clustered <multi_vs_mono>` less unbalanced, but slightly harder than the previous one,

+ :ref:`EMF versus multi-clustered <emf_vs_multi>` is the most difficult task, and the most relevant biologically, so we will analyze two reb-balancing methods for this task: `sub-sampling <sub_sampling>`_ and `over-sampling <oversampling>`_

.. note::
    For each task, a config file for **SuMMIT** is available in the ``./config_files`` directory.
    The hyper-parameters have been fixed to the best of our knowledge, but we are aware that they still can be optimized.
    We used 75% of the dataset to train each classifier and the remaining 25% to test the classifiers.

Results summary
<<<<<<<<<<<<<<<

+--------------------------+----------------+----------------+---------------+--------------------+
| Task                     | Max. Accuracy  | Max. f1-score  |  Algorithm    |  View              |
+==========================+================+================+===============+====================+
| EMF vs Mono              | 0.99           | **0.82**       | MuCombo       | **Multiview**      |
+--------------------------+----------------+----------------+---------------+--------------------+
| Multi vs Mono            | 0.87           | **0.76**       | Random Forest | PPINetwork_topology|
+--------------------------+----------------+----------------+---------------+--------------------+
| EMF vs Multi base        | 0.89           | **0.37**       | Late Fusion   | **Multiview**      |
+--------------------------+----------------+----------------+---------------+--------------------+
| EMF vs Multi sub-sample  | 0.78           | **0.80**       | Decision Tree | PPINetwork_topology|
+--------------------------+----------------+----------------+---------------+--------------------+