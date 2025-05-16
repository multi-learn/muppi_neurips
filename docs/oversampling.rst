Over-sampling
============================================


Method
------

Here, instead of using a sub-set of the majority class, we generate new examples from the minority class using `SMOTE <https://imbalanced-learn.readthedocs.io/en/stable/generated/imblearn.over_sampling.SMOTE.html>`_.

For the following results, the train set of each classifier is consists of real and SMOTE-generated examples, but the test set is restricted to the real examples. This allows us to test if, by learning on generated examples, we are able to improve the classification on the real ones.

Accuracy
--------

As explained earlier the test set of this task consists only on real examples, so it is as un-balanced as the original dataset. The accuracy here is barely relevant, as the majority of the dataset is comprised of mulit clustered proteins. 

.. raw:: html
    :file: ./images/over_sampling_results/multi_vs_emf_over_acc.html






F1-score
--------

The F1-score is far more relevant in this situation, and we can see on the Figure below that the best algorithm is the early fusion.

.. raw:: html
    :file: ./images/over_sampling_results/multi_vs_emf_over_f1.html



Comparison with un-balanced data
-------------------------------

Performance
<<<<<<<<<<<

The performance improvement with SMOTE is not as impressive as with sub-sampling, the multiview algorithms are able to extract more information from the generated examples than from the base dataset.

Duration
<<<<<<<<

Running **SuMMIT** on the over-sampled dataset takes more time than on the raw data, as SMOTE generates more examples. However, depending on the algorithmic complexity of the used algorithm, this difference can vary from slight (with a decision tree, for example) to considerable for algorithms using quadratic program solvers. Moreover, as the over-sampled dataset is unique, is is still faster to run **SuMMIT** on it than on the 10 sub-sampled ones.


