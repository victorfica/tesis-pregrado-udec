# -*- coding: utf-8 -*-


from sklearn.feature_selection.base import SelectorMixin
from sklearn.base import BaseEstimator
import numpy as np



class ImportanceSelect(BaseEstimator, SelectorMixin):
      
    def __init__(self, model, n=1):
        self.model = model
        self.n = n
    def fit(self, *args, **kwargs):
        self.model.fit(*args, **kwargs)
        return self
    def _get_support_mask(self):
        #check_is_fitted(self, 'feature_importances_')

        if self.n == 'all':
            return np.ones(self.model.feature_importances_.shape, dtype=bool)
        elif self.n == 0:
            return np.zeros(self.model.feature_importances_.shape, dtype=bool)
        else:
            scores = self.model.feature_importances_
            mask = np.zeros(scores.shape, dtype=bool)

            # Request a stable sort. Mergesort takes more memory (~40MB per
            # megafeature on x86-64).
            mask[np.argsort(scores, kind="mergesort")[-self.n:]] = 1
        return mask
            
    
    def transform(self, X):
        return X[:,self.model.feature_importances_.argsort()[::-1][:self.n]]