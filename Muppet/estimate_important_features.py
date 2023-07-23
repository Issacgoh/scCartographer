
import warnings
import itertools
import math
import pandas as pd

class estimate_important_features:
    def __init__(self, model, top_n):
        print('Estimating feature importance')
        classes =  list(model.classes_)
        try:
            model_features = list(itertools.chain(*list(model.features)))
        except:
            warnings.warn('no features recorded in data, naming features by position')
            print('if low-dim lr was submitted, run linear decoding function to obtain true feature set')
            model_features = list(range(0,model.coef_.shape[1]))
            model.features = model_features
        print('Calculating the Euler number to the power of coefficients')
        impt_ = pow(math.e,model.coef_)
        try:
            self.euler_pow_mat = pd.DataFrame(impt_,columns = list(itertools.chain(*list(model.features))),index = list(model.classes_))
        except:
            self.euler_pow_mat = pd.DataFrame(impt_,columns = list(model.features),index = list(model.classes_))
        self.top_n_features = pd.DataFrame(index = list(range(0,top_n)))
        print('Estimating feature importance for each class')
        mat = self.euler_pow_mat
        for class_pred_pos in list(range(0,len(mat.T.columns))):
            class_pred = list(mat.T.columns)[class_pred_pos]
            temp_mat =  pd.DataFrame(mat.T[class_pred])
            temp_mat['coef'] = model.coef_[class_pred_pos]
            temp_mat = temp_mat.sort_values(by = [class_pred], ascending=False)
            temp_mat = temp_mat.reset_index()
            temp_mat.columns = ['feature','e^coef','coef']
            temp_mat = temp_mat[['feature','e^coef','coef']]
            temp_mat.columns =str(class_pred)+ "_" + temp_mat.columns
            self.top_n_features = pd.concat([self.top_n_features,temp_mat.head(top_n)], join="inner",ignore_index = False, axis=1)
            self.to_n_features_long = model_feature_sf(long_format_features(self.top_n_features),'e^coef')
