#!/usr/bin/env python3

"""
INTRO:
    1)calculate_probability is used to calculate decision scores
    2)print_scores is used to print scores, such as precision score
    3)plot_roc_curve is used to plot roc curve
"""


def calculate_probability(input_classifier, x, target):
    
    '''
    A faster way to get decision_scores for multiple kinds of classifiers including lnl

    -input_classifier classifier used to calculate decision score
    -x standard features
    -target true label
    '''

    from cleanlab.latent_estimation import estimate_confident_joint_and_cv_pred_proba
    # proba = input_classifier.predict_proba(x)[:, 1]
    # decision_scores = np.log(proba / (1 - proba))
    try:
        y_scores = input_classifier.decision_function(x)
    except:
        confident_joint, y_scores = estimate_confident_joint_and_cv_pred_proba(
            X=x, 
            s=target,
            clf=input_classifier,
            cv_n_folds= 2)

        y_scores = y_scores[:,1]
    return y_scores

def print_scores(y_true, y_predicted):

    """
    Print precision/recall/f1/confusion matrix/roc_auc
    """

    import pandas as pd
    from sklearn.metrics import precision_score,recall_score,f1_score, accuracy_score, roc_curve, roc_auc_score,confusion_matrix

    print('Info: Accuracy score: ', accuracy_score(y_true, y_predicted))
    print('Info: Precision score: ',precision_score(y_true, y_predicted))
    print('Info: Recall score: ', recall_score(y_true, y_predicted))
    print('Info: f1 score: ', f1_score(y_true, y_predicted))
    print('Confusion Matrix:')
    conf_mtrix = pd.DataFrame(confusion_matrix(y_true, y_predicted), columns=['Predicted Negtive', 'Predicted Positive'], index=['Actual Negtive', 'Actual Postive'])
    print(conf_mtrix.to_string())
    roc_auc = roc_auc_score(y_true, y_predicted)
    print('Info: ROC AUC score:', roc_auc)

    
def plot_roc_curve(y_true, decision_scores, label=None):
    """
    Plot one roc curve when giving fpr and tpr
    This is only for models generation
    """
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve

    fpr, tpr, thresholds = roc_curve(y_true, decision_scores)
    plt.plot(fpr, tpr, linewidth=2, label=label)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel = 'False Positive Rate'
    plt.ylabel = 'True Positive Rate'
    plt.legend()