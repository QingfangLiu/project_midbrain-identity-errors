import os
import numpy as np
import nibabel as nib
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score

def extract_beta_data(sub, local_dir, modelname, nconds, nallruns, lin_index):
    """
    Extract beta values for one subject.
    Returns: X (samples x voxels), y (labels)
    """
    subno = f'sub{sub}'
    print(f"Processing {subno}...")

    all_data = []
    all_labels = []

    for r in range(nallruns):
        run_data = np.zeros((nconds, len(lin_index[0])))  # shape: (nconds, nvox)

        for n in range(nconds):
            beta_idx = n + nconds * r + 1  # assumes 6 betas per run
            beta_path = os.path.join(local_dir, 'MultivariateGLM', subno, modelname, f'beta_{beta_idx:04d}.nii')

            vol = nib.load(beta_path).get_fdata()
            run_data[n, :] = vol[lin_index]

        all_data.append(run_data)
        all_labels.append(np.arange(1, nconds + 1))

    X = np.vstack(all_data)       # shape: (nconds * nallruns, nvox)
    y = np.hstack(all_labels)     # shape: (nconds * nallruns,)
    runs = np.repeat(np.arange(nallruns), nconds)

    return subno, X, y, runs


def decode_pairwise(X, y, runs, cond_pair, scale='none', clf_name='SVM', C=1.0):
    """
    Decode between two conditions using leave-one-run-out CV.
    Supports multiple classifiers via `clf_name`.
    """
    # Filter data for the condition pair
    mask = np.isin(y, cond_pair)
    X_pair = X[mask]
    y_pair = y[mask]
    runs_pair = runs[mask]

    if scale == 'demean':
        X_pair = X_pair - X_pair.mean(axis=1, keepdims=True)

    # Select classifier
    if clf_name == 'SVM':
        clf_base = SVC(kernel='linear', C=C)
    elif clf_name == 'LogReg':
        clf_base = LogisticRegression(C=C, solver='liblinear')
    elif clf_name == 'LDA':
        clf_base = LinearDiscriminantAnalysis()
    elif clf_name == 'Ridge':
        clf_base = RidgeClassifier(alpha=1.0)
    elif clf_name == 'RF':
        clf_base = RandomForestClassifier(n_estimators=100, max_depth=3, random_state=42)
    else:
        raise ValueError(f"Unknown classifier name: {clf_name}")

    # Leave-One-Run-Out CV
    logo = LeaveOneGroupOut()
    accs = []

    for train_idx, test_idx in logo.split(X_pair, y_pair, runs_pair):
        X_train, X_test = X_pair[train_idx], X_pair[test_idx]
        y_train, y_test = y_pair[train_idx], y_pair[test_idx]

        clf = clf_base.__class__(**clf_base.get_params())  # fresh clone per fold
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)

        accs.append(accuracy_score(y_test, y_pred))

    return np.mean(accs)



