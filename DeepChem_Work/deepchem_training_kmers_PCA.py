"""
deepchem_training_kmers_PCA.py

Trains a DeepChem regression model to predict molecular weight of target molecules
using aptamer sequence features. The pipeline performs the following:

1. Extract k-mers from aptamer sequences
2. One-hot encode k-mers
3. Apply PCA for dimensionality reduction
4. Combine numerical features (GC content, length, binding affinity)
5. Train a multitask regression model using DeepChem
6. Evaluate using Pearson R² score on train / validation / test sets

Recommended execution environment: Google Colab with GPU runtime.
"""

import numpy as np
import pandas as pd
from itertools import product
from sklearn.preprocessing import OneHotEncoder
from sklearn.decomposition import PCA
import deepchem as dc


# ---------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------

def get_kmers(sequence, k=3):
    """Generate overlapping k-mers from a nucleotide sequence."""
    return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]


def preprocess_dataframe(filepath):
    """
    Load dataset, clean column names, generate k-mers, apply one-hot encoding,
    PCA reduction, and prepare datasets for DeepChem.
    """
    data = pd.read_excel(filepath)
    data.columns = data.columns.str.strip()

    k = 3
    bases = ['A', 'C', 'G', 'T']
    kmers = [''.join(p) for p in product(bases, repeat=k)]

    data['Kmers'] = data['Sequence'].astype(str).apply(lambda seq: get_kmers(seq, k))

    encoder = OneHotEncoder(sparse_output=False, handle_unknown='ignore', categories=[kmers])
    encoded_vectors = [
        encoder.fit_transform(np.array(kseq).reshape(-1, 1)).flatten()
        for kseq in data['Kmers']
    ]

    max_length = max(len(v) for v in encoded_vectors)
    padded = np.array([np.pad(v, (0, max_length - len(v))) for v in encoded_vectors])

    pca = PCA(n_components=500)
    reduced = pca.fit_transform(padded)

    numeric_columns = ['GC Content', 'Length', 'Binding Affinity']
    other_features = data[numeric_columns].apply(pd.to_numeric, errors='coerce').fillna(0)
    features = np.hstack([other_features.to_numpy(), reduced])

    targets = data[['Molecular weight']].apply(pd.to_numeric, errors='coerce').fillna(0)
    labels = data['Name'].astype(str).to_numpy()

    dataset = dc.data.NumpyDataset(X=features.astype(np.float32),
                                   y=targets.to_numpy().astype(np.float32),
                                   ids=labels)

    return dataset, features.shape[1], targets.shape[1]


# ---------------------------------------------------------------------
# Main procedure
# ---------------------------------------------------------------------

def main():
    filepath = "/content/xlfornumpy_cleaned.xlsx"  # Update before running

    dataset, n_features, n_tasks = preprocess_dataframe(filepath)

    splitter = dc.splits.RandomSplitter()
    train_dataset, valid_dataset, test_dataset = splitter.train_valid_test_split(
        dataset, frac_train=0.8, frac_valid=0.1, frac_test=0.1
    )

    model = dc.models.MultitaskRegressor(
        n_tasks=n_tasks,
        n_features=n_features,
        layer_sizes=[1000, 500],
        dropouts=[0.2, 0.2],
        learning_rate=0.001
    )

    print("Training model...")
    model.fit(train_dataset, nb_epoch=50)
    print("Training completed.\n")

    metric = dc.metrics.Metric(dc.metrics.pearson_r2_score)

    train_score = model.evaluate(train_dataset, metrics=[metric])
    valid_score = model.evaluate(valid_dataset, metrics=[metric])
    test_score = model.evaluate(test_dataset, metrics=[metric])

    print("Evaluation Results:")
    print("Train:", train_score)
    print("Validation:", valid_score)
    print("Test:", test_score)


if __name__ == "__main__":
    main()
