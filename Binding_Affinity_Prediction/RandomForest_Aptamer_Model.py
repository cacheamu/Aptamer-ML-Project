"""
RandomForest_Aptamer_Model.py

Trains a Random Forest regression model to predict binding affinity of aptamers
using sequence-derived features including:

• k-mer composition (k = 2)
• GC content
• molecular weight
• sequence length

The script:
1. Loads aptamer data from Excel
2. Extracts numerical + sequence features
3. Trains a RandomForestRegressor
4. Evaluates using RMSE
5. Saves the trained model and DictVectorizer for future predictions
"""

import pandas as pd
import numpy as np
from collections import Counter
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from sklearn.feature_extraction import DictVectorizer
import joblib


# ---------------------------------------------------------------------
# Feature extraction functions
# ---------------------------------------------------------------------

def get_kmers(seq, k=2):
    seq = seq.replace("5'-", "").replace("-3'", "").replace(" ", "").upper()
    return [seq[i:i + k] for i in range(len(seq) - k + 1)]


def gc_content(seq):
    seq = seq.upper()
    if len(seq) == 0:
        return 0
    return (seq.count('G') + seq.count('C')) / len(seq)


def molecular_weight(seq):
    weights = {'A': 313.21, 'T': 304.2, 'G': 329.21, 'C': 289.18}
    return sum(weights.get(base, 0) for base in seq.upper())


# ---------------------------------------------------------------------
# Main ML pipeline
# ---------------------------------------------------------------------

def main():
    filepath = "/content/finalexcelsheetapta19mar(5).xlsx"  # Update before running
    sheet = "main aptamers"

    df = pd.read_excel(filepath, sheet_name=sheet)

    df = df.dropna(subset=["Sequence", "Binding Affinity"])
    df["Sequence"] = df["Sequence"].astype(str)
    df["Binding Affinity"] = pd.to_numeric(df["Binding Affinity"], errors="coerce")
    df = df.dropna(subset=["Binding Affinity"])

    df["Length"] = df["Sequence"].apply(len)
    df["GC Content"] = df["Sequence"].apply(gc_content)
    df["Molecular weight"] = df["Sequence"].apply(molecular_weight)

    kmer_features = []
    for seq in df["Sequence"]:
        kmers = get_kmers(seq, k=2)
        kmer_counts = dict(Counter(kmers))
        kmer_features.append(kmer_counts)

    vectorizer = DictVectorizer(sparse=False)
    X_kmers = vectorizer.fit_transform(kmer_features)

    X_other = df[["Length", "GC Content", "Molecular weight"]].to_numpy()
    X = np.hstack([X_kmers, X_other])

    y = df["Binding Affinity"].to_numpy()

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )

    model = RandomForestRegressor(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))
    print(f"RMSE: {rmse:.4f}")

    joblib.dump(model, "aptamer_model.pkl")
    joblib.dump(vectorizer, "dict_vectorizer.pkl")
    print("Model saved to aptamer_model.pkl")
    print("Vectorizer saved to dict_vectorizer.pkl")


if __name__ == "__main__":
    main()
