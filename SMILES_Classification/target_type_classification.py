"""
target_type_classification.py

Classifies small molecules based on their SMILES notation into defined
functional categories using a machine learning classifier (RandomForest).

This script illustrates a complete SMILES-based classification workflow:

1. Convert SMILES strings to RDKit molecular fingerprints
2. Create feature vectors usable in machine learning
3. Train/test split
4. Train RandomForest classifier
5. Evaluate using accuracy score

Note:
    - Dataset must contain columns: "SMILES" and "Target_Type"
    - This script is modular and can be used for other molecular classification tasks
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score


def smiles_to_fingerprint(smiles, radius=2, n_bits=2048):
    """
    Convert a SMILES string into a molecular fingerprint.
    Returns a NumPy array of bits.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    return list(fp)


def main():
    filepath = "/content/smiles_dataset.xlsx"  # Update before use

    df = pd.read_excel(filepath)
    df = df.dropna(subset=["SMILES", "Target_Type"])

    # Convert SMILES to fingerprints
    fingerprints = []
    valid_labels = []
    for smiles, label in zip(df["SMILES"], df["Target_Type"]):
        fp = smiles_to_fingerprint(smiles)
        if fp is not None:
            fingerprints.append(fp)
            valid_labels.append(label)

    X = fingerprints
    y = valid_labels

    # Train/test split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )

    # Train classifier
    model = RandomForestClassifier(n_estimators=200, random_state=42)
    model.fit(X_train, y_train)

    # Test performance
    y_pred = model.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    print(f"Classification Accuracy: {accuracy:.4f}")


if __name__ == "__main__":
    main()
