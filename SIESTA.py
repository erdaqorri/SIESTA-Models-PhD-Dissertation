import argparse
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier
import joblib

# Initialize argument parser
parser = argparse.ArgumentParser(description="SIESTA")
parser.add_argument("--features_csv", type=str, required=True, help="Path to the input features, in CSV format.")
parser.add_argument("--output_model", type=str, required=True, help="Path to save the trained SIESTA model.")
args = parser.parse_args()

def protein_aware_split(df, protein_col, test_size=0.2, random_state=None):
    """
    Split dataset into train and test sets while ensuring that
    variants of the same protein are in the same split.
    """
    if random_state is not None:
        np.random.seed(random_state)

    unique_proteins = df[protein_col].unique()
    n_test = int(len(unique_proteins) * test_size)
    test_proteins = np.random.choice(unique_proteins, size=n_test, replace=False)

    test_df = df[df[protein_col].isin(test_proteins)]
    train_df = df[~df[protein_col].isin(test_proteins)]

    # Print dataset statistics
    print(f"Total proteins: {len(unique_proteins)}")
    print(f"Train proteins: {len(train_df[protein_col].unique())}, samples: {len(train_df)}")
    print(f"Test proteins: {len(test_df[protein_col].unique())}, samples: {len(test_df)}")

    return train_df, test_df

# Load the features
try:
    features_df = pd.read_csv(args.features_csv)
except Exception as e:
    print(f"Error loading features CSV: {e}")
    exit(1)

# Perform the protein-aware split
train_df, test_df = protein_aware_split(
    features_df,
    protein_col='uniprotID',
    test_size=0.2,
    random_state=27
)

# Check unique proteins in train and test sets
train_proteins = train_df['uniprotID'].unique()
print("Train Proteins:", train_proteins)
test_proteins = test_df['uniprotID'].unique()
print("Test Proteins:", test_proteins)

# Specify feature columns
features = ['Cα-Dist', 'dRMS Local', 'ΔSASA Normalized', 'ΔCα-pLDDT', 'MJ Potential Mutant', 'Entropy', 'PSSM Native', 'Hydrophobicity', 'substitutionMatrix']

# Prepare training and testing data
x_train = train_df[features]
y_train = train_df['labels']

x_test = test_df[features]
y_test = test_df['labels']

# Standardize features
scaler = StandardScaler()
x_train = scaler.fit_transform(x_train)
x_test = scaler.transform(x_test)

# Initialize and fit the XGBoost classifier
xgb = XGBClassifier(
    objective="binary:logistic",
    eta=0.3,
    gamma=0,
    max_depth=6,
    min_child_weight=1,
    subsample=1,
    colsample_bytree=1,
    reg_lambda=1,
    alpha=0,
    scale_pos_weight=1,
    base_score=0.5,
    booster='gbtree',
    nthread=6,
    random_state=27
)

# Fit the model
xgb.fit(x_train, y_train.to_numpy().ravel())

# Save the trained model
joblib.dump(xgb, args.output_model)

print(f"SIESTA trained and saved to {args.output_model}")
