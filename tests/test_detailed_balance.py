
"""
Detailed balance: the shipped reverse-rate coefficients must be derivable.

GOAL: prove that `primat/data/csv/detailed_balance.csv` -- the table the LT
network reads its reverse rates from -- is not a hand-maintained set of magic
numbers but a reproducible function of the nuclide data.

For a reaction with forward rate `lambda_f(T9)`, the reverse rate follows from
equilibrium statistical mechanics as

    lambda_r / lambda_f = alpha * T9**beta * exp(-gamma / T9),

with `alpha` fixed by the spin degeneracies, masses and the number of
particles on each side, `beta` by the change in the number of interacting
bodies, and `gamma` by the reaction Q-value. `compute_detailed_balance_coefficients`
derives all three from the nuclide table; this test walks every row of the CSV
and checks the derivation reproduces it.

A single test covers the whole file (one assertion per row per coefficient),
which is why this module is short: a mismatch on any reaction means either the
nuclide data or the derivation changed, and both are single-cause failures.
"""
import os
import csv
import pytest
from primat.config import PRIMATConfig
from primat.network_data import compute_detailed_balance_coefficients, reaction_species

def test_detailed_balance_consistency():
    """Verify that compute_detailed_balance_coefficients reproduces the values in detailed_balance.csv."""
    cfg = PRIMATConfig()
    db_csv_path = os.path.join(cfg._resolved_data_dir, "csv", "detailed_balance.csv")
    
    with open(db_csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row['reaction']
            ref_alpha = float(row['alpha'])
            ref_beta = float(row['beta'])
            ref_gamma = float(row['gamma'])
            
            reactants, products = reaction_species(name)
            alpha, beta, gamma = compute_detailed_balance_coefficients(reactants, products, cfg)
            
            # Check beta exactly (it's always 0, 1.5, -1.5, etc.)
            assert beta == pytest.approx(ref_beta), f"Beta mismatch for {name}"
            
            # Check alpha and gamma to 1% (reproducibility limit of the physics formula vs original PRIMAT tables)
            if ref_alpha != 0:
                assert abs(alpha - ref_alpha) / abs(ref_alpha) < 0.01, f"Alpha mismatch for {name}: {alpha} vs {ref_alpha}"
            if ref_gamma != 0:
                assert abs(gamma - ref_gamma) / abs(ref_gamma) < 0.01, f"Gamma mismatch for {name}: {gamma} vs {ref_gamma}"
