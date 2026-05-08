import pytest
import numpy as np
from QMzyme.utils import rmsd, compute_translation_and_rotation, kabsch_transform

def test_kabsch_and_rmsd_logic():
    # Setup independent coordinate sets
    target = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0]
    ], dtype=float)

    # Mobile has the coordinates rotated and shifted
    mobile = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [-1.0, 0.0, 0.0]
    ], dtype=float) + 10.0

    # Test compute_translation_and_rotation
    t, r = compute_translation_and_rotation(mobile, target)
    
    # Verify the rotation matrix is orthogonal (R * R_T = Identity)
    identity_matrix = np.eye(3)
    actual_identity = np.dot(r, r.T)
    
    # Compare the entire matrix at once
    assert actual_identity.flatten() == pytest.approx(identity_matrix.flatten(), abs=1e-6)

    # Test kabsch_transform
    transformed_mobile = kabsch_transform(mobile, t, r)
    
    # Compare flattened arrays for a clean pytest assertion
    assert transformed_mobile.flatten() == pytest.approx(target.flatten(), abs=1e-6)

    # Unaligned RMSD should be large
    raw_rmsd = rmsd(mobile, target, align=False)
    assert raw_rmsd > 10.0

    # Aligned RMSD should be zero
    aligned_rmsd = rmsd(mobile, target, align=True)
    assert aligned_rmsd == pytest.approx(0.0, abs=1e-6)

def test_rmsd_identity():
    coords = np.random.rand(5, 3)
    assert rmsd(coords, coords) == pytest.approx(0.0)

def test_rmsd_manual_delta():
    xyz1 = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
    xyz2 = np.array([[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]])
    
    # Each atom is sqrt(1^2 + 1^2 + 1^2) = sqrt(3) away
    # Mean of (sqrt(3)^2) is 3, sqrt(3) is ~1.732
    expected_rmsd = 1.73205081
    assert rmsd(xyz1, xyz2, align=False) == pytest.approx(expected_rmsd)