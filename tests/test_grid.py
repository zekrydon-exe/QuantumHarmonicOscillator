import numpy as np
import pytest

from src.hamiltonian import create_grid


def test_grid_has_correct_number_of_points():
    q, dq = create_grid(-4.0, 4.0, 9)

    assert len(q) == 9


def test_grid_has_correct_spacing():
    q, dq = create_grid(-4.0, 4.0, 9)

    assert np.isclose(dq, 0.8)


def test_grid_excludes_boundaries():
    q, dq = create_grid(-4.0, 4.0, 9)

    assert q[0] > -4.0
    assert q[-1] < 4.0


def test_grid_is_uniformly_spaced():
    q, dq = create_grid(-4.0, 4.0, 9)

    differences = np.diff(q)

    assert np.allclose(differences, dq)


def test_symmetric_grid_includes_origin():
    q, dq = create_grid(-4.0, 4.0, 9)

    assert np.isclose(q[4], 0.0)
    assert np.allclose(q, -q[::-1])

def test_grid_rejects_non_integer_point_count():
    with pytest.raises(TypeError):
        create_grid(-4.0, 4.0, 9.5)


def test_grid_rejects_boolean_point_count():
    with pytest.raises(TypeError):
        create_grid(-4.0, 4.0, True)


def test_grid_rejects_nonpositive_point_count():
    with pytest.raises(ValueError):
        create_grid(-4.0, 4.0, 0)


def test_grid_rejects_reversed_boundaries():
    with pytest.raises(ValueError):
        create_grid(4.0, -4.0, 9)


def test_grid_rejects_equal_boundaries():
    with pytest.raises(ValueError):
        create_grid(4.0, 4.0, 9)


def test_grid_rejects_nonfinite_boundaries():
    invalid_boundaries = [
        (float("nan"), 4.0),
        (-4.0, float("inf")),
        (float("-inf"), 4.0),
    ]

    for q_min, q_max in invalid_boundaries:
        with pytest.raises(ValueError):
            create_grid(q_min, q_max, 9)


def test_grid_rejects_invalid_boundary_types():
    with pytest.raises(TypeError):
        create_grid("left", 4.0, 9)

    with pytest.raises(TypeError):
        create_grid(-4.0, None, 9)


def test_grid_rejects_boolean_boundaries():
    with pytest.raises(TypeError):
        create_grid(True, 4.0, 9)