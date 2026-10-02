"""Local Joukowsky geometry utilities for the rock-mask descriptor."""

import numpy as np


def polygon_area(points):
    x = points[:, 0]
    y = points[:, 1]
    return float(0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def polygon_perimeter(points):
    return float(np.linalg.norm(np.roll(points, -1, axis=0) - points, axis=1).sum())


def anisotropy(points, eps=1e-12):
    centered = points - points.mean(axis=0)
    vals = np.linalg.eigvalsh(centered.T @ centered / max(len(points), 1))
    return float((vals[-1] + eps) / (vals[0] + eps))


def joukowsky_metrics_v2(points_2d, a=0.25):
    """Rotation/reflection/order-invariant Joukowsky metrics in physical units."""
    points = np.asarray(points_2d, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or len(points) < 3:
        raise ValueError("O polígono precisa ter pelo menos três pontos 2D.")
    if not np.isfinite(points).all() or not np.isfinite(a) or a < 0:
        raise ValueError("Coordenadas e a devem ser finitos; a >= 0.")
    centered = points - points.mean(axis=0)
    z = centered[:, 0] + 1j * centered[:, 1]
    radius = float(np.mean(np.abs(z)))
    if radius <= 1e-10:
        raise ValueError("Polígono degenerado.")
    normalized = z / radius
    edges = np.roll(z, -1) - z
    if np.any(np.abs(edges) <= 1e-10):
        raise ValueError("Aresta degenerada.")
    measures = []
    for edge in edges:
        oriented = normalized * np.conj(edge / abs(edge))
        transformed = radius * (oriented + a * a / oriented)
        coords = np.column_stack((transformed.real, transformed.imag))
        measures.append((polygon_area(coords), polygon_perimeter(coords), anisotropy(coords)))
    mean = np.mean(measures, axis=0)
    return {
        "mean_radius": radius,
        "j_area": float(mean[0]),
        "j_perimeter": float(mean[1]),
        "j_anisotropy": float(mean[2]),
    }