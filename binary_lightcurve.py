"""
Simple eclipsing binary light curve simulator.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import List, Sequence, Tuple

SIGMA_SB = 5.670374419e-8  # Stefan-Boltzmann constant (W/m^2/K^4)
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)


@dataclass
class Star:
    mass: float  # kilograms
    radius: float  # meters
    temperature: float  # Kelvin

    @property
    def luminosity(self) -> float:
        """Return bolometric luminosity using Stefan-Boltzmann law."""
        return 4 * math.pi * self.radius ** 2 * SIGMA_SB * self.temperature ** 4

    @property
    def surface_brightness(self) -> float:
        return SIGMA_SB * self.temperature ** 4


@dataclass
class BinarySystem:
    primary: Star
    secondary: Star
    semi_major_axis: float  # meters
    eccentricity: float
    inclination: float  # radians
    period: float  # seconds

    def _mean_to_true_anomaly(self, mean_anomaly: float) -> float:
        """Solve Kepler's equation for the true anomaly."""
        # Newton-Raphson on eccentric anomaly
        ecc_anomaly = mean_anomaly
        for _ in range(50):
            f = ecc_anomaly - self.eccentricity * math.sin(ecc_anomaly) - mean_anomaly
            fp = 1 - self.eccentricity * math.cos(ecc_anomaly)
            delta = -f / fp
            ecc_anomaly += delta
            if abs(delta) < 1e-10:
                break
        cos_nu = (math.cos(ecc_anomaly) - self.eccentricity) / (
            1 - self.eccentricity * math.cos(ecc_anomaly)
        )
        sin_nu = (math.sqrt(1 - self.eccentricity ** 2) * math.sin(ecc_anomaly)) / (
            1 - self.eccentricity * math.cos(ecc_anomaly)
        )
        return math.atan2(sin_nu, cos_nu)

    def _orbital_separation(self, true_anomaly: float) -> float:
        return self.semi_major_axis * (1 - self.eccentricity ** 2) / (
            1 + self.eccentricity * math.cos(true_anomaly)
        )

    def _projected_positions(self, true_anomaly: float) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
        r = self._orbital_separation(true_anomaly)
        # Position of secondary relative to primary along orbital plane
        x_orb = r * math.cos(true_anomaly)
        y_orb = r * math.sin(true_anomaly)

        # Center-of-mass scaling
        m1, m2 = self.primary.mass, self.secondary.mass
        x2, y2 = (m1 / (m1 + m2)) * x_orb, (m1 / (m1 + m2)) * y_orb
        x1, y1 = -(m2 / (m1 + m2)) * x_orb, -(m2 / (m1 + m2)) * y_orb

        # Rotate by inclination about x-axis
        cos_i, sin_i = math.cos(self.inclination), math.sin(self.inclination)
        y1_proj, z1_proj = y1 * cos_i, y1 * sin_i
        y2_proj, z2_proj = y2 * cos_i, y2 * sin_i

        star1 = (x1, y1_proj, z1_proj)
        star2 = (x2, y2_proj, z2_proj)
        return star1, star2

    def _circle_overlap(self, d: float, r1: float, r2: float) -> float:
        """Return area of overlap of two circles with separation d."""
        if d >= r1 + r2:
            return 0.0
        if d <= abs(r1 - r2):
            return math.pi * min(r1, r2) ** 2

        r1_sq, r2_sq = r1 ** 2, r2 ** 2
        alpha = math.acos((d ** 2 + r1_sq - r2_sq) / (2 * d * r1))
        beta = math.acos((d ** 2 + r2_sq - r1_sq) / (2 * d * r2))
        area = r1_sq * alpha + r2_sq * beta - d * r1 * math.sin(alpha)
        return area

    def _flux_with_occultation(
        self, projected_separation: float, front: Star, back: Star, front_radius: float, back_radius: float
    ) -> float:
        overlap = self._circle_overlap(projected_separation, front_radius, back_radius)
        # area of back visible after being partially blocked
        back_area = math.pi * back_radius ** 2 - overlap
        return front.surface_brightness * math.pi * front_radius ** 2 + back.surface_brightness * back_area

    def light_curve(self, phase_count: int = 500) -> Tuple[List[float], List[float]]:
        phases = _linspace(0.0, 1.0, phase_count)
        fluxes: List[float] = []

        total_area_flux = (
            self.primary.surface_brightness * math.pi * self.primary.radius ** 2
            + self.secondary.surface_brightness * math.pi * self.secondary.radius ** 2
        )

        for phase in phases:
            mean_anomaly = 2 * math.pi * phase
            true_anomaly = self._mean_to_true_anomaly(mean_anomaly)
            star1, star2 = self._projected_positions(true_anomaly)
            separation = _euclidean_distance_2d(star1[:2], star2[:2])

            # Determine which star is in front (smaller z means closer to observer)
            if star1[2] < star2[2]:
                flux = self._flux_with_occultation(separation, self.primary, self.secondary, self.primary.radius, self.secondary.radius)
            else:
                flux = self._flux_with_occultation(separation, self.secondary, self.primary, self.secondary.radius, self.primary.radius)

            fluxes.append(flux / total_area_flux)

        return phases, fluxes


def period_from_parameters(semi_major_axis: float, mass1: float, mass2: float) -> float:
    return 2 * math.pi * math.sqrt(semi_major_axis ** 3 / (G * (mass1 + mass2)))


def _linspace(start: float, stop: float, num: int) -> List[float]:
    if num <= 1:
        return [start]
    step = (stop - start) / (num - 1)
    return [start + i * step for i in range(num)]


def _euclidean_distance_2d(a: Sequence[float], b: Sequence[float]) -> float:
    return math.hypot(a[0] - b[0], a[1] - b[1])
