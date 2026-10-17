"""Oblique Stereographic Cartographic Projection (Stereo 70).

This module implements the Oblique Stereographic projection equations on the WGS84
reference ellipsoid as used by the Romanian national coordinate system (Stereo 70).

The transformation utilizes the Gauss double projection method:
1. Conformal mapping of the reference ellipsoid onto a conformal sphere.
2. Perspective stereographic projection from the conformal sphere to the projection plane.

Mathematical Formulation:
-------------------------
Ellipsoidal Parameters:
    - Semi-major axis: :math:`a = 6378137.0\\text{ m}`
    - Inverse flattening: :math:`1/f = 298.257223563`
    - First eccentricity squared: :math:`e^2 = 2f - f^2`
    - Meridian curvature radius: :math:`M = \\frac{a(1 - e^2)}{(1 - e^2 \\sin^2\\varphi)^{3/2}}`
    - Prime vertical curvature radius: :math:`N = \\frac{a}{\\sqrt{1 - e^2 \\sin^2\\varphi}}`

Conformal Sphere Parameters:
    - Sphere radius: :math:`R = \\sqrt{M_0 N_0}` at natural origin latitude :math:`\\varphi_0 = 46^\\circ\\text{N}`
    - Scale constant: :math:`n = \\sqrt{1 + \\frac{e^2 \\cos^4\\varphi_0}{1 - e^2}}`
    - Integration constant :math:`c` and conformal origin latitude :math:`\\chi_0`

Stereographic Forward Projection :math:`(\\varphi, \\lambda) \\to (N, E)`:
    Conformal sphere coordinates :math:`(\\chi, \\Lambda)`:
    .. math::
        \\chi = 2 \\arctan\\left[ c \\left( \\frac{1 + \\sin\\varphi}{1 - \\sin\\varphi} \\right)^{n/2}
        \\left( \\frac{1 - e\\sin\\varphi}{1 + e\\sin\\varphi} \\right)^{ne/2} \\right] - \\frac{\\pi}{2}

    .. math::
        \\Lambda = n(\\lambda - \\lambda_0) + \\lambda_0

    Projected Plane Coordinates:
    .. math::
        B = 1 + \\sin\\chi \\sin\\chi_0 + \\cos\\chi \\cos\\chi_0 \\cos(\\Lambda - \\lambda_0)

    .. math::
        N = N_0 + \\frac{2 R k_0 (\\sin\\chi \\cos\\chi_0 - \\cos\\chi \\sin\\chi_0 \\cos(\\Lambda - \\lambda_0))}{B}

    .. math::
        E = E_0 + \\frac{2 R k_0 \\cos\\chi \\sin(\\Lambda - \\lambda_0)}{B}

Stereographic Inverse Projection :math:`(N, E) \\to (\\varphi, \\lambda)`:
    Inverts planar coordinates to conformal sphere :math:`(\\chi, \\Lambda)`, derives isometric
    latitude :math:`q`, and iteratively resolves geodetic latitude :math:`\\varphi` using
    Newton-Raphson iteration until convergence threshold :math:`\\epsilon < 10^{-11}\\text{ rad}`.
"""

import math
from typing import Sequence, Tuple

from pytransdatro import constants, exceptions


class StereoProj:
    """Provides coordinate conversions for the Oblique Stereographic projection (Stereo 70).

    Supports:
        - Forward projection: :math:`(\\text{Lat}, \\text{Long}) \\to (N, E)`
        - Inverse projection: :math:`(N, E) \\to (\\text{Lat}, \\text{Long})`
    """

    # Reference ellipsoid constants (WGS84)
    WGS84_ELL_SMAJ_AXIS: float = constants.WGS84_ELL_SMAJ_AXIS
    WGS84_ELL_INV_FLAT: float = constants.WGS84_ELL_INV_FLAT

    # Projection parameters (Natural origin 46° N, 25° E)
    ORIG_LAT: float = constants.STEREO70_ORIG_LAT
    ORIG_LON: float = constants.STEREO70_ORIG_LON
    FALSE_N: float = constants.STEREO70_FALSE_N
    FALSE_E: float = constants.STEREO70_FALSE_E
    PROJ_SCALE: float = constants.STEREO70_SCALE_FACTOR

    def __init__(self) -> None:
        """Initialize projection engine with WGS84 ellipsoid and conformal sphere parameters."""
        self._ell: StereoProj._Ellipsoid = self._Ellipsoid(
            StereoProj.WGS84_ELL_SMAJ_AXIS,
            StereoProj.WGS84_ELL_INV_FLAT,
        )
        self._orig_lat: float = StereoProj.ORIG_LAT
        self._orig_lon: float = StereoProj.ORIG_LON
        self._f_n: float = StereoProj.FALSE_N
        self._f_e: float = StereoProj.FALSE_E
        self._s: float = StereoProj.PROJ_SCALE
        self.__conf_sphere: StereoProj._ConfSph = self._ConfSph(self)

    def to_grid(
        self,
        lat: Sequence[float],
        lon: Sequence[float],
    ) -> Tuple[list[float], list[float]]:
        """Converts geodetic (geographic) coordinates in radians to projected grid coordinates (meters).

        :math:`(\\text{Lat}, \\text{Long}) \\to (N, E)`

        :param lat: Sequence of geodetic latitudes in radians
        :type lat: Sequence[float]
        :param lon: Sequence of geodetic longitudes in radians
        :type lon: Sequence[float]
        :return: Tuple of ``(n_out, e_out)`` lists in meters
        :rtype: Tuple[list[float], list[float]]
        """
        lat_c_arr, lon_c_arr = self.__conf_sphere.get_latlon_from_geo(lat, lon)
        lat_c0 = self.__conf_sphere.orig_lat
        lon_c0 = self.__conf_sphere.orig_lon
        sin_lat_c0 = math.sin(lat_c0)
        cos_lat_c0 = math.cos(lat_c0)
        f_n = self._f_n
        f_e = self._f_e
        r2s = 2 * self.__conf_sphere.r * self._s
        _sin = math.sin
        _cos = math.cos

        n_out: list[float] = []
        e_out: list[float] = []
        for lat_c, lon_c in zip(lat_c_arr, lon_c_arr):
            sin_lat_c = _sin(lat_c)
            cos_lat_c = _cos(lat_c)
            dlon = lon_c - lon_c0
            b = 1 + sin_lat_c * sin_lat_c0 + cos_lat_c * cos_lat_c0 * _cos(dlon)
            n_out.append(
                f_n + r2s * (sin_lat_c * cos_lat_c0 - cos_lat_c * sin_lat_c0 * _cos(dlon)) / b
            )
            e_out.append(f_e + r2s * cos_lat_c * _sin(dlon) / b)
        return n_out, e_out

    def to_geo(
        self,
        n: Sequence[float],
        e: Sequence[float],
    ) -> Tuple[list[float], list[float]]:
        """Converts projected grid coordinates in meters to geodetic (geographic) coordinates in radians.

        :math:`(N, E) \\to (\\text{Lat}, \\text{Long})`

        :param n: Sequence of northings in meters
        :type n: Sequence[float]
        :param e: Sequence of eastings in meters
        :type e: Sequence[float]
        :return: Tuple of ``(lat_out, lon_out)`` lists in radians
        :rtype: Tuple[list[float], list[float]]
        """
        lat_c_arr, lon_c_arr = self.__conf_sphere.get_latlon_from_grid(n, e)

        e_first = self._ell.first_ecc
        e2 = self._ell.first_ecc2
        c_conf = self.__conf_sphere.c
        n_conf = self.__conf_sphere.n
        orig_lon = self._orig_lon

        _sin = math.sin
        _cos = math.cos
        _log = math.log
        _exp = math.exp
        _atan = math.atan
        _tan = math.tan
        _pow = math.pow
        pi_2 = math.pi / 2
        pi_4 = math.pi / 4
        tol = constants.STEREO70_CONV_TOLERANCE
        max_iter = constants.STEREO70_MAX_ITERATIONS

        lat_out: list[float] = []
        lon_out: list[float] = []
        for lat_c, lon_c in zip(lat_c_arr, lon_c_arr):
            sin_lat_c = _sin(lat_c)
            iso_lat = (0.5 * _log((1 + sin_lat_c) / (c_conf * (1 - sin_lat_c)))) / n_conf
            r_lat = 2 * _atan(_exp(iso_lat)) - pi_2
            diff = tol + 1
            i = 1
            while diff >= tol and i < max_iter:
                r_lat_before_next_iter = r_lat
                sin_rlat = _sin(r_lat)
                iso_lat_i = _log(
                    _tan(r_lat / 2 + pi_4)
                    * _pow((1 - e_first * sin_rlat) / (1 + e_first * sin_rlat), e_first / 2)
                )
                r_lat = r_lat - (
                    (iso_lat_i - iso_lat) * _cos(r_lat) * (1 - e2 * sin_rlat * sin_rlat) / (1 - e2)
                )
                diff = abs(r_lat_before_next_iter - r_lat)
                i += 1
            lat_out.append(r_lat)
            lon_out.append(orig_lon + (lon_c - orig_lon) / n_conf)
        return lat_out, lon_out

    class _Ellipsoid:
        """Internal helper calculating geometric elements of the reference ellipsoid."""

        def __init__(self, smaj_axis: float, inv_flat: float) -> None:
            """Initialize reference ellipsoid dimensions.

            :param smaj_axis: Semi-major axis :math:`a` in meters
            :type smaj_axis: float
            :param inv_flat: Inverse flattening :math:`1/f`
            :type inv_flat: float
            """
            self.smaj_axis: float = smaj_axis
            self.inv_flat: float = inv_flat
            self.flat: float = 1.0 / inv_flat
            self.smin_axis: float = smaj_axis * (1.0 - 1.0 / inv_flat)
            self.first_ecc2: float = 2.0 * self.flat - self.flat * self.flat
            self.first_ecc: float = math.sqrt(self.first_ecc2)
            self.sec_ecc2: float = self.first_ecc2 / (1.0 - self.first_ecc2)
            self.sec_ecc: float = math.sqrt(self.sec_ecc2)

        def get_rad_M_N(self, lat: float) -> Tuple[float, float]:
            """Calculates meridian curvature radius :math:`M` and prime vertical radius :math:`N`.

            .. math::
                M = \\frac{a(1 - e^2)}{(1 - e^2 \\sin^2\\varphi)^{3/2}}, \\quad
                N = \\frac{a}{\\sqrt{1 - e^2 \\sin^2\\varphi}}

            :param lat: Geodetic latitude in radians
            :type lat: float
            :return: Tuple of ``(M, N)`` curvature radii in meters
            :rtype: Tuple[float, float]
            """
            tmp_v = math.sqrt(1.0 - self.first_ecc2 * math.pow(math.sin(lat), 2))
            m = self.smaj_axis * (1.0 - self.first_ecc2) / math.pow(tmp_v, 3)
            n = self.smaj_axis / tmp_v
            return m, n

    class _ConfSph:
        """Internal helper defining the conformal sphere for the Gauss double projection."""

        def __init__(self, proj: "StereoProj") -> None:
            """Compute conformal sphere constants based on ellipsoid and natural origin latitude.

            :param proj: Parent StereoProj instance
            :type proj: StereoProj
            """
            r_m, r_n = proj._ell.get_rad_M_N(proj._orig_lat)
            e2 = proj._ell.first_ecc2
            e = proj._ell.first_ecc
            sin_orig_lat = math.sin(proj._orig_lat)
            s1 = (1 + sin_orig_lat) / (1 - sin_orig_lat)
            s2 = (1 - e * sin_orig_lat) / (1 + e * sin_orig_lat)
            r = math.sqrt(r_m * r_n)
            n = math.sqrt(1 + (e2 * math.pow(math.cos(proj._orig_lat), 4) / (1 - e2)))
            w1 = math.pow(s1, n) * math.pow(s2, e * n)
            sin_chi0 = (w1 - 1) / (w1 + 1)
            c = (n + sin_orig_lat) * (1 - sin_chi0) / ((n - sin_orig_lat) * (1 + sin_chi0))
            w2 = c * w1

            self.__proj: StereoProj = proj
            self.r: float = r
            self.n: float = n
            self.c: float = c
            self.orig_lat: float = math.asin((w2 - 1) / (w2 + 1))
            self.orig_lon: float = proj._orig_lon

        def get_latlon_from_grid(
            self,
            n: Sequence[float],
            e: Sequence[float],
        ) -> Tuple[list[float], list[float]]:
            """Computes equivalent conformal latitude and longitude on sphere from grid coordinates.

            :param n: Sequence of northing coordinates in meters
            :type n: Sequence[float]
            :param e: Sequence of easting coordinates in meters
            :type e: Sequence[float]
            :return: Tuple of ``(lat_out, lon_out)`` lists of conformal coordinates on sphere
            :rtype: Tuple[list[float], list[float]]
            """
            sf = self.__proj._s
            r = self.r
            g = 2 * r * sf * math.tan((math.pi / 4) - (self.orig_lat / 2))
            h = 4 * r * sf * math.tan(self.orig_lat) + g
            f_n = self.__proj._f_n
            f_e = self.__proj._f_e
            orig_lat = self.orig_lat
            orig_lon = self.orig_lon
            r2s = 2 * r * sf
            _atan = math.atan
            _tan = math.tan

            lat_out: list[float] = []
            lon_out: list[float] = []
            for n_val, e_val in zip(n, e):
                dn = n_val - f_n
                de = e_val - f_e
                i = _atan(de / (h + dn))
                j = _atan(de / (g - dn)) - i
                lat_out.append(orig_lat + 2 * _atan((dn - de * _tan(j / 2)) / r2s))
                lon_out.append(j + 2 * i + orig_lon)
            return lat_out, lon_out

        def get_latlon_from_geo(
            self,
            lat: Sequence[float],
            lon: Sequence[float],
        ) -> Tuple[list[float], list[float]]:
            """Computes equivalent conformal coordinates on sphere from ellipsoidal coordinates.

            :param lat: Sequence of geodetic latitudes in radians
            :type lat: Sequence[float]
            :param lon: Sequence of geodetic longitudes in radians
            :type lon: Sequence[float]
            :return: Tuple of ``(lat_out, lon_out)`` lists on conformal sphere
            :rtype: Tuple[list[float], list[float]]
            """
            e = self.__proj._ell.first_ecc
            n = self.n
            orig_lon = self.orig_lon
            c = self.c
            _sin = math.sin
            _asin = math.asin
            _pow = math.pow

            lat_out: list[float] = []
            lon_out: list[float] = []
            half_pi = math.pi / 2.0
            for lat_val, lon_val in zip(lat, lon):
                if abs(lat_val) >= half_pi or abs(lon_val) > math.pi:
                    raise exceptions.InvalidCoordinateErr(
                        "Geographic coordinates must satisfy |lat| < pi/2 and |lon| <= pi radians."
                    )
                sin_lat = _sin(lat_val)
                if abs(sin_lat) >= 1.0:
                    raise exceptions.InvalidCoordinateErr(
                        "Latitude cannot lie at the +/-90 degree pole singularity."
                    )
                sa = (1 + sin_lat) / (1 - sin_lat)
                sb = (1 - e * sin_lat) / (1 + e * sin_lat)
                w = c * _pow((sa * _pow(sb, e)), n)
                lat_out.append(_asin((w - 1) / (w + 1)))
                lon_out.append(n * (lon_val - orig_lon) + orig_lon)
            return lat_out, lon_out
