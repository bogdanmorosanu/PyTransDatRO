"""4-Parameter 2D Conformal Helmert Coordinate Transformation.

This module provides the Helmert2D class for 2D conformal Helmert transformations
between the Romanian national projection (Stereo 70) and the intermediate
conformal projection on the GRS80/WGS84 ellipsoid (StereoGRS80).

Transformation parameters are loaded dynamically from the active SPG binary grid
via SpgReader upon instantiation.
"""

import math
from typing import Optional, Sequence, Tuple

from pytransdatro.spg_reader import SpgReader


class Helmert2D:
    """Calculates the 2D conformal Helmert transformation between Stereo70 and StereoGRS80.

    Mathematical formulation:
        For forward transformation (:math:`\\text{sign} = +1`):
        .. math::
            \\begin{pmatrix} N' \\\\ E' \\end{pmatrix} =
            (1 + \\delta m)
            \\begin{pmatrix}
                \\cos\\theta & \\sin\\theta \\\\
                -\\sin\\theta & \\cos\\theta
            \\end{pmatrix}
            \\begin{pmatrix} N \\\\ E \\end{pmatrix} +
            \\begin{pmatrix} t_N \\\\ t_E \\end{pmatrix}

        For inverse transformation (:math:`\\text{sign} = -1`):
        .. math::
            \\begin{pmatrix} N' \\\\ E' \\end{pmatrix} =
            (1 - \\delta m)
            \\begin{pmatrix}
                \\cos\\theta & -\\sin\\theta \\\\
                \\sin\\theta & \\cos\\theta
            \\end{pmatrix}
            \\begin{pmatrix} N \\\\ E \\end{pmatrix} -
            \\begin{pmatrix} t_N \\\\ t_E \\end{pmatrix}

    where:
        - :math:`t_N, t_E`: Translations along Northing and Easting axes (meters)
        - :math:`\\theta`: Rotation angle (radians)
        - :math:`\\delta m = \\text{ppm} \\cdot 10^{-6}`: Dimensionless scale distortion factor
    """

    def __init__(
        self,
        tn: Optional[float] = None,
        te: Optional[float] = None,
        ppm: Optional[float] = None,
        rot_rad: Optional[float] = None,
    ) -> None:
        """Initialize Helmert2D with explicit parameters or values from active SPG grid.

        :param tn: Translation in Northing (meters), defaults to SPG grid parameter
        :type tn: Optional[float]
        :param te: Translation in Easting (meters), defaults to SPG grid parameter
        :type te: Optional[float]
        :param ppm: Scale distortion in parts-per-million, defaults to SPG grid parameter
        :type ppm: Optional[float]
        :param rot_rad: Rotation in radians, defaults to SPG grid parameter
        :type rot_rad: Optional[float]
        """
        reader = SpgReader()
        self.__tn: float = reader.helmert_tn if tn is None else float(tn)
        self.__te: float = reader.helmert_te if te is None else float(te)
        self.__ppm: float = reader.helmert_ppm if ppm is None else float(ppm)
        self.__r: float = reader.helmert_rot_rad if rot_rad is None else float(rot_rad)

        # Precompute constants to optimize speed for batch calculations
        self.__sinr: float = math.sin(self.__r)
        self.__cosr: float = math.cos(self.__r)
        self.__scale_factor: float = self.__ppm * 1e-6

    @property
    def tn(self) -> float:
        """Northing translation parameter in meters."""
        return self.__tn

    @property
    def te(self) -> float:
        """Easting translation parameter in meters."""
        return self.__te

    @property
    def ppm(self) -> float:
        """Scale distortion parameter in parts-per-million."""
        return self.__ppm

    @property
    def r(self) -> float:
        """Rotation angle parameter in radians."""
        return self.__r

    def trans(
        self,
        n: Sequence[float],
        e: Sequence[float],
        sign: int,
    ) -> Tuple[list[float], list[float]]:
        """Transforms sequences of northing and easting coordinates: :math:`(N, E) \\to (N', E')`.

        :param n: Sequence of northing coordinates (meters)
        :type n: Sequence[float]
        :param e: Sequence of easting coordinates (meters)
        :type e: Sequence[float]
        :param sign: Transformation direction: ``+1`` for forward (Stereo70 to StereoGRS80)
            or ``-1`` for inverse (StereoGRS80 to Stereo70)
        :type sign: int
        :return: Tuple of ``(n_out, e_out)`` lists of transformed coordinates
        :rtype: Tuple[list[float], list[float]]
        """
        scale = 1.0 + sign * self.__scale_factor
        cosr_scale = scale * self.__cosr
        sinr_scale_sign = scale * sign * self.__sinr
        tn_sign = sign * self.__tn
        te_sign = sign * self.__te

        n_out: list[float] = []
        e_out: list[float] = []
        for n_val, e_val in zip(n, e):
            n_out.append(n_val * cosr_scale + e_val * sinr_scale_sign + tn_sign)
            e_out.append(e_val * cosr_scale - n_val * sinr_scale_sign + te_sign)
        return n_out, e_out
